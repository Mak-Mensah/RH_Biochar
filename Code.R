############################################################
## GOLD-STANDARD META-ANALYSIS SCRIPT
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(metafor)
  library(clubSandwich)
  library(janitor)
  library(readr)
  library(fs)
  library(glue)
  library(flextable)
  library(officer)
  library(stringr)
  library(patchwork)
  library(ggplot2)
  library(grid)
})

## =========================================================
## 1) PATHS
## =========================================================
input_csv  <- "CWPSdata.csv"
## If you want to hard-point to uploaded path instead:
## input_csv <- "/mnt/data/CWPSdata.csv"

output_dir <- "meta_outputs_cwps"
dir_create(output_dir)

if (!file.exists(input_csv)) {
  stop("Cannot find '", input_csv, "'. Put it in your working directory or set input_csv to a full path.")
}

## =========================================================
## 2) HELPERS
## =========================================================
pct <- function(x) (exp(x) - 1) * 100

sig_code <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}

fmt_ci <- function(est, se) {
  lo <- est - 1.96 * se
  hi <- est + 1.96 * se
  glue("{round(est,3)} [{round(lo,3)}, {round(hi,3)}]")
}

fmt_pct_ci <- function(est_lnrr, se_lnrr) {
  lo <- est_lnrr - 1.96 * se_lnrr
  hi <- est_lnrr + 1.96 * se_lnrr
  glue("{sprintf('%+.1f', pct(est_lnrr))}% [{sprintf('%+.1f', pct(lo))}, {sprintf('%+.1f', pct(hi))}]")
}

z_na <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

safe_file <- function(x) {
  x <- as.character(x)
  x <- gsub("[/\\\\]", "-", x)
  x <- gsub("[:*?\"<>|]", "", x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("_+", "_", x)
  trimws(x)
}

theme_pub <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = rel(0.95)),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

## robust numeric parsing (handles "15,000", etc.)
parse_num <- function(x) readr::parse_number(as.character(x), locale = locale(grouping_mark = ","))

fit_mv <- function(df, mods = NULL) {
  if (is.null(mods)) {
    rma.mv(yi, vi,
           random = ~ 1 | study_id/comparison_id,
           data = df, method = "REML")
  } else {
    rma.mv(yi, vi,
           mods = mods,
           random = ~ 1 | study_id/comparison_id,
           data = df, method = "REML")
  }
}

cr2_test <- function(m, cluster) {
  clubSandwich::coef_test(m, cluster = cluster, vcov = "CR2")
}

i2_mv_approx <- function(m) {
  tau2 <- sum(m$sigma2, na.rm = TRUE)
  100 * tau2 / (tau2 + mean(m$vi, na.rm = TRUE))
}

tau2_total <- function(m) sum(m$sigma2, na.rm = TRUE)

R2_tau <- function(m_null, m_mod) {
  t0 <- tau2_total(m_null)
  t1 <- tau2_total(m_mod)
  if (is.na(t0) || t0 <= 0) return(NA_real_)
  pmax(0, (t0 - t1) / t0)
}

save_docx_table <- function(df, title, path) {
  ft <- flextable(df) |>
    autofit() |>
    align(align = "center", part = "all") |>
    add_header_lines(values = title)
  save_as_docx(setNames(list(ft), title), path = path)
}

plot_residual_funnel <- function(resid, sei, title, subtitle = NULL) {
  d <- tibble(resid = resid, sei = sei) |> filter(is.finite(resid), is.finite(sei))
  if (nrow(d) < 5) return(NULL)
  
  ggplot(d, aes(x = resid, y = sei)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.6) +
    scale_y_reverse() +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Residuals (yi - fitted)",
      y = "SE (sqrt(vi), reversed)"
    ) +
    theme_pub(13)
}

pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf

save_plot_robust <- function(p, filename, width = 9, height = 6, dpi = 450, bg = "white") {
  out <- file.path(output_dir, filename)
  
  ok <- tryCatch({
    ggsave(out, plot = p, width = width, height = height, dpi = dpi, bg = bg)
    TRUE
  }, error = function(e) {
    message("ggsave failed for ", filename, ": ", conditionMessage(e))
    FALSE
  })
  
  if (!ok) {
    tryCatch({
      png(out, width = width * dpi, height = height * dpi, res = dpi)
      print(p)
      dev.off()
      message("Saved via png() fallback: ", filename)
    }, error = function(e) {
      message("png() fallback failed for ", filename, ": ", conditionMessage(e))
      try(dev.off(), silent = TRUE)
      stop("Failed to save plot: ", filename)
    })
  }
}

fit_mv_stats <- function(d, fml) {
  m0 <- fit_mv(d)
  m1 <- fit_mv(d, mods = fml)
  rob <- cr2_test(m1, d$study_id)
  
  aic0 <- tryCatch(AIC(m0), error = function(e) NA_real_)
  aic1 <- tryCatch(AIC(m1), error = function(e) NA_real_)
  
  list(
    m0 = m0,
    m1 = m1,
    rob = rob,
    AIC_null = aic0,
    AIC_mod  = aic1,
    Delta_AIC = aic1 - aic0,
    R2_tau = R2_tau(m0, m1),
    tau2_null = tau2_total(m0),
    tau2_mod  = tau2_total(m1)
  )
}

ensure_complete <- function(df, mods, min_k = 10, min_studies = 5) {
  mods <- mods[mods %in% names(df)]
  if (length(mods) == 0) stop("No moderators supplied (none found in df).")
  
  while (length(mods) > 0) {
    d <- df |> dplyr::filter(dplyr::if_all(dplyr::all_of(mods), ~ !is.na(.x)))
    
    if (nrow(d) >= min_k && dplyr::n_distinct(d$study_id) >= min_studies) {
      return(list(mods = mods, data = d))
    }
    
    na_rate <- vapply(mods, function(m) mean(is.na(df[[m]])), numeric(1))
    if (all(is.na(na_rate)) || length(na_rate) == 0) break
    
    drop_one <- names(which.max(na_rate))[1]
    miss_pct <- round(100 * max(na_rate, na.rm = TRUE), 1)
    
    message(sprintf(
      "Dropping '%s' (missingness=%.1f%%) to retain enough complete cases.",
      drop_one, miss_pct
    ))
    
    mods <- setdiff(mods, drop_one)
  }
  
  stop("No set of moderators yields enough complete cases. Use fewer moderators or fit separate models.")
}

drop_high_cor <- function(df, vars, cutoff = 0.70) {
  vars <- vars[vars %in% names(df)]
  if (length(vars) < 2) return(vars)
  M <- df |> select(all_of(vars)) |> as.matrix()
  cm <- suppressWarnings(cor(M, use = "pairwise.complete.obs"))
  
  drop <- character(0)
  while (TRUE) {
    cm2 <- cm; diag(cm2) <- 0
    if (all(is.na(cm2)) || max(abs(cm2), na.rm = TRUE) < cutoff) break
    idx <- which(abs(cm2) == max(abs(cm2), na.rm = TRUE), arr.ind = TRUE)[1, ]
    a <- colnames(cm2)[idx[1]]; b <- colnames(cm2)[idx[2]]
    na_a <- sum(is.na(df[[a]])); na_b <- sum(is.na(df[[b]]))
    drop_one <- if (na_a > na_b) a else if (na_b > na_a) b else b
    drop <- union(drop, drop_one)
    keep <- setdiff(colnames(cm), drop)
    if (length(keep) < 2) break
    cm <- cm[keep, keep, drop = FALSE]
  }
  setdiff(vars, drop)
}

## collapse sparse/high-card factors to Top-N + Other (MV-only)
collapse_topN <- function(df, var, topN = 10, other = "Other", missing = "Unknown") {
  if (!var %in% names(df)) return(df)
  x <- as.character(df[[var]])
  x[is.na(x) | x == ""] <- missing
  tab <- sort(table(x), decreasing = TRUE)
  keep <- names(tab)[seq_len(min(topN, length(tab)))]
  x[!x %in% keep] <- other
  df[[var]] <- factor(x)
  df
}

## robust MV fit with optimizer cascade (+ simplified random as last resort)
fit_mv_safe <- function(df, mods, random = ~ 1 | study_id/comparison_id, method = "REML") {
  
  tries <- list(
    list(control = list(optimizer = "nlminb", eval.max = 20000, iter.max = 20000)),
    list(control = list(optimizer = "optim", optmethod = "BFGS", maxit = 200000)),
    list(control = list(optimizer = "optim", optmethod = "Nelder-Mead", maxit = 200000))
  )
  
  for (i in seq_along(tries)) {
    m <- tryCatch(
      metafor::rma.mv(
        yi, vi,
        mods   = mods,
        random = random,
        data   = df,
        method = method,
        control = tries[[i]]$control
      ),
      error = function(e) NULL
    )
    if (!is.null(m)) return(m)
  }
  
  ## last resort: simplify random structure
  for (i in seq_along(tries)) {
    m <- tryCatch(
      metafor::rma.mv(
        yi, vi,
        mods   = mods,
        random = ~ 1 | study_id,
        data   = df,
        method = method,
        control = tries[[i]]$control
      ),
      error = function(e) NULL
    )
    if (!is.null(m)) return(m)
  }
  
  NULL
}

read_plot_rds <- function(path) {
  if (!file.exists(path)) stop("Missing panel RDS: ", path)
  readRDS(path)
}

vif_manual <- function(X) {
  X <- as.matrix(X)
  v <- rep(NA_real_, ncol(X))
  names(v) <- colnames(X)
  for (j in seq_len(ncol(X))) {
    yj <- X[, j]
    xj <- X[, -j, drop = FALSE]
    if (sd(yj, na.rm = TRUE) == 0) { v[j] <- NA_real_; next }
    fit <- lm(yj ~ xj)
    r2 <- summary(fit)$r.squared
    v[j] <- 1 / (1 - r2)
  }
  tibble(Term = names(v), VIF = as.numeric(v)) |> arrange(desc(VIF))
}

## =========================================================
## 3) READ + CLEAN (CWPS)
## =========================================================
dat_raw <- read_csv(input_csv, locale = locale(encoding = "latin1"),
                    show_col_types = FALSE) |>
  clean_names() |>
  mutate(across(where(is.character), ~na_if(.x, "")))

## Harmonize names for moderators you listed
if ("soil_ec_d_s_m" %in% names(dat_raw) && !"soil_ec_ds_m" %in% names(dat_raw)) {
  dat_raw <- dat_raw |> rename(soil_ec_ds_m = soil_ec_d_s_m)
}
if ("soil_ph" %in% names(dat_raw) && !"soil_p_h" %in% names(dat_raw)) {
  dat_raw <- dat_raw |> rename(soil_p_h = soil_ph)
}
if ("biochar_ph" %in% names(dat_raw) && !"biochar_p_h" %in% names(dat_raw)) {
  dat_raw <- dat_raw |> rename(biochar_p_h = biochar_ph)
}
if ("climate_zone_koppen" %in% names(dat_raw) && !"climate_zone" %in% names(dat_raw)) {
  dat_raw <- dat_raw |> rename(climate_zone = climate_zone_koppen)
}

## Required columns for lnRR (ROM)
req <- c("study_id",
         "outcome_variable",
         "treatment_mean","treatment_sd","treatment_n",
         "control_mean","control_sd","control_n")
miss <- setdiff(req, names(dat_raw))
if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

## Numeric parsing for outcomes + numeric moderators
num_cols <- intersect(
  c("soil_organic_matter_percent","soil_ec_ds_m","soil_p_h",
    "biochar_rate_t_ha","biochar_p_h","biochar_cec","biochar_toc_percent",
    "treatment_mean","treatment_sd","treatment_n",
    "control_mean","control_sd","control_n"),
  names(dat_raw)
)

dat <- dat_raw |>
  mutate(across(all_of(num_cols), parse_num)) |>
  mutate(
    study_id = factor(study_id),
    
    outcome_type = factor(ifelse(is.na(outcome_variable) | outcome_variable == "",
                                 "Unknown", outcome_variable)),
    
    group = case_when(
      str_detect(as.character(outcome_type), regex("yield|biomass|establishment", TRUE)) ~ "Plant performance",
      str_detect(as.character(outcome_type), regex("water use efficiency", TRUE)) ~ "Water relations",
      str_detect(as.character(outcome_type), regex("physiological|nutrient", TRUE)) ~ "Physiology / nutrients",
      str_detect(as.character(outcome_type), regex("root", TRUE)) ~ "Root traits",
      str_detect(as.character(outcome_type), regex("soil", TRUE)) ~ "Soil properties",
      TRUE ~ "Other"
    ) |> factor(),
    
    ## study_type standard (from experiment_type)
    study_type = if ("experiment_type" %in% names(dat_raw)) factor(experiment_type) else factor("Unknown"),
    
    comparison_id = row_number()
  ) |>
  filter(
    treatment_mean > 0, control_mean > 0,
    treatment_n >= 2, control_n >= 2,
    treatment_sd > 0, control_sd > 0
  )

if ("Field" %in% levels(dat$study_type)) {
  dat$study_type <- relevel(dat$study_type, ref = "Field")
}

write_csv(dat, file.path(output_dir, "data_after_cleaning.csv"))

## =========================================================
## 4) EFFECT SIZES (lnRR via ROM) + sei + z-scores
## =========================================================
dat_es <- escalc(
  measure = "ROM",
  m1i = treatment_mean, sd1i = treatment_sd, n1i = treatment_n,
  m2i = control_mean,   sd2i = control_sd,   n2i = control_n,
  data = dat
) |>
  as_tibble() |>
  mutate(sei = sqrt(vi))

## ===================== EDIT APPLIED HERE =====================
## Ensure key categorical moderators are factors so factor plots export.
key_fac_mods <- c("climate_zone", "drought_level", "drought_timing", "soil_texture")
for (v in intersect(key_fac_mods, names(dat_es))) {
  dat_es[[v]] <- dat_es[[v]] |>
    as.character() |>
    na_if("") |>
    replace_na("Unknown") |>
    factor()
}
## =============================================================

exclude_from_z <- c(
  "yi","vi","sei","zi","pval","ci_lb","ci_ub",
  "treatment_mean","treatment_sd","treatment_n",
  "control_mean","control_sd","control_n"
)
num_mods <- names(dat_es)[sapply(dat_es, is.numeric)]
num_mods <- setdiff(num_mods, exclude_from_z)
for (nm in num_mods) dat_es[[paste0(nm, "_z")]] <- z_na(dat_es[[nm]])

write_csv(dat_es, file.path(output_dir, "effect_sizes_lnRR.csv"))

## =========================================================
## 4B) MODERATOR FIGURES (NUMERIC bubbles + FACTOR coefficient plots)
## =========================================================
min_k_total <- 10
min_studies <- 5
min_k_per_level <- 3
max_levels_plot <- 15

exclude_vars <- c(
  "study_id","comparison_id",
  "outcome_variable","outcome_unit","outcome_type","group",
  "treatment_mean","treatment_sd","treatment_n",
  "control_mean","control_sd","control_n",
  "yi","vi","sei"
)
cand <- setdiff(names(dat_es), exclude_vars)

## ---- NUMERIC bubble plots ----
num_cand <- cand[vapply(dat_es[cand], is.numeric, logical(1))]
num_cand <- setdiff(num_cand, "vi")

for (xvar in num_cand) {
  
  d <- dat_es |> filter(is.finite(yi), is.finite(vi), is.finite(.data[[xvar]]))
  if (nrow(d) < min_k_total || n_distinct(d$study_id) < min_studies) next
  
  fml <- as.formula(paste0("~ ", xvar))
  stats <- tryCatch(fit_mv_stats(d, fml), error = function(e) NULL)
  if (is.null(stats)) next
  rob <- stats$rob
  
  idx <- match(xvar, rownames(rob))
  if (is.na(idx)) next
  
  slope_est <- rob$beta[idx]
  slope_se  <- rob$SE[idx]
  slope_p   <- rob$p_Satt[idx]
  
  lab <- glue(
    "CR2 slope: β = {round(slope_est, 4)} (SE = {round(slope_se, 4)}), p = {signif(slope_p, 3)}\n",
    "R²(τ²) = {ifelse(is.na(stats$R2_tau), 'NA', sprintf('%.3f', stats$R2_tau))}; ΔAIC = {ifelse(is.na(stats$Delta_AIC), 'NA', sprintf('%.2f', stats$Delta_AIC))}\n",
    "k = {nrow(d)}, studies = {n_distinct(d$study_id)}"
  )
  
  p <- ggplot(d, aes(x = .data[[xvar]], y = yi, size = 1/vi)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_point(alpha = 0.65) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.0) +
    annotate("text", x = Inf, y = Inf, label = lab,
             hjust = 1.05, vjust = 1.1, size = 3.4, fontface = "italic") +
    labs(
      title = glue("Meta-regression: lnRR vs {xvar}"),
      subtitle = "Point size = inverse-variance weight (1/vi). Line = fitted trend with 95% CI band.",
      x = xvar,
      y = "ln Response Ratio (lnRR)",
      size = "Weight (1/vi)"
    ) +
    theme_pub(13)
  
  saveRDS(p, file.path(output_dir, paste0("figure_bubble_", safe_file(xvar), ".rds")))
  ggsave(file.path(output_dir, paste0("figure_bubble_", safe_file(xvar), ".png")),
         p, width = 8, height = 5.2, dpi = 600, bg = "white")
}

## ---- FACTOR coefficient plots ----
fac_all <- cand[vapply(dat_es[cand], is.factor, logical(1))]

## ===================== EDIT APPLIED HERE =====================
## Treat requested mods as high-cardinality candidates so they plot robustly.
high_card_vars <- intersect(
  c("crop_type","potassium_source","biochar_feedstock",
    "climate_zone","drought_level","drought_timing","soil_texture"),
  fac_all
)
## =============================================================

topN_highcard <- 10

fac_cand <- union(
  fac_all[vapply(fac_all, function(v) nlevels(dat_es[[v]]) <= max_levels_plot, logical(1))],
  high_card_vars
)

factor_fig_tbl <- tibble()

for (fvar in fac_cand) {
  
  d <- dat_es |> filter(!is.na(.data[[fvar]]), is.finite(yi), is.finite(vi))
  if (nrow(d) < min_k_total || n_distinct(d$study_id) < min_studies) next
  
  lvl_counts <- d |>
    count(.data[[fvar]], name = "k_level") |>
    mutate(level_chr = as.character(.data[[fvar]])) |>
    filter(!is.na(level_chr), level_chr != "") |>
    filter(k_level >= min_k_per_level)
  
  used_topN <- FALSE
  if (fvar %in% high_card_vars) {
    used_topN <- TRUE
    lvl_counts <- lvl_counts |>
      arrange(desc(k_level)) |>
      slice_head(n = topN_highcard)
  }
  
  keep_lvls <- lvl_counts |> pull(level_chr)
  d <- d |> filter(as.character(.data[[fvar]]) %in% keep_lvls)
  if (n_distinct(d[[fvar]]) < 2) next
  
  fml <- as.formula(paste0("~ 0 + ", fvar))
  m <- tryCatch(
    rma.mv(yi, vi, mods = fml,
           random = ~ 1 | study_id/comparison_id,
           data = d, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(m)) next
  
  rob <- tryCatch(cr2_test(m, d$study_id), error = function(e) NULL)
  if (is.null(rob)) next
  
  tt <- tibble(
    term = rownames(rob),
    estimate = rob$beta,
    se = rob$SE,
    p = rob$p_Satt
  ) |>
    mutate(
      level = str_replace(term, paste0("^", fvar), ""),
      level = str_replace_all(level, "^\\s*|\\s*$", ""),
      lo = estimate - 1.96*se,
      hi = estimate + 1.96*se,
      pct_est = pct(estimate),
      pct_lo  = pct(lo),
      pct_hi  = pct(hi),
      sig = sig_code(p),
      k_total = nrow(d),
      studies = n_distinct(d$study_id),
      moderator = fvar,
      topN_flag = used_topN
    ) |>
    arrange(pct_est)
  
  factor_fig_tbl <- bind_rows(factor_fig_tbl, tt)
  
  sub_txt <- if (used_topN) {
    glue("Top-{topN_highcard} most frequent levels (k per level ≥ {min_k_per_level}); k={nrow(d)}, studies={n_distinct(d$study_id)}.")
  } else {
    glue("k={nrow(d)}, studies={n_distinct(d$study_id)}. Levels with k<{min_k_per_level} omitted.")
  }
  
  p_fac <- ggplot(tt, aes(x = pct_est, y = reorder(level, pct_est))) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(size = 2.8) +
    geom_errorbarh(aes(xmin = pct_lo, xmax = pct_hi), height = 0.22, linewidth = 0.6) +
    labs(
      title = glue("Subgroup effects by {fvar} (multilevel; CR2)"),
      subtitle = sub_txt,
      x = "Percent change vs control",
      y = NULL
    ) +
    theme_pub(12)
  
  saveRDS(p_fac, file.path(output_dir, paste0("figure_factor_", safe_file(fvar), ".rds")))
  ggsave(file.path(output_dir, paste0("figure_factor_", safe_file(fvar), ".png")),
         p_fac, width = 9.5, height = 6.0, dpi = 600, bg = "white")
}

write_csv(factor_fig_tbl, file.path(output_dir, "table_factor_moderator_level_effects_CR2.csv"))

## =========================================================
## 5) OVERALL MULTILEVEL MODELS (by outcome_type + by group)
## =========================================================
overall_outcome_tbl <- tibble()
for (ot in levels(dat_es$outcome_type)) {
  d <- dat_es |> filter(outcome_type == ot)
  if (nrow(d) < 3) next
  m <- fit_mv(d)
  rob <- cr2_test(m, d$study_id)
  overall_outcome_tbl <- bind_rows(
    overall_outcome_tbl,
    tibble(
      Outcome_Type = as.character(ot),
      k = nrow(d),
      Studies = n_distinct(d$study_id),
      lnRR = as.numeric(coef(m)),
      CI_low = m$ci.lb,
      CI_high = m$ci.ub,
      Percent = pct(as.numeric(coef(m))),
      p_model = m$pval,
      p_CR2 = rob$p_Satt[1],
      tau2_total = sum(m$sigma2),
      I2_approx = i2_mv_approx(m)
    )
  )
}
write_csv(overall_outcome_tbl, file.path(output_dir, "table_overall_by_outcome_type.csv"))

overall_group_tbl <- tibble()
for (g in levels(dat_es$group)) {
  d <- dat_es |> filter(group == g)
  if (nrow(d) < 3) next
  m <- fit_mv(d)
  rob <- cr2_test(m, d$study_id)
  overall_group_tbl <- bind_rows(
    overall_group_tbl,
    tibble(
      Group = as.character(g),
      k = nrow(d),
      Studies = n_distinct(d$study_id),
      lnRR = as.numeric(coef(m)),
      CI_low = m$ci.lb,
      CI_high = m$ci.ub,
      Percent = pct(as.numeric(coef(m))),
      p_model = m$pval,
      p_CR2 = rob$p_Satt[1],
      tau2_total = sum(m$sigma2),
      I2_approx = i2_mv_approx(m)
    )
  )
}
write_csv(overall_group_tbl, file.path(output_dir, "table_overall_by_group.csv"))

## =========================================================
## 6) FIGURE: GROUP SUMMARY
## =========================================================
group_plot <- overall_group_tbl |>
  mutate(
    SE = (CI_high - CI_low) / (2 * 1.96),
    pct_est = pct(lnRR),
    pct_lo  = pct(CI_low),
    pct_hi  = pct(CI_high),
    sig = sig_code(p_CR2),
    label = glue("{sprintf('%+.1f', pct_est)}%{sig}\n",
                 "{fmt_pct_ci(lnRR, SE)}\n",
                 "(k={k}, studies={Studies})"),
    Group = factor(Group, levels = rev(Group))
  ) |>
  mutate(
    x_pad   = (max(pct_hi, na.rm = TRUE) - min(pct_lo, na.rm = TRUE)) * 0.04,
    label_x = pct_hi + x_pad
  )

p_group <- ggplot(group_plot, aes(y = Group)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = pct_lo, xmax = pct_hi), height = 0.22, linewidth = 0.6) +
  geom_point(aes(x = pct_est), size = 3.1) +
  geom_text(aes(x = label_x, label = label),
            hjust = 0, size = 3.1, lineheight = 0.95) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Pooled percent change (treatment vs control) by broad group",
    subtitle = "95% CI bars. Significance from CR2 robust p-values.",
    x = "Percent change vs control",
    y = NULL
  ) +
  theme_pub(base_size = 12) +
  theme(plot.margin = margin(10, 55, 10, 10))

ggsave(file.path(output_dir, "figure_overall_by_group.png"),
       p_group, width = 10.5, height = 6.5, dpi = 600, bg = "white")

## =========================================================
## 6B) CLASSIC PUBLICATION BIAS (study-aggregated by group)
## =========================================================
agg <- dat_es |>
  group_by(group, study_id) |>
  summarise(
    yi = mean(yi, na.rm = TRUE),
    vi = mean(vi, na.rm = TRUE),
    sei = sqrt(mean(vi, na.rm = TRUE)),
    .groups = "drop"
  )

egger_tbl    <- tibble()
begg_tbl     <- tibble()
trimfill_tbl <- tibble()

for (g in unique(agg$group)) {
  
  dg <- agg |> filter(group == g)
  if (nrow(dg) < 5) next
  
  dg2 <- dg |> arrange(yi)
  m2  <- rma(yi, vi, data = dg2, method = "REML")
  
  tau2 <- as.numeric(m2$tau2)
  I2   <- as.numeric(m2$I2)
  QE_p <- as.numeric(m2$QEp)
  
  eg <- NULL
  if (nrow(dg2) >= 10) {
    eg <- tryCatch(regtest(m2, model = "rma", predictor = "sei"),
                   error = function(e) NULL)
  }
  eg_z <- if (!is.null(eg)) as.numeric(eg$zval) else NA_real_
  eg_p <- if (!is.null(eg)) as.numeric(eg$pval) else NA_real_
  eg_note <- if (nrow(dg2) < 10) "Egger not run (k<10)" else if (is.na(eg_p)) "Egger failed" else glue("Egger: z={round(eg_z,2)}, p={signif(eg_p,3)}")
  
  egger_tbl <- bind_rows(
    egger_tbl,
    tibble(Group = as.character(g), k = nrow(dg2), tau2 = tau2, I2 = I2, QEp = QE_p,
           Egger_z = eg_z, Egger_p = eg_p, Egger_note = eg_note)
  )
  
  bg <- tryCatch(ranktest(m2), error = function(e) NULL)
  bg_tau <- if (!is.null(bg)) as.numeric(bg$tau) else NA_real_
  bg_p   <- if (!is.null(bg)) as.numeric(bg$pval) else NA_real_
  bg_note <- if (is.na(bg_p)) "Begg failed" else glue("Begg: tau={round(bg_tau,2)}, p={signif(bg_p,3)}")
  
  begg_tbl <- bind_rows(
    begg_tbl,
    tibble(Group = as.character(g), k = nrow(dg2), Begg_tau = bg_tau, Begg_p = bg_p, Begg_note = bg_note)
  )
  
  g_file <- safe_file(g)
  
  png(file.path(output_dir, paste0("forest_", g_file, ".png")),
      width = 2200, height = 1500, res = 250)
  par(mar = c(5, 8, 4, 2))
  forest(m2, slab = dg2$study_id, xlab = "Effect size (lnRR)", main = paste("Forest plot –", g), cex = 0.9)
  mtext(glue("k={nrow(dg2)}; τ²={round(tau2,3)}; I²={round(I2,1)}%; Q-test p={signif(QE_p,3)}"),
        side = 3, line = 0.5, adj = 0, cex = 0.9)
  dev.off()
  
  png(file.path(output_dir, paste0("funnel_", g_file, ".png")),
      width = 1400, height = 1050, res = 220)
  funnel(m2, main = paste("Funnel plot –", g), xlab = "lnRR")
  abline(v = coef(m2), lty = 2)
  mtext(glue("k={nrow(dg2)}; dashed = pooled effect"), side = 3, line = 0.5, adj = 0, cex = 0.85)
  mtext(eg_note, side = 3, line = 1.4, adj = 0, cex = 0.85)
  mtext(bg_note, side = 3, line = 2.3, adj = 0, cex = 0.85)
  dev.off()
  
  taf <- tryCatch(trimfill(m2), error = function(e) NULL)
  if (!is.null(taf)) {
    k0 <- tryCatch(as.numeric(taf$k0), error = function(e) NA_real_)
    
    est0 <- as.numeric(coef(m2));  se0 <- as.numeric(m2$se)
    lo0 <- est0 - 1.96 * se0;      hi0 <- est0 + 1.96 * se0
    
    est1 <- as.numeric(coef(taf)); se1 <- as.numeric(taf$se)
    lo1 <- est1 - 1.96 * se1;      hi1 <- est1 + 1.96 * se1
    
    trimfill_tbl <- bind_rows(
      trimfill_tbl,
      tibble(
        Group = as.character(g),
        k_observed = nrow(dg2),
        k0_imputed = k0,
        est_lnRR_original = est0,
        lo_lnRR_original  = lo0,
        hi_lnRR_original  = hi0,
        est_lnRR_trimfill = est1,
        lo_lnRR_trimfill  = lo1,
        hi_lnRR_trimfill  = hi1,
        delta_lnRR = est1 - est0
      )
    )
    
    png(file.path(output_dir, paste0("funnel_trimfill_", g_file, ".png")),
        width = 1400, height = 1050, res = 220)
    funnel(taf, main = paste("Funnel plot (Trim-and-Fill) –", g), xlab = "lnRR")
    abline(v = coef(m2),  lty = 2)
    abline(v = coef(taf), lty = 3)
    mtext(glue("Observed k={nrow(dg2)}; imputed k0={ifelse(is.na(k0),'NA',k0)}"),
          side = 3, line = 0.5, adj = 0, cex = 0.85)
    mtext(glue("Dashed=original ({round(est0,3)}); dotted=trimfill ({round(est1,3)})"),
          side = 3, line = 1.4, adj = 0, cex = 0.85)
    dev.off()
  }
}

write_csv(egger_tbl,    file.path(output_dir, "table_egger_tests_by_group.csv"))
write_csv(begg_tbl,     file.path(output_dir, "table_begg_ranktest_by_group.csv"))
write_csv(trimfill_tbl, file.path(output_dir, "table_trimfill_by_group.csv"))

## =========================================================
## 6C) MULTILEVEL Egger + PET-PEESE + residual funnels (by group)
## =========================================================
egger_mv_tbl <- tibble()
petpeese_tbl <- tibble()

for (g in levels(dat_es$group)) {
  
  d <- dat_es |> filter(group == g, is.finite(yi), is.finite(vi), is.finite(sei))
  if (nrow(d) < 10 || n_distinct(d$study_id) < 5) next
  
  ## multilevel Egger (sei)
  m_egger <- tryCatch(
    rma.mv(yi, vi, mods = ~ sei,
           random = ~ 1 | study_id/comparison_id,
           data = d, method = "REML"),
    error = function(e) NULL
  )
  if (!is.null(m_egger)) {
    rob <- tryCatch(cr2_test(m_egger, d$study_id), error = function(e) NULL)
    if (!is.null(rob) && nrow(rob) >= 2) {
      egger_mv_tbl <- bind_rows(
        egger_mv_tbl,
        tibble(
          Group = as.character(g),
          k = nrow(d),
          Studies = n_distinct(d$study_id),
          Egger_ML_Intercept = rob$beta[1],
          Egger_ML_Int_SE = rob$SE[1],
          Egger_ML_Int_p = rob$p_Satt[1],
          Egger_ML_Slope_sei = rob$beta[2],
          Egger_ML_Slope_SE = rob$SE[2],
          Egger_ML_Slope_p = rob$p_Satt[2]
        )
      )
    }
  }
  
  ## residual funnel from null MV
  m_null <- tryCatch(fit_mv(d), error = function(e) NULL)
  if (!is.null(m_null)) {
    p_rf <- plot_residual_funnel(
      resid = as.numeric(residuals(m_null)),
      sei   = sqrt(d$vi),
      title = paste("Residual funnel (multilevel) –", g),
      subtitle = "Residuals from multilevel RE model; y-axis SE reversed."
    )
    if (!is.null(p_rf)) {
      save_plot_robust(p_rf, paste0("funnel_residual_multilevel_", safe_file(g), ".png"),
                       width = 7.6, height = 5.4, dpi = 450)
    }
  }
  
  ## PET (sei) and PEESE (vi)
  m_pet <- tryCatch(
    rma.mv(yi, vi, mods = ~ sei,
           random = ~ 1 | study_id/comparison_id,
           data = d, method = "REML"),
    error = function(e) NULL
  )
  m_peese <- tryCatch(
    rma.mv(yi, vi, mods = ~ vi,
           random = ~ 1 | study_id/comparison_id,
           data = d, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(m_pet) || is.null(m_peese)) next
  
  r_pet   <- tryCatch(cr2_test(m_pet, d$study_id),   error = function(e) NULL)
  r_peese <- tryCatch(cr2_test(m_peese, d$study_id), error = function(e) NULL)
  if (is.null(r_pet) || is.null(r_peese)) next
  
  pet_int   <- r_pet$beta[1];   pet_int_se <- r_pet$SE[1];   pet_int_p <- r_pet$p_Satt[1]
  pet_slp_p <- if (nrow(r_pet) >= 2) r_pet$p_Satt[2] else NA_real_
  
  peese_int <- r_peese$beta[1]; peese_int_se <- r_peese$SE[1]; peese_int_p <- r_peese$p_Satt[1]
  peese_slp_p <- if (nrow(r_peese) >= 2) r_peese$p_Satt[2] else NA_real_
  
  use_peese <- is.finite(pet_slp_p) && pet_slp_p < 0.05
  chosen_est <- if (use_peese) peese_int else pet_int
  chosen_se  <- if (use_peese) peese_int_se else pet_int_se
  chosen     <- if (use_peese) "PEESE" else "PET"
  
  petpeese_tbl <- bind_rows(
    petpeese_tbl,
    tibble(
      Group = as.character(g),
      k = nrow(d),
      Studies = n_distinct(d$study_id),
      
      PET_int = pet_int,
      PET_int_SE = pet_int_se,
      PET_int_p = pet_int_p,
      PET_slope_p = pet_slp_p,
      
      PEESE_int = peese_int,
      PEESE_int_SE = peese_int_se,
      PEESE_int_p = peese_int_p,
      PEESE_slope_p = peese_slp_p,
      
      PETPEESE_choice = chosen,
      PETPEESE_est = chosen_est,
      PETPEESE_SE = chosen_se,
      PETPEESE_CI = fmt_ci(chosen_est, chosen_se)
    )
  )
}

write_csv(egger_mv_tbl, file.path(output_dir, "table_egger_multilevel_CR2_by_group.csv"))
write_csv(petpeese_tbl, file.path(output_dir, "table_pet_peese_multilevel_CR2_by_group.csv"))

## =========================================================
## 6D) Yang-style sensitivity: FE + VCV rho-imputed + CR2 (if available)
## =========================================================
do_vcv <- "impute_covariance_matrix" %in% getNamespaceExports("clubSandwich")
fe_vcv_tbl <- tibble()

if (!do_vcv) {
  message("NOTE: clubSandwich::impute_covariance_matrix not available -> FE+VCV sensitivity skipped.")
} else {
  
  rho_grid <- c(0.2, 0.5, 0.8)
  
  for (g in levels(dat_es$group)) {
    d <- dat_es |> filter(group == g, is.finite(yi), is.finite(vi))
    if (nrow(d) < 10 || n_distinct(d$study_id) < 5) next
    
    for (rho in rho_grid) {
      V <- tryCatch(
        clubSandwich::impute_covariance_matrix(vi = d$vi, cluster = d$study_id, r = rho, smooth_vi = TRUE),
        error = function(e) NULL
      )
      if (is.null(V)) next
      
      m_fe <- tryCatch(
        rma.mv(yi, V = V, mods = ~ 1, random = ~ 0 | study_id, data = d, method = "FE"),
        error = function(e) NULL
      )
      if (is.null(m_fe)) next
      
      rob <- tryCatch(cr2_test(m_fe, d$study_id), error = function(e) NULL)
      if (is.null(rob)) next
      
      fe_vcv_tbl <- bind_rows(
        fe_vcv_tbl,
        tibble(
          Group = as.character(g),
          rho = rho,
          k = nrow(d),
          Studies = n_distinct(d$study_id),
          FEVCV_est = rob$beta[1],
          FEVCV_SE_CR2 = rob$SE[1],
          FEVCV_p_CR2 = rob$p_Satt[1],
          FEVCV_CI = fmt_ci(rob$beta[1], rob$SE[1])
        )
      )
    }
  }
  
  write_csv(fe_vcv_tbl, file.path(output_dir, "table_FEV_CV_CR2_by_group_rho_sensitivity.csv"))
}

## =========================================================
## 6E) ORCHARD PLOTS (optional; if orchaRd installed)
## =========================================================
orch_log <- file.path(output_dir, "orchard_log.txt")
cat("ORCHARD LOG\n", file = orch_log)

p_orch_g <- NULL
p_orch_o <- NULL

if (!requireNamespace("orchaRd", quietly = TRUE)) {
  cat("Package orchaRd NOT installed; orchard plots skipped.\n", file = orch_log, append = TRUE)
} else {
  
  cat("Package orchaRd detected.\n", file = orch_log, append = TRUE)
  
  ## Group model + orchard
  m_group_mod <- tryCatch(
    metafor::rma.mv(
      yi, vi,
      mods   = ~ 0 + group,
      random = ~ 1 | study_id/comparison_id,
      data   = dat_es,
      method = "REML"
    ),
    error = function(e) {
      cat("m_group_mod FAILED: ", conditionMessage(e), "\n", file = orch_log, append = TRUE)
      NULL
    }
  )
  
  if (!is.null(m_group_mod)) {
    p_orch_g <- tryCatch(
      orchaRd::orchard_plot(
        object = m_group_mod,
        mod    = "group",
        group  = "study_id",
        xlab   = "ln Response Ratio (lnRR)"
      ) +
        theme_pub(13) +
        labs(title = "Orchard plot – Group effects (multilevel model)"),
      error = function(e) {
        cat("orchard_plot(group) FAILED: ", conditionMessage(e), "\n", file = orch_log, append = TRUE)
        NULL
      }
    )
    
    if (!is.null(p_orch_g)) {
      save_plot_robust(p_orch_g, "figure_orchard_group.png", width = 9.2, height = 6.0, dpi = 450)
      cat("Saved: figure_orchard_group.png\n", file = orch_log, append = TRUE)
    }
  }
  
  ## Outcome model + orchard
  m_out_mod <- tryCatch(
    metafor::rma.mv(
      yi, vi,
      mods   = ~ 0 + outcome_type,
      random = ~ 1 | study_id/comparison_id,
      data   = dat_es,
      method = "REML"
    ),
    error = function(e) {
      cat("m_out_mod FAILED: ", conditionMessage(e), "\n", file = orch_log, append = TRUE)
      NULL
    }
  )
  
  if (!is.null(m_out_mod)) {
    p_orch_o <- tryCatch(
      orchaRd::orchard_plot(
        object = m_out_mod,
        mod    = "outcome_type",
        group  = "study_id",
        xlab   = "ln Response Ratio (lnRR)"
      ) +
        theme_pub(13) +
        labs(title = "Orchard plot – Outcome type effects (multilevel model)"),
      error = function(e) {
        cat("orchard_plot(outcome_type) FAILED: ", conditionMessage(e), "\n", file = orch_log, append = TRUE)
        NULL
      }
    )
    
    if (!is.null(p_orch_o)) {
      save_plot_robust(p_orch_o, "figure_orchard_outcome_type.png", width = 10.5, height = 7.0, dpi = 450)
      cat("Saved: figure_orchard_outcome_type.png\n", file = orch_log, append = TRUE)
    }
  }
  
  cat("Done.\n", file = orch_log, append = TRUE)
}

## =========================================================
## 7) MULTIVARIABLE MODERATOR MODEL (AUTO-REDUCE UNTIL FITS)
## =========================================================

non_mods <- c(
  "study_id","comparison_id",
  "outcome_type","group","outcome_variable","outcome_unit",
  "treatment_mean","treatment_sd","treatment_n",
  "control_mean","control_sd","control_n",
  "yi","vi","sei"
)

preferred_mods <- c(
  "climate_zone",
  "study_type",                 # from experiment_type
  "crop_type",
  "growth_stage_treated",
  "drought_level",
  "drought_description",
  "drought_timing",
  "soil_texture",
  "soil_organic_matter_percent_z",
  "soil_ec_ds_m_z",
  "soil_p_h_z",
  "biochar_applied",
  "biochar_feedstock",
  "biochar_rate_t_ha_z",
  "biochar_p_h_z",
  "biochar_cec_z",
  "biochar_toc_percent_z",
  "potassium_source"
)

preferred_mods <- preferred_mods[preferred_mods %in% names(dat_es)]
preferred_mods <- setdiff(preferred_mods, non_mods)
if (length(preferred_mods) == 0) stop("No preferred moderators found after cleaning/harmonization.")

## ---- MV-only recode/collapse helpers ----
collapse_for_mv <- function(df, var, topN = 8, min_n = 2, other = "Other", missing = "Unknown") {
  if (!var %in% names(df)) return(df)
  
  x <- as.character(df[[var]])
  x[is.na(x) | x == ""] <- missing
  x <- factor(x)
  
  ## topN lump, then lump tiny levels
  x <- forcats::fct_lump_n(x, n = topN, other_level = other)
  x <- forcats::fct_lump_min(x, min = min_n, other_level = other)
  
  df[[var]] <- droplevels(x)
  df
}

fixed_df_burden <- function(df, mods) {
  mods <- mods[mods %in% names(df)]
  if (length(mods) == 0) return(0L)
  
  sum(vapply(mods, function(v) {
    if (is.factor(df[[v]])) max(0L, nlevels(df[[v]]) - 1L) else 1L
  }, integer(1)))
}

fit_mv_safe2 <- function(df, fml, method = "REML") {
  m <- fit_mv_safe(df, mods = fml)
  if (!is.null(m)) return(m)
  
  ## fallback: ML sometimes converges when REML doesn’t
  m2 <- tryCatch(fit_mv_safe(df, mods = fml, method = "ML"), error = function(e) NULL)
  m2
}

## ---- univariate screen (optional, but keep robust) ----
screen_tbl <- tibble(
  Moderator = character(),
  Term = character(),
  beta = numeric(),
  SE = numeric(),
  p_CR2 = numeric(),
  k = integer(),
  Studies = integer()
)

for (v in preferred_mods) {
  d <- dat_es |> filter(!is.na(.data[[v]]))
  if (nrow(d) < 10) next
  if (n_distinct(d$study_id) < 5) next
  if (length(unique(d[[v]])) < 2) next
  
  fml <- as.formula(paste0("~ ", v))
  
  tmp <- tryCatch({
    m  <- fit_mv(d, mods = fml)
    rb <- cr2_test(m, d$study_id)
    tibble(
      Moderator = v,
      Term = rownames(rb),
      beta = rb$beta,
      SE = rb$SE,
      p_CR2 = rb$p_Satt,
      k = nrow(d),
      Studies = n_distinct(d$study_id)
    )
  }, error = function(e) NULL)
  
  if (!is.null(tmp)) screen_tbl <- bind_rows(screen_tbl, tmp)
}

write_csv(screen_tbl, file.path(output_dir, "table_moderator_screen_univariate_CR2.csv"))

mods_keep <- NULL
if (nrow(screen_tbl) > 0) {
  mods_keep <- screen_tbl |>
    filter(Term != "Intrcpt") |>
    group_by(Moderator) |>
    summarise(min_p = min(p_CR2, na.rm = TRUE), .groups = "drop") |>
    filter(min_p < 0.20) |>
    pull(Moderator)
}

if (is.null(mods_keep) || length(mods_keep) == 0) {
  ## sane fallback given your CWPS missingness profile
  mods_keep <- intersect(
    preferred_mods,
    c("study_type","drought_level","biochar_applied",
      "soil_p_h_z","soil_ec_ds_m_z","soil_organic_matter_percent_z","climate_zone")
  )
  if (length(mods_keep) == 0) mods_keep <- preferred_mods[1:min(6, length(preferred_mods))]
}

## Separate numeric vs factors; drop highly correlated numerics
is_num <- function(v) vapply(v, function(nm) is.numeric(dat_es[[nm]]), logical(1))
num_mods_keep <- mods_keep[is_num(mods_keep)]
fac_mods      <- setdiff(mods_keep, num_mods_keep)
num_mods2     <- drop_high_cor(dat_es, num_mods_keep, cutoff = 0.70)

mods_final_mv <- c(fac_mods, num_mods2)
mods_final_mv <- mods_final_mv[mods_final_mv %in% names(dat_es)]

## ---- MV-only collapsed copy (this is key) ----
dat_es_mv <- dat_es

## Collapse the worst offenders (and sparse ones)
dat_es_mv <- collapse_for_mv(dat_es_mv, "crop_type",            topN = 6, min_n = 2)
dat_es_mv <- collapse_for_mv(dat_es_mv, "drought_description",  topN = 6, min_n = 2)
dat_es_mv <- collapse_for_mv(dat_es_mv, "climate_zone",         topN = 8, min_n = 2)
dat_es_mv <- collapse_for_mv(dat_es_mv, "potassium_source",     topN = 6, min_n = 2)
dat_es_mv <- collapse_for_mv(dat_es_mv, "biochar_feedstock",    topN = 5, min_n = 2)
dat_es_mv <- collapse_for_mv(dat_es_mv, "soil_texture",         topN = 6, min_n = 2)
dat_es_mv <- collapse_for_mv(dat_es_mv, "growth_stage_treated", topN = 5, min_n = 2)
dat_es_mv <- collapse_for_mv(dat_es_mv, "drought_timing",       topN = 5, min_n = 2)

## ---- AUTO-REDUCE LOOP: drop moderators until the model fits ----
missing_rate <- function(df, v) mean(is.na(df[[v]]))

attempt_fit <- function(mods) {
  mods <- mods[mods %in% names(dat_es_mv)]
  if (length(mods) == 0) return(NULL)
  
  ## ensure enough complete cases
  res <- tryCatch(
    ensure_complete(dat_es_mv, mods, min_k = 12, min_studies = 6),
    error = function(e) NULL
  )
  if (is.null(res)) return(NULL)
  
  d_mv <- res$data
  mods <- res$mods
  
  ## drop factor moderators that collapse to 1 level after CC
  bad_fac <- mods[vapply(mods, function(v) is.factor(d_mv[[v]]) && n_distinct(d_mv[[v]]) < 2, logical(1))]
  if (length(bad_fac) > 0) {
    mods <- setdiff(mods, bad_fac)
    d_mv <- dat_es_mv |> filter(if_all(all_of(mods), ~ !is.na(.x)))
  }
  
  ## df burden check: keep fixed df comfortably below k
  df_burden <- fixed_df_burden(d_mv, mods)
  k_cc <- nrow(d_mv)
  
  ## be conservative: require df_burden <= k_cc - 8
  if (df_burden > (k_cc - 8)) return(list(ok = FALSE, reason = "too_many_df", mods = mods, d_mv = d_mv))
  
  mv_formula <- as.formula(paste0("~ ", paste(mods, collapse = " + ")))
  m <- fit_mv_safe2(d_mv, mv_formula, method = "REML")
  if (is.null(m)) return(list(ok = FALSE, reason = "fit_failed", mods = mods, d_mv = d_mv))
  
  list(ok = TRUE, mods = mods, d_mv = d_mv, mv_formula = mv_formula, m = m)
}

mods_working <- mods_final_mv

## As a hard ceiling, never keep more than 6 factor moderators at once (small k)
reduce_factor_count <- function(mods, df, max_fac = 4) {
  fac <- mods[mods %in% names(df) & vapply(mods, function(v) is.factor(df[[v]]), logical(1))]
  num <- setdiff(mods, fac)
  if (length(fac) <= max_fac) return(mods)
  
  ## keep factor mods with least missingness first
  fac2 <- fac[order(vapply(fac, function(v) missing_rate(df, v), numeric(1)))]
  c(fac2[1:max_fac], num)
}

mods_working <- reduce_factor_count(mods_working, dat_es_mv, max_fac = 4)

fit_res <- NULL
for (iter in 1:30) {
  fit_res <- attempt_fit(mods_working)
  if (!is.null(fit_res) && isTRUE(fit_res$ok)) break
  
  if (is.null(fit_res)) {
    ## if even CC fails, drop the highest-missing moderator
    mr <- vapply(mods_working, function(v) missing_rate(dat_es_mv, v), numeric(1))
    drop_one <- names(which.max(mr))[1]
    message("Auto-reduce: dropping '", drop_one, "' (missingness=", round(100*max(mr),1), "%) because CC failed.")
    mods_working <- setdiff(mods_working, drop_one)
    next
  }
  
  d_mv_tmp <- fit_res$d_mv
  mods_tmp <- fit_res$mods
  
  ## choose drop candidate:
  ## 1) biggest df (factor with most levels) within CC data
  df_each <- vapply(mods_tmp, function(v) {
    if (is.factor(d_mv_tmp[[v]])) max(0L, nlevels(d_mv_tmp[[v]]) - 1L) else 1L
  }, integer(1))
  
  ## prioritize dropping big-df factors; otherwise drop highest missingness
  fac_candidates <- names(df_each)[df_each > 1]
  if (length(fac_candidates) > 0) {
    drop_one <- fac_candidates[which.max(df_each[fac_candidates])]
    message("Auto-reduce: dropping '", drop_one, "' (df=", df_each[drop_one], ") to reduce fixed-effect burden.")
  } else {
    mr <- vapply(mods_tmp, function(v) missing_rate(dat_es_mv, v), numeric(1))
    drop_one <- names(which.max(mr))[1]
    message("Auto-reduce: dropping '", drop_one, "' (missingness=", round(100*mr[drop_one],1), "%) after fit failure.")
  }
  
  mods_working <- setdiff(mods_working, drop_one)
  mods_working <- reduce_factor_count(mods_working, dat_es_mv, max_fac = 4)
  
  if (length(mods_working) < 2) break
}

if (is.null(fit_res) || !isTRUE(fit_res$ok)) {
  ## last-ditch core model (should fit for CWPS)
  message("Auto-reduce failed; fitting last-ditch core MV model.")
  core <- intersect(
    names(dat_es_mv),
    c("study_type","drought_level","biochar_applied","soil_p_h_z","soil_ec_ds_m_z","soil_organic_matter_percent_z","climate_zone")
  )
  core <- core[core %in% names(dat_es_mv)]
  res2 <- ensure_complete(dat_es_mv, core, min_k = 12, min_studies = 6)
  d_mv <- res2$data
  core <- res2$mods
  mv_formula <- as.formula(paste0("~ ", paste(core, collapse = " + ")))
  m_mv <- fit_mv_safe2(d_mv, mv_formula, method = "REML")
  if (is.null(m_mv)) stop("Multivariable model failed even in last-ditch core form.")
  mods_final_mv <- core
} else {
  mods_final_mv <- fit_res$mods
  d_mv <- fit_res$d_mv
  mv_formula <- fit_res$mv_formula
  m_mv <- fit_res$m
}

message("Final MV moderators: ", paste(mods_final_mv, collapse = " + "))
message("MV complete-case k = ", nrow(d_mv), "; studies = ", n_distinct(d_mv$study_id))
message("Fixed-effect df burden ≈ ", fixed_df_burden(d_mv, mods_final_mv))

## Null model on the same CC data for R2_tau / AIC
m0_mv <- fit_mv(d_mv)

rob_mv <- cr2_test(m_mv, d_mv$study_id)

aic0_mv <- tryCatch(AIC(m0_mv), error = function(e) NA_real_)
aic1_mv <- tryCatch(AIC(m_mv),  error = function(e) NA_real_)
dAIC_mv <- aic1_mv - aic0_mv
r2_mv   <- R2_tau(m0_mv, m_mv)

mv_tbl <- tibble(
  Term = rownames(rob_mv),
  Estimate = rob_mv$beta,
  SE_CR2 = rob_mv$SE,
  p_CR2 = rob_mv$p_Satt,
  Sig = sig_code(rob_mv$p_Satt),
  CI = fmt_ci(rob_mv$beta, rob_mv$SE),
  k = nrow(d_mv),
  Studies = n_distinct(d_mv$study_id),
  AIC_null = aic0_mv,
  AIC_mod  = aic1_mv,
  Delta_AIC = dAIC_mv,
  R2_tau = r2_mv,
  tau2_null = tau2_total(m0_mv),
  tau2_mod  = tau2_total(m_mv)
)

write_csv(mv_tbl, file.path(output_dir, "table_moderator_multivariable_CR2.csv"))
save_docx_table(
  mv_tbl |> select(Term, CI, p_CR2, Sig, k, Studies, R2_tau, Delta_AIC),
  "Moderator effects (MULTIVARIABLE; multilevel; CR2 robust) + fit stats",
  file.path(output_dir, "table_moderator_multivariable_CR2.docx")
)

## VIF on MV design matrix (if feasible)
X <- model.matrix(mv_formula, data = d_mv)
X2 <- if ("(Intercept)" %in% colnames(X)) X[, colnames(X) != "(Intercept)", drop = FALSE] else X
vif_tbl <- if (ncol(X2) >= 2) vif_manual(X2) else tibble(Term = colnames(X2), VIF = NA_real_)
write_csv(vif_tbl, file.path(output_dir, "table_vif_multivariable.csv"))

## =========================================================
## 8) LEAVE-ONE-STUDY-OUT (OVERALL)
## =========================================================
loco_df <- tibble()
for (sid in unique(dat_es$study_id)) {
  dsub <- dat_es |> filter(study_id != sid)
  if (nrow(dsub) < 3) next
  msub <- fit_mv(dsub)
  loco_df <- bind_rows(loco_df, tibble(left_out = as.character(sid), lnRR = as.numeric(coef(msub))))
}

m_full <- fit_mv(dat_es)
full_est <- as.numeric(coef(m_full))
full_se  <- sqrt(as.numeric(vcov(m_full)))
full_lo  <- full_est - 1.96 * full_se
full_hi  <- full_est + 1.96 * full_se

p_loco <- ggplot(loco_df, aes(x = reorder(left_out, lnRR), y = lnRR)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = full_lo, ymax = full_hi),
            inherit.aes = FALSE, alpha = 0.08) +
  geom_hline(yintercept = full_est, linetype = "dashed") +
  geom_point(size = 2) +
  coord_flip() +
  labs(
    title = "Leave-one-study-out sensitivity",
    subtitle = glue("Dashed = pooled lnRR; shaded band = 95% CI [{round(full_lo,3)}, {round(full_hi,3)}]"),
    x = "Study removed",
    y = "Pooled lnRR"
  ) +
  theme_pub(base_size = 13)

ggsave(file.path(output_dir, "figure_leave_one_out.png"),
       p_loco, width = 7.8, height = 6.2, dpi = 600, bg = "white")

loco_summary <- tibble(
  pooled_lnRR = full_est,
  pooled_lnRR_lo = full_lo,
  pooled_lnRR_hi = full_hi,
  pooled_percent = pct(full_est),
  pooled_percent_lo = pct(full_lo),
  pooled_percent_hi = pct(full_hi),
  loo_min_lnRR = min(loco_df$lnRR, na.rm = TRUE),
  loo_max_lnRR = max(loco_df$lnRR, na.rm = TRUE),
  loo_range_lnRR = max(loco_df$lnRR, na.rm = TRUE) - min(loco_df$lnRR, na.rm = TRUE),
  loo_min_percent = pct(min(loco_df$lnRR, na.rm = TRUE)),
  loo_max_percent = pct(max(loco_df$lnRR, na.rm = TRUE))
)
write_csv(loco_summary, file.path(output_dir, "table_leave_one_out_summary.csv"))

## =========================================================
## 9) PATCHWORK COMPOSITES
##   (A) Auto composites from all .rds panels (batched)
##   (B) Requested Figure_4–Figure_8 composites with alphabet tags
## =========================================================

## ---- (A) AUTO composites from RDS ----
make_composite_from_rds <- function(files, out_stub, ncol = 2, width = 12, height = 9) {
  if (length(files) == 0) return(invisible(NULL))
  plots <- lapply(files, read_plot_rds)
  p <- wrap_plots(plots, ncol = ncol) + plot_annotation(tag_levels = "a")
  
  ggsave(fs::path(output_dir, paste0(out_stub, ".pdf")),
         plot = p, width = width, height = height, units = "in",
         device = pdf_device)
  
  ggsave(fs::path(output_dir, paste0(out_stub, ".png")),
         plot = p, width = width, height = height, units = "in",
         dpi = 600, bg = "white")
  
  invisible(p)
}

rds_factor <- sort(fs::dir_ls(output_dir, glob = "figure_factor_*.rds"))
rds_bubble <- sort(fs::dir_ls(output_dir, glob = "figure_bubble_*.rds"))

batch <- function(x, k = 4) split(x, ceiling(seq_along(x) / k))

factor_batches <- batch(rds_factor, k = 4)
if (length(factor_batches) > 0) {
  i <- 1L
  for (files in factor_batches) {
    make_composite_from_rds(files, out_stub = glue("composite_factor_panels_{i}"),
                            ncol = 2, width = 12, height = 9)
    i <- i + 1L
  }
}

bubble_batches <- batch(rds_bubble, k = 4)
if (length(bubble_batches) > 0) {
  i <- 1L
  for (files in bubble_batches) {
    make_composite_from_rds(files, out_stub = glue("composite_bubble_panels_{i}"),
                            ncol = 2, width = 12, height = 9)
    i <- i + 1L
  }
}

## ---- (B) Requested composites: Figure_4–Figure_8 ----
read_rds_stub <- function(stub) {
  # stub without extension, in output_dir
  path <- fs::path(output_dir, paste0(stub, ".rds"))
  if (!file.exists(path)) stop("Missing RDS: ", path)
  readRDS(path)
}

read_png_stub <- function(stub) {
  # stub without extension, in output_dir; returns patchwork element
  path <- fs::path(output_dir, paste0(stub, ".png"))
  if (!file.exists(path)) stop("Missing PNG: ", path)
  
  if (requireNamespace("magick", quietly = TRUE)) {
    img <- magick::image_read(path)
    grob <- grid::rasterGrob(as.raster(img), interpolate = TRUE)
    return(patchwork::wrap_elements(full = grob))
  }
  if (requireNamespace("png", quietly = TRUE)) {
    img <- png::readPNG(path)
    grob <- grid::rasterGrob(img, interpolate = TRUE)
    return(patchwork::wrap_elements(full = grob))
  }
  stop("Need either {magick} or {png} installed to compose PNG panels. Install one of them.")
}

save_composite <- function(p, out_stub, width, height) {
  ggsave(fs::path(output_dir, paste0(out_stub, ".pdf")),
         plot = p, width = width, height = height, units = "in",
         device = pdf_device)
  ggsave(fs::path(output_dir, paste0(out_stub, ".png")),
         plot = p, width = width, height = height, units = "in",
         dpi = 600, bg = "white")
}

## Figure 4: factor drought_timing + drought_level
Figure_4 <- wrap_plots(
  read_rds_stub("figure_factor_drought_timing"),
  read_rds_stub("figure_factor_drought_level"),
  ncol = 1
) + plot_annotation(tag_levels = "a")
save_composite(Figure_4, "Figure_4", width = 10.5, height = 11)

## Figure 5: factor climate_zone + soil_texture
Figure_5 <- wrap_plots(
  read_rds_stub("figure_factor_climate_zone"),
  read_rds_stub("figure_factor_soil_texture"),
  ncol = 1
) + plot_annotation(tag_levels = "a")
save_composite(Figure_5, "Figure_5", width = 10.5, height = 11)

## Figure 6: bubble soil_organic_matter_percent_z + soil_p_h_z + soil_ec_ds_m_z
Figure_6 <- wrap_plots(
  read_rds_stub("figure_bubble_soil_organic_matter_percent_z"),
  read_rds_stub("figure_bubble_soil_p_h_z"),
  read_rds_stub("figure_bubble_soil_ec_ds_m_z"),
  ncol = 1
) + plot_annotation(tag_levels = "a")
save_composite(Figure_6, "Figure_6", width = 10.5, height = 15)

## Figure 7: bubble biochar_rate_t_ha_z + biochar_p_h_z + biochar_cec_z
Figure_7 <- wrap_plots(
  read_rds_stub("figure_bubble_biochar_rate_t_ha_z"),
  read_rds_stub("figure_bubble_biochar_p_h_z"),
  read_rds_stub("figure_bubble_biochar_cec_z"),
  ncol = 1
) + plot_annotation(tag_levels = "a")
save_composite(Figure_7, "Figure_7", width = 10.5, height = 15)

## Figure 8: leave-one-out + funnels for Other + Water_relations (PNG panels)
Figure_8 <- wrap_plots(
  read_png_stub("figure_leave_one_out"),
  read_png_stub("funnel_Other"),
  read_png_stub("funnel_residual_multilevel_Other"),
  read_png_stub("funnel_trimfill_Other"),
  read_png_stub("funnel_trimfill_Water_relations"),
  read_png_stub("funnel_Water_relations"),
  ncol = 2
) + plot_annotation(tag_levels = "a")
save_composite(Figure_8, "Figure_8", width = 13, height = 10)

## =========================================================
## 10) SOFTWARE VERSIONS + SESSION INFO
## =========================================================
pkgs_used <- c("tidyverse", "metafor", "clubSandwich", "janitor", "readr", "fs", "glue",
               "flextable", "officer", "stringr", "patchwork")
pkg_tbl <- tibble(
  R = R.version.string,
  Package = pkgs_used,
  Version = map_chr(pkgs_used, ~ as.character(packageVersion(.x)))
)
write_csv(pkg_tbl, file.path(output_dir, "table_software_versions.csv"))
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "sessionInfo.txt"))

message("✅ Meta-analysis completed successfully.")
message("✅ Outputs saved in: ", output_dir)
message("✅ Requested composites saved: Figure_4–Figure_8 (PDF + PNG) in ", output_dir)
message("✅ Orchard plots attempted; see orchard_log.txt for details.")
