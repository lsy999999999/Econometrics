options(stringsAsFactors = FALSE)

this_file <- tryCatch(normalizePath(sys.frames()[[1]]$ofile), error = function(e) NA_character_)
if (is.na(this_file)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) this_file <- normalizePath(sub("^--file=", "", file_arg[1]))
}

project_root <- normalizePath(file.path(dirname(this_file), ".."))
step3_data_file <- file.path(project_root, "output", "step3", "step3_cross_section_data.csv")
out_dir <- file.path(project_root, "output", "step4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(step3_data_file)) {
  stop("Missing Step 3 data. Run code/step3_multivariate_cross_section.R first.")
}

df <- read.csv(step3_data_file)

base_formula <- CAR ~ RateVulnerability_z + MktBeta_z + SMBBeta_z + HMLBeta_z +
  RMWBeta_z + CMABeta_z + HistVol_z + PreMomentum_z

base_terms <- c(
  "RateVulnerability_z", "MktBeta_z", "SMBBeta_z", "HMLBeta_z",
  "RMWBeta_z", "CMABeta_z", "HistVol_z", "PreMomentum_z"
)

control_terms <- setdiff(base_terms, "RateVulnerability_z")

fit_ols <- lm(base_formula, data = df)
df$ols_fitted <- fitted(fit_ols)
df$ols_resid <- resid(fit_ols)
df$ols_abs_resid <- abs(df$ols_resid)
df$ols_resid_sq <- df$ols_resid^2
df$hat <- hatvalues(fit_ols)
df$cooks_d <- cooks.distance(fit_ols)
df$dfbetas_rate_vulnerability <- dfbetas(fit_ols)[, "RateVulnerability_z"]
df$rv_quartile <- cut(
  df$RateVulnerability_z,
  breaks = quantile(df$RateVulnerability_z, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("Q1 lowest", "Q2", "Q3", "Q4 highest")
)

weighted_vcov <- function(model, type = "const") {
  x <- model.matrix(model)
  model_weights <- weights(model)
  if (is.null(model_weights)) model_weights <- rep(1, nrow(x))
  residual <- residuals(model)
  bread <- solve(t(x) %*% (model_weights * x))
  if (type == "const") return(vcov(model))

  hat_values <- hatvalues(model)
  adjustment <- switch(
    type,
    HC0 = rep(1, length(residual)),
    HC1 = rep(nrow(x) / (nrow(x) - ncol(x)), length(residual)),
    HC2 = 1 / pmax(1 - hat_values, .Machine$double.eps),
    HC3 = 1 / pmax(1 - hat_values, .Machine$double.eps)^2,
    stop("Unsupported covariance type: ", type)
  )
  meat <- t(x) %*% ((model_weights^2 * residual^2 * adjustment) * x)
  bread %*% meat %*% bread
}

coefficient_table <- function(model, model_name, se_type = "const", term_filter = NULL) {
  vc <- weighted_vcov(model, se_type)
  estimates <- coef(model)
  se <- sqrt(diag(vc))
  t_value <- estimates / se
  df_resid <- df.residual(model)
  p_value <- 2 * pt(abs(t_value), df = df_resid, lower.tail = FALSE)
  ci_low <- estimates + qt(0.025, df = df_resid) * se
  ci_high <- estimates + qt(0.975, df = df_resid) * se
  out <- data.frame(
    model = model_name,
    se_type = se_type,
    term = names(estimates),
    estimate = unname(estimates),
    std_error = unname(se),
    t_value = unname(t_value),
    p_value = unname(p_value),
    conf_low_95 = unname(ci_low),
    conf_high_95 = unname(ci_high),
    row.names = NULL
  )
  if (!is.null(term_filter)) out <- out[out$term %in% term_filter, ]
  out
}

lm_test <- function(aux_model, test_name) {
  n <- nobs(aux_model)
  r2 <- summary(aux_model)$r.squared
  df_test <- aux_model$rank - 1
  statistic <- n * r2
  data.frame(
    test = test_name,
    statistic = statistic,
    df = df_test,
    p_value = pchisq(statistic, df = df_test, lower.tail = FALSE),
    method = "LM = n * auxiliary R-squared",
    row.names = NULL
  )
}

bp_test <- function(model, data, var_terms, test_name) {
  data$u2 <- resid(model)^2
  aux_formula <- as.formula(paste("u2 ~", paste(var_terms, collapse = " + ")))
  lm_test(lm(aux_formula, data = data), test_name)
}

white_fitted_test <- function(model, data, fitted_name, test_name) {
  data$u2 <- resid(model)^2
  data[[fitted_name]] <- fitted(model)
  aux_formula <- as.formula(paste("u2 ~", fitted_name, "+ I(", fitted_name, "^2)"))
  lm_test(lm(aux_formula, data = data), test_name)
}

make_white_full_data <- function(data, terms) {
  out <- data
  for (term in terms) out[[paste0(term, "_sq")]] <- out[[term]]^2
  if (length(terms) >= 2) {
    for (i in seq_len(length(terms) - 1)) {
      for (j in (i + 1):length(terms)) {
        name <- paste0(terms[i], "_x_", terms[j])
        out[[name]] <- out[[terms[i]]] * out[[terms[j]]]
      }
    }
  }
  out
}

white_full_test <- function(model, data, terms, test_name) {
  white_data <- make_white_full_data(data, terms)
  white_data$u2 <- resid(model)^2
  square_terms <- paste0(terms, "_sq")
  interaction_terms <- names(white_data)[grepl("_x_", names(white_data), fixed = TRUE)]
  rhs <- paste(c(terms, square_terms, interaction_terms), collapse = " + ")
  aux <- lm(as.formula(paste("u2 ~", rhs)), data = white_data)
  lm_test(aux, test_name)
}

robust_wald_test <- function(model, terms, vcov_type = "HC3", test_name) {
  coefficients <- coef(model)
  vc <- weighted_vcov(model, vcov_type)
  missing_terms <- setdiff(terms, names(coefficients))
  if (length(missing_terms) > 0) stop("Terms not found: ", paste(missing_terms, collapse = ", "))

  r_matrix <- matrix(0, nrow = length(terms), ncol = length(coefficients))
  colnames(r_matrix) <- names(coefficients)
  for (i in seq_along(terms)) r_matrix[i, terms[i]] <- 1

  restrictions <- as.vector(r_matrix %*% coefficients)
  restricted_vcov <- r_matrix %*% vc %*% t(r_matrix)
  wald_chi2 <- as.numeric(t(restrictions) %*% solve(restricted_vcov) %*% restrictions)
  q <- length(terms)
  f_stat <- wald_chi2 / q
  data.frame(
    test = test_name,
    vcov_type = vcov_type,
    terms = paste(terms, collapse = " ; "),
    chi_square = wald_chi2,
    chi_square_df = q,
    chi_square_p = pchisq(wald_chi2, df = q, lower.tail = FALSE),
    F_statistic = f_stat,
    F_df1 = q,
    F_df2 = df.residual(model),
    F_p_value = pf(f_stat, df1 = q, df2 = df.residual(model), lower.tail = FALSE),
    row.names = NULL
  )
}

fit_m2_hc3_rate <- function(data, model_name) {
  model <- lm(base_formula, data = data)
  coefficient_table(model, model_name, "HC3", "RateVulnerability_z")
}

se_types <- c("const", "HC0", "HC1", "HC2", "HC3")
ols_standard_vs_hc <- do.call(
  rbind,
  lapply(se_types, function(se_type) {
    coefficient_table(fit_ols, "OLS_M2", se_type, term_filter = "RateVulnerability_z")
  })
)

hetero_tests <- do.call(
  rbind,
  list(
    bp_test(fit_ols, df, base_terms, "OLS M2 Breusch-Pagan: u^2 on all regressors"),
    white_fitted_test(fit_ols, df, "ols_fitted", "OLS M2 modified White: u^2 on fitted and fitted^2"),
    white_full_test(fit_ols, df, base_terms, "OLS M2 full White: u^2 on regressors, squares, and interactions")
  )
)

robust_wald_tests <- do.call(
  rbind,
  list(
    robust_wald_test(fit_ols, control_terms, "HC3", "OLS M2 controls jointly zero"),
    robust_wald_test(fit_ols, c("SMBBeta_z", "HMLBeta_z", "RMWBeta_z", "CMABeta_z"), "HC3", "OLS M2 Fama-French style betas jointly zero")
  )
)

wls_weights <- 1 / pmax(df$HistVol^2, .Machine$double.eps)
fit_wls <- lm(base_formula, data = df, weights = wls_weights)

eps <- 1e-8
df$log_ols_resid_sq <- log(df$ols_resid^2 + eps)
variance_aux <- lm(log_ols_resid_sq ~ RateVulnerability_z + HistVol_z + PreMomentum_z, data = df)
df$fgls_h_hat <- exp(fitted(variance_aux))
fgls_weights <- 1 / pmax(df$fgls_h_hat, .Machine$double.eps)
fit_fgls <- lm(base_formula, data = df, weights = fgls_weights)
df$fgls_fitted <- fitted(fit_fgls)

variance_aux_white <- lm(log_ols_resid_sq ~ ols_fitted + I(ols_fitted^2), data = df)
df$h_white <- exp(fitted(variance_aux_white))
df$w_white <- 1 / pmax(df$h_white, .Machine$double.eps)
df$w_white <- df$w_white / mean(df$w_white)
fit_fgls_white <- lm(base_formula, data = df, weights = w_white)
df$fgls_white_fitted <- fitted(fit_fgls_white)

quartile_var <- aggregate(ols_resid ~ rv_quartile, data = df, FUN = function(x) var(x, na.rm = TRUE))
names(quartile_var) <- c("rv_quartile", "residual_variance")
quartile_var$n <- as.integer(table(df$rv_quartile)[as.character(quartile_var$rv_quartile)])
df <- merge(df, quartile_var[, c("rv_quartile", "residual_variance")], by = "rv_quartile", all.x = TRUE, sort = FALSE)
names(df)[names(df) == "residual_variance"] <- "h_group"
df$w_group <- 1 / pmax(df$h_group, .Machine$double.eps)
df$w_group <- df$w_group / mean(df$w_group)
df <- df[order(match(df$industry, read.csv(step3_data_file)$industry)), ]
fit_wls_group <- lm(base_formula, data = df, weights = w_group)
df$wls_group_fitted <- fitted(fit_wls_group)

wls_fgls_comparison <- do.call(
  rbind,
  list(
    coefficient_table(fit_ols, "OLS_M2", "const", "RateVulnerability_z"),
    coefficient_table(fit_ols, "OLS_M2", "HC3", "RateVulnerability_z"),
    coefficient_table(fit_wls, "WLS_weight_1_over_HistVol_sq", "const", "RateVulnerability_z"),
    coefficient_table(fit_wls, "WLS_weight_1_over_HistVol_sq", "HC3", "RateVulnerability_z"),
    coefficient_table(fit_fgls, "FGLS_log_variance_function", "const", "RateVulnerability_z"),
    coefficient_table(fit_fgls, "FGLS_log_variance_function", "HC3", "RateVulnerability_z"),
    coefficient_table(fit_fgls_white, "FGLS_White_fitted_and_fitted_sq", "const", "RateVulnerability_z"),
    coefficient_table(fit_fgls_white, "FGLS_White_fitted_and_fitted_sq", "HC3", "RateVulnerability_z"),
    coefficient_table(fit_wls_group, "Group_WLS_RV_quartile_variance", "const", "RateVulnerability_z"),
    coefficient_table(fit_wls_group, "Group_WLS_RV_quartile_variance", "HC3", "RateVulnerability_z")
  )
)

post_weight_tests <- do.call(
  rbind,
  list(
    white_fitted_test(fit_wls, df, "wls_fitted", "WLS modified White: u^2 on fitted and fitted^2"),
    white_fitted_test(fit_fgls, df, "fgls_fitted", "FGLS modified White: u^2 on fitted and fitted^2"),
    white_fitted_test(fit_fgls_white, df, "fgls_white_fitted", "FGLS-White post-test: u^2 on fitted and fitted^2"),
    white_fitted_test(fit_wls_group, df, "wls_group_fitted", "Group-WLS post-test: u^2 on fitted and fitted^2")
  )
)
hetero_tests <- rbind(hetero_tests, post_weight_tests)

model_decision <- c(
  "Step 4 fixes the Step 3 M2 full-controls model and focuses on inference under possible heteroskedasticity.",
  "In this cross-section, industry portfolios differ in pre-event volatility and macro sensitivity, so homoskedasticity is not guaranteed.",
  "The main benchmark remains OLS point estimates with heteroskedasticity-robust standard errors, especially HC3 because N = 49 is small.",
  "WLS uses 1 / HistVol^2 as an economically interpretable precision weight; it is treated as a robustness check, not the main model.",
  "FGLS estimates log residual variance using RateVulnerability_z, HistVol_z, and PreMomentum_z, then applies 1 / h_hat weights.",
  "FGLS-White estimates log residual variance using fitted values and fitted values squared, matching the modified White diagnostic logic.",
  "Groupwise WLS uses rate-vulnerability quartile residual variances as weights, matching the residual variance plot.",
  "If WLS/FGLS post-tests still suggest heteroskedasticity, use robust WLS/FGLS standard errors and do not replace the OLS-HC conclusion."
)

write.csv(df, file.path(out_dir, "step4_data_with_residuals.csv"), row.names = FALSE)
write.csv(hetero_tests, file.path(out_dir, "step4_heteroskedasticity_tests.csv"), row.names = FALSE)
write.csv(ols_standard_vs_hc, file.path(out_dir, "step4_ols_standard_vs_hc.csv"), row.names = FALSE)
write.csv(robust_wald_tests, file.path(out_dir, "step4_robust_wald_tests.csv"), row.names = FALSE)
write.csv(wls_fgls_comparison, file.path(out_dir, "step4_wls_fgls_comparison.csv"), row.names = FALSE)
write.csv(quartile_var, file.path(out_dir, "step4_residual_variance_by_quartile.csv"), row.names = FALSE)
write.csv(coefficient_table(variance_aux, "FGLS_variance_auxiliary", "const"), file.path(out_dir, "step4_fgls_variance_auxiliary.csv"), row.names = FALSE)
write.csv(coefficient_table(variance_aux_white, "FGLS_White_variance_auxiliary", "const"), file.path(out_dir, "step4_fgls_white_variance_auxiliary.csv"), row.names = FALSE)
influence_table <- df[order(df$cooks_d, decreasing = TRUE), c(
  "industry", "CAR", "RateVulnerability_z", "ols_resid", "hat", "cooks_d", "dfbetas_rate_vulnerability"
)]
write.csv(influence_table, file.path(out_dir, "step4_influence_diagnostics.csv"), row.names = FALSE)

cook_top3 <- head(influence_table$industry, 3)
influence_sensitivity <- do.call(
  rbind,
  list(
    fit_m2_hc3_rate(df, "M2_full_sample_HC3"),
    fit_m2_hc3_rate(df[df$industry != "Gold", ], "M2_excluding_Gold_HC3"),
    fit_m2_hc3_rate(df[!(df$industry %in% cook_top3), ], "M2_excluding_Cooks_top3_HC3")
  )
)
influence_sensitivity$excluded_industries <- c(
  "None",
  "Gold",
  paste(cook_top3, collapse = "; ")
)
write.csv(influence_sensitivity, file.path(out_dir, "step4_influence_sensitivity.csv"), row.names = FALSE)

leave_one_out_m2_hc3 <- do.call(
  rbind,
  lapply(df$industry, function(industry) {
    temp <- fit_m2_hc3_rate(df[df$industry != industry, ], paste0("drop_", industry))
    temp$dropped_industry <- industry
    temp
  })
)
leave_one_out_m2_hc3 <- leave_one_out_m2_hc3[, c(
  "dropped_industry", "model", "se_type", "term", "estimate", "std_error",
  "t_value", "p_value", "conf_low_95", "conf_high_95"
)]
write.csv(leave_one_out_m2_hc3, file.path(out_dir, "step4_leave_one_out_M2_HC3.csv"), row.names = FALSE)

writeLines(model_decision, file.path(out_dir, "step4_model_decision_notes.txt"))

png(file.path(out_dir, "step4_residuals_vs_fitted.png"), width = 1200, height = 850, res = 150)
plot(df$ols_fitted, df$ols_resid, pch = 19, col = "#1f77b4", xlab = "Fitted CAR, percent", ylab = "OLS residual", main = "Residuals vs fitted values")
abline(h = 0, col = "#d62728", lwd = 2)
grid(col = "gray85")
text(df$ols_fitted, df$ols_resid, labels = df$industry, pos = 3, cex = 0.55, col = "gray30")
dev.off()

png(file.path(out_dir, "step4_abs_residuals_vs_fitted.png"), width = 1200, height = 850, res = 150)
plot(df$ols_fitted, df$ols_abs_resid, pch = 19, col = "#1f77b4", xlab = "Fitted CAR, percent", ylab = "|OLS residual|", main = "Absolute residuals vs fitted values")
grid(col = "gray85")
text(df$ols_fitted, df$ols_abs_resid, labels = df$industry, pos = 3, cex = 0.55, col = "gray30")
dev.off()

png(file.path(out_dir, "step4_residuals_vs_rate_vulnerability.png"), width = 1200, height = 850, res = 150)
plot(df$RateVulnerability_z, df$ols_resid, pch = 19, col = "#1f77b4", xlab = "Rate vulnerability (standardized)", ylab = "OLS residual", main = "Residuals vs rate vulnerability")
abline(h = 0, col = "#d62728", lwd = 2)
grid(col = "gray85")
text(df$RateVulnerability_z, df$ols_resid, labels = df$industry, pos = 3, cex = 0.55, col = "gray30")
dev.off()

png(file.path(out_dir, "step4_residuals_by_rate_vulnerability_quartile.png"), width = 1200, height = 850, res = 150)
barplot(quartile_var$residual_variance, names.arg = quartile_var$rv_quartile, col = "#1f77b4", ylab = "Residual variance", main = "Residual variance by rate-vulnerability quartile")
grid(nx = NA, ny = NULL, col = "gray85")
dev.off()

png(file.path(out_dir, "step4_scale_location_plot.png"), width = 1200, height = 850, res = 150)
plot(df$ols_fitted, sqrt(abs(df$ols_resid)), pch = 19, col = "#1f77b4", xlab = "Fitted CAR, percent", ylab = "sqrt(|OLS residual|)", main = "Scale-location plot")
grid(col = "gray85")
text(df$ols_fitted, sqrt(abs(df$ols_resid)), labels = df$industry, pos = 3, cex = 0.55, col = "gray30")
dev.off()

cat("Done.\n")
cat("Output directory:", out_dir, "\n\n")
cat("Heteroskedasticity tests:\n")
print(hetero_tests, row.names = FALSE)
cat("\nOLS standard vs HC for RateVulnerability_z:\n")
print(ols_standard_vs_hc, row.names = FALSE)
cat("\nRobust Wald tests:\n")
print(robust_wald_tests, row.names = FALSE)
cat("\nWLS/FGLS comparison:\n")
print(wls_fgls_comparison, row.names = FALSE)
cat("\nInfluence sensitivity:\n")
print(influence_sensitivity, row.names = FALSE)
cat("\nLeave-one-out HC3 summary:\n")
print(
  data.frame(
    negative_coefficients = sum(leave_one_out_m2_hc3$estimate < 0),
    significant_at_5pct = sum(leave_one_out_m2_hc3$p_value < 0.05),
    significant_at_10pct = sum(leave_one_out_m2_hc3$p_value < 0.10),
    max_p_value = max(leave_one_out_m2_hc3$p_value),
    weakest_drop = leave_one_out_m2_hc3$dropped_industry[which.max(leave_one_out_m2_hc3$p_value)],
    total = nrow(leave_one_out_m2_hc3)
  ),
  row.names = FALSE
)
