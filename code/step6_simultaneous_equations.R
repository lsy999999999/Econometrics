options(stringsAsFactors = FALSE)

this_file <- tryCatch(normalizePath(sys.frames()[[1]]$ofile), error = function(e) NA_character_)
if (is.na(this_file)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) this_file <- normalizePath(sub("^--file=", "", file_arg[1]))
}

project_root <- normalizePath(file.path(dirname(this_file), ".."))
in_file <- file.path(project_root, "output", "step5", "step5_data.csv")
out_dir <- file.path(project_root, "output", "step6")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_file)) stop("Missing ", in_file, ". Run code/step5_endogeneity_iv.R first.")

df <- read.csv(in_file)
df$CAR_3day <- df$CAR

required_vars <- c(
  "industry", "CAR_3day", "RateVulnerability_z",
  "Z_precovid_2Y_z", "Z_early_nonoverlap_2Y_z", "Z_verylong_2Y_z", "Z_termspread_z",
  "MktBeta_z", "SMBBeta_z", "HMLBeta_z", "RMWBeta_z", "CMABeta_z",
  "HistVol_z", "PreMomentum_z"
)

missing_variables <- data.frame(
  variable = required_vars,
  available = required_vars %in% names(df),
  note = ifelse(required_vars %in% names(df), "available", "missing from step5_data.csv")
)
write.csv(missing_variables, file.path(out_dir, "step6_missing_variables.csv"), row.names = FALSE)
if (any(!missing_variables$available)) {
  stop("Step 6 missing required variables. See output/step6/step6_missing_variables.csv")
}

controls_return <- c("MktBeta_z", "SMBBeta_z", "HMLBeta_z", "RMWBeta_z", "CMABeta_z", "HistVol_z", "PreMomentum_z")
historical_ivs <- c("Z_precovid_2Y_z", "Z_early_nonoverlap_2Y_z", "Z_verylong_2Y_z")
historical_ivs_overid <- list(
  precovid_plus_verylong = c("Z_precovid_2Y_z", "Z_verylong_2Y_z"),
  precovid_plus_early_nonoverlap = c("Z_precovid_2Y_z", "Z_early_nonoverlap_2Y_z"),
  precovid_plus_termspread = c("Z_precovid_2Y_z", "Z_termspread_z")
)
controls_vulnerability <- historical_ivs
reverse_excluded <- controls_return

complete_vars <- unique(c("industry", "CAR_3day", "RateVulnerability_z", controls_return, historical_ivs, "Z_termspread_z"))
sem_data <- df[complete.cases(df[, complete_vars]), complete_vars]
write.csv(sem_data, file.path(out_dir, "step6_sem_data.csv"), row.names = FALSE)

hc_vcov <- function(model, type = "HC3") {
  x <- model.matrix(model)
  residual <- residuals(model)
  bread <- solve(crossprod(x))
  hat_values <- hatvalues(model)
  adjustment <- switch(
    type,
    HC0 = rep(1, length(residual)),
    HC1 = rep(nrow(x) / (nrow(x) - ncol(x)), length(residual)),
    HC2 = 1 / pmax(1 - hat_values, .Machine$double.eps),
    HC3 = 1 / pmax(1 - hat_values, .Machine$double.eps)^2,
    stop("Unsupported vcov type")
  )
  meat <- t(x) %*% ((residual^2 * adjustment) * x)
  bread %*% meat %*% bread
}

coef_row_lm <- function(model, equation, method, key_term, note = "") {
  vc <- hc_vcov(model, "HC3")
  estimates <- coef(model)
  se <- sqrt(diag(vc))
  t_value <- estimates / se
  p_value <- 2 * pt(abs(t_value), df = df.residual(model), lower.tail = FALSE)
  idx <- match(key_term, names(estimates))
  data.frame(
    equation = equation,
    method = method,
    key_term = key_term,
    estimate = unname(estimates[idx]),
    robust_se = unname(se[idx]),
    statistic = unname(t_value[idx]),
    p_value = unname(p_value[idx]),
    N = nobs(model),
    note = note,
    row.names = NULL
  )
}

manual_2sls <- function(data, y, x_endog, controls, instruments) {
  rhs_x <- c(x_endog, controls)
  rhs_z <- c(instruments, controls)
  y_vec <- as.matrix(data[, y, drop = FALSE])
  x_matrix <- model.matrix(as.formula(paste("~", paste(rhs_x, collapse = " + "))), data = data)
  z_matrix <- model.matrix(as.formula(paste("~", paste(rhs_z, collapse = " + "))), data = data)
  pz <- z_matrix %*% MASS::ginv(crossprod(z_matrix)) %*% t(z_matrix)
  x_hat <- pz %*% x_matrix
  bread <- MASS::ginv(t(x_hat) %*% x_matrix)
  beta <- bread %*% t(x_hat) %*% y_vec
  rownames(beta) <- colnames(x_matrix)
  fitted_values <- x_matrix %*% beta
  residual <- as.vector(y_vec - fitted_values)
  n <- nrow(x_matrix)
  k <- ncol(x_matrix)
  hat_diag <- diag(x_hat %*% bread %*% t(x_hat))
  adjustment <- 1 / pmax(1 - hat_diag, .Machine$double.eps)^2
  meat <- t(x_hat) %*% ((residual^2 * adjustment) * x_hat)
  vcov_hc3 <- bread %*% meat %*% t(bread)
  se <- sqrt(diag(vcov_hc3))
  t_value <- as.vector(beta) / se
  p_value <- 2 * pt(abs(t_value), df = n - k, lower.tail = FALSE)
  list(
    coefficients = as.vector(beta),
    names = rownames(beta),
    vcov_hc3 = vcov_hc3,
    se_hc3 = se,
    t_value = t_value,
    p_value = p_value,
    residuals = residual,
    fitted = as.vector(fitted_values),
    x_matrix = x_matrix,
    z_matrix = z_matrix,
    n = n,
    k = k,
    instruments = instruments,
    controls = controls,
    x_endog = x_endog,
    y = y
  )
}

coef_row_iv <- function(iv_fit, equation, method, note = "") {
  idx <- match(iv_fit$x_endog, iv_fit$names)
  data.frame(
    equation = equation,
    method = method,
    key_term = iv_fit$x_endog,
    estimate = iv_fit$coefficients[idx],
    robust_se = iv_fit$se_hc3[idx],
    statistic = iv_fit$t_value[idx],
    p_value = iv_fit$p_value[idx],
    N = iv_fit$n,
    note = note,
    row.names = NULL
  )
}

first_stage_table <- function(data, equation, x_endog, included_exog, excluded_exog) {
  full_rhs <- c(excluded_exog, included_exog)
  full <- lm(as.formula(paste(x_endog, "~", paste(full_rhs, collapse = " + "))), data = data)
  if (length(included_exog) > 0) {
    restricted <- lm(as.formula(paste(x_endog, "~", paste(included_exog, collapse = " + "))), data = data)
  } else {
    restricted <- lm(as.formula(paste(x_endog, "~ 1")), data = data)
  }
  a <- anova(restricted, full)
  partial_r2 <- (sum(resid(restricted)^2) - sum(resid(full)^2)) / sum(resid(restricted)^2)
  coef_mat <- summary(full)$coefficients
  rows <- data.frame(
    equation = equation,
    endogenous_rhs_variable = x_endog,
    excluded_IV = excluded_exog,
    coefficient = unname(coef(full)[excluded_exog]),
    std_error = coef_mat[excluded_exog, "Std. Error"],
    p_value = coef_mat[excluded_exog, "Pr(>|t|)"],
    joint_first_stage_F = a$F[2],
    joint_first_stage_p = a$`Pr(>F)`[2],
    partial_R2 = partial_r2,
    weak_IV_flag_F_below_10 = isTRUE(a$F[2] < 10),
    N = nobs(full),
    row.names = NULL
  )
  rows
}

dwh_test <- function(data, equation, y, x_endog, controls, instruments) {
  first <- lm(as.formula(paste(x_endog, "~", paste(c(instruments, controls), collapse = " + "))), data = data)
  data$first_stage_resid <- resid(first)
  augmented <- lm(as.formula(paste(y, "~", paste(c(x_endog, controls, "first_stage_resid"), collapse = " + "))), data = data)
  coef_mat <- summary(augmented)$coefficients
  data.frame(
    equation = equation,
    test = "Durbin-Wu-Hausman residual inclusion",
    endogenous_rhs_variable = x_endog,
    instruments = paste(instruments, collapse = "; "),
    residual_coefficient = coef(augmented)["first_stage_resid"],
    p_value = coef_mat["first_stage_resid", "Pr(>|t|)"],
    reject_exogeneity_5pct = coef_mat["first_stage_resid", "Pr(>|t|)"] < 0.05,
    N = nobs(augmented),
    row.names = NULL
  )
}

sargan_test <- function(iv_fit, equation, label) {
  overid_df <- ncol(iv_fit$z_matrix) - ncol(iv_fit$x_matrix)
  if (overid_df <= 0) {
    return(data.frame(
      equation = equation,
      IV_set = label,
      test_type = "Sargan",
      statistic = NA_real_,
      df = overid_df,
      p_value = NA_real_,
      reject_at_5pct = NA,
      interpretation = "Exactly identified; overidentification test unavailable.",
      row.names = NULL
    ))
  }
  aux <- lm(iv_fit$residuals ~ iv_fit$z_matrix[, -1, drop = FALSE])
  statistic <- iv_fit$n * summary(aux)$r.squared
  data.frame(
    equation = equation,
    IV_set = label,
    test_type = "Sargan",
    statistic = statistic,
    df = overid_df,
    p_value = pchisq(statistic, df = overid_df, lower.tail = FALSE),
    reject_at_5pct = pchisq(statistic, df = overid_df, lower.tail = FALSE) < 0.05,
    interpretation = "Tests overidentifying restrictions under homoskedasticity; non-rejection is not proof of exogeneity.",
    row.names = NULL
  )
}

hansen_j_test <- function(iv_fit, equation, label) {
  overid_df <- ncol(iv_fit$z_matrix) - ncol(iv_fit$x_matrix)
  if (overid_df <= 0) {
    return(data.frame(
      equation = equation,
      IV_set = label,
      test_type = "Hansen_style_J",
      statistic = NA_real_,
      df = overid_df,
      p_value = NA_real_,
      reject_at_5pct = NA,
      interpretation = "Exactly identified; robust overidentification test unavailable.",
      row.names = NULL
    ))
  }
  z <- iv_fit$z_matrix
  u <- iv_fit$residuals
  n <- nrow(z)
  moments <- z * as.numeric(u)
  gbar <- colMeans(moments)
  s_matrix <- crossprod(scale(moments, center = FALSE, scale = FALSE)) / n
  statistic <- as.numeric(n * t(gbar) %*% MASS::ginv(s_matrix) %*% gbar)
  p_value <- pchisq(statistic, df = overid_df, lower.tail = FALSE)
  data.frame(
    equation = equation,
    IV_set = label,
    test_type = "Hansen_style_J",
    statistic = statistic,
    df = overid_df,
    p_value = p_value,
    reject_at_5pct = p_value < 0.05,
    interpretation = "Robust moment-based J using 2SLS residuals; interpret cautiously with N = 49.",
    row.names = NULL
  )
}

identify_equation <- function(equation, dependent_endogenous, rhs_endogenous, included_exog, excluded_exog, economic_exclusion_plausible) {
  rhs_count <- length(rhs_endogenous)
  excluded_count <- length(excluded_exog)
  order_condition <- excluded_count >= rhs_count
  status <- if (!order_condition) {
    "underidentified_by_order_condition"
  } else if (excluded_count == rhs_count) {
    "just_identified_by_order_condition"
  } else {
    "overidentified_by_order_condition"
  }
  data.frame(
    equation = equation,
    dependent_endogenous = dependent_endogenous,
    rhs_endogenous = paste(rhs_endogenous, collapse = "; "),
    rhs_endogenous_count = rhs_count,
    included_exogenous = paste(included_exog, collapse = "; "),
    excluded_exogenous = paste(excluded_exog, collapse = "; "),
    excluded_exogenous_count = excluded_count,
    order_condition_satisfied = order_condition,
    identification_status_order_condition = status,
    economic_exclusion_plausible = economic_exclusion_plausible,
    note = ifelse(
      economic_exclusion_plausible,
      "The exclusion restriction is economically defensible enough for auxiliary IV evidence.",
      "The order condition may hold mechanically, but the exclusion restriction is economically contestable."
    ),
    row.names = NULL
  )
}

variable_roles <- data.frame(
  variable = c("CAR_3day", "RateVulnerability_z", historical_ivs, "Z_termspread_z", controls_return),
  role = c(
    "endogenous variable",
    "endogenous or predetermined exposure variable depending on specification",
    rep("excluded exogenous / predetermined IV for return equation", length(historical_ivs)),
    "excluded exogenous / predetermined IV candidate for return equation",
    rep("included controls in return equation; questionable excluded IVs in reverse equation", length(controls_return))
  ),
  equation = c(
    "return equation; exploratory reverse equation",
    "return equation; exploratory reverse equation",
    rep("excluded from return equation", length(historical_ivs)),
    "excluded from return equation in overidentification check",
    rep("included in return equation; excluded from exploratory vulnerability equation", length(controls_return))
  ),
  lecture6_concept = c(
    "endogenous variable",
    "endogenous variable in SEM; also a pre-event constructed regressor in the main empirical design",
    rep("exogenous / predetermined variable", length(historical_ivs) + 1),
    rep("exogenous controls; exclusion restriction must be justified if used as IVs", length(controls_return))
  ),
  row.names = NULL
)
write.csv(variable_roles, file.path(out_dir, "step6_sem_variable_roles.csv"), row.names = FALSE)

identification_conditions <- rbind(
  identify_equation(
    "return_equation",
    "CAR_3day",
    "RateVulnerability_z",
    controls_return,
    historical_ivs,
    TRUE
  ),
  identify_equation(
    "vulnerability_equation_exploratory",
    "RateVulnerability_z",
    "CAR_3day",
    controls_vulnerability,
    reverse_excluded,
    FALSE
  )
)
write.csv(identification_conditions, file.path(out_dir, "step6_identification_conditions.csv"), row.names = FALSE)

first_stage_relevance <- rbind(
  first_stage_table(sem_data, "return_equation_system_IV_set", "RateVulnerability_z", controls_return, historical_ivs),
  first_stage_table(sem_data, "return_equation_primary_IV", "RateVulnerability_z", controls_return, "Z_precovid_2Y_z"),
  first_stage_table(sem_data, "vulnerability_equation_exploratory", "CAR_3day", controls_vulnerability, reverse_excluded)
)
write.csv(first_stage_relevance, file.path(out_dir, "step6_first_stage_relevance.csv"), row.names = FALSE)
write.csv(first_stage_relevance, file.path(out_dir, "step6_weak_iv_tests.csv"), row.names = FALSE)

ols_return <- lm(CAR_3day ~ RateVulnerability_z + MktBeta_z + SMBBeta_z + HMLBeta_z + RMWBeta_z + CMABeta_z + HistVol_z + PreMomentum_z, data = sem_data)
ols_vulnerability <- lm(RateVulnerability_z ~ CAR_3day + Z_precovid_2Y_z + Z_early_nonoverlap_2Y_z + Z_verylong_2Y_z, data = sem_data)

iv_return_primary <- manual_2sls(sem_data, "CAR_3day", "RateVulnerability_z", controls_return, "Z_precovid_2Y_z")
iv_return_system_set <- manual_2sls(sem_data, "CAR_3day", "RateVulnerability_z", controls_return, historical_ivs)
iv_vulnerability <- manual_2sls(sem_data, "RateVulnerability_z", "CAR_3day", controls_vulnerability, reverse_excluded)

ols_vs_2sls_return <- rbind(
  coef_row_lm(ols_return, "return_equation", "OLS_HC3", "RateVulnerability_z", "Baseline single-equation estimate."),
  coef_row_iv(iv_return_primary, "return_equation", "2SLS_manual_HC3_primary_Z_precovid", "Primary IV uses the strongest predetermined historical exposure."),
  coef_row_iv(iv_return_system_set, "return_equation", "2SLS_manual_HC3_system_IV_set", "Uses the SEM historical IV set.")
)
write.csv(ols_vs_2sls_return, file.path(out_dir, "step6_ols_vs_2sls_return_equation.csv"), row.names = FALSE)

ols_vs_2sls_vulnerability <- rbind(
  coef_row_lm(ols_vulnerability, "vulnerability_equation_exploratory", "OLS_HC3", "CAR_3day", "Reverse equation is exploratory because CAR occurs after the historical RV measure."),
  coef_row_iv(iv_vulnerability, "vulnerability_equation_exploratory", "2SLS_manual_HC3", "Technically estimable, but exclusion of factor betas from the RV equation is economically contestable.")
)
write.csv(ols_vs_2sls_vulnerability, file.path(out_dir, "step6_ols_vs_2sls_vulnerability_equation.csv"), row.names = FALSE)

overid_tests <- do.call(rbind, lapply(names(historical_ivs_overid), function(label) {
  fit <- manual_2sls(sem_data, "CAR_3day", "RateVulnerability_z", controls_return, historical_ivs_overid[[label]])
  rbind(
    sargan_test(fit, "return_equation", label),
    hansen_j_test(fit, "return_equation", label)
  )
}))
overid_reverse <- rbind(
  sargan_test(iv_vulnerability, "vulnerability_equation_exploratory", "factor_betas_as_excluded_IVs"),
  hansen_j_test(iv_vulnerability, "vulnerability_equation_exploratory", "factor_betas_as_excluded_IVs")
)
overid_tests <- rbind(overid_tests, overid_reverse)
write.csv(overid_tests, file.path(out_dir, "step6_overid_tests.csv"), row.names = FALSE)

endogeneity_tests <- rbind(
  dwh_test(sem_data, "return_equation_primary_IV", "CAR_3day", "RateVulnerability_z", controls_return, "Z_precovid_2Y_z"),
  dwh_test(sem_data, "return_equation_system_IV_set", "CAR_3day", "RateVulnerability_z", controls_return, historical_ivs),
  dwh_test(sem_data, "vulnerability_equation_exploratory", "RateVulnerability_z", "CAR_3day", controls_vulnerability, reverse_excluded)
)
write.csv(endogeneity_tests, file.path(out_dir, "step6_endogeneity_tests.csv"), row.names = FALSE)

simple_car_rf <- lm(CAR_3day ~ Z_precovid_2Y_z, data = sem_data)
simple_rv_rf <- lm(RateVulnerability_z ~ Z_precovid_2Y_z, data = sem_data)
simple_iv_no_controls <- manual_2sls(sem_data, "CAR_3day", "RateVulnerability_z", character(), "Z_precovid_2Y_z")
ils_demo <- data.frame(
  item = c("reduced_form_CAR_on_Z", "reduced_form_RV_on_Z", "ILS_Wald_ratio_no_controls", "2SLS_no_controls_same_IV"),
  estimate = c(
    coef(simple_car_rf)["Z_precovid_2Y_z"],
    coef(simple_rv_rf)["Z_precovid_2Y_z"],
    coef(simple_car_rf)["Z_precovid_2Y_z"] / coef(simple_rv_rf)["Z_precovid_2Y_z"],
    simple_iv_no_controls$coefficients[match("RateVulnerability_z", simple_iv_no_controls$names)]
  ),
  note = c(
    "Reduced form: external predetermined variable to endogenous CAR.",
    "Reduced form: external predetermined variable to endogenous RV.",
    "Teaching-only ILS/Wald ratio in a just-identified no-control simplification.",
    "Matches the just-identified IV estimate in the same no-control simplification."
  ),
  row.names = NULL
)
write.csv(ils_demo, file.path(out_dir, "step6_ils_reduced_form_demo.csv"), row.names = FALSE)

system_results <- data.frame()
system_residual_correlation <- data.frame(
  method = character(),
  residual_correlation = numeric(),
  p_value = numeric(),
  note = character()
)

system_summary_path <- file.path(out_dir, "step6_system_2sls_summary.txt")
system3_summary_path <- file.path(out_dir, "step6_system_3sls_summary.txt")

systemfit_available <- requireNamespace("systemfit", quietly = TRUE)
if (systemfit_available) {
  eq_return <- CAR_3day ~ RateVulnerability_z + MktBeta_z + SMBBeta_z + HMLBeta_z + RMWBeta_z + CMABeta_z + HistVol_z + PreMomentum_z
  eq_vulnerability <- RateVulnerability_z ~ CAR_3day + Z_precovid_2Y_z + Z_early_nonoverlap_2Y_z + Z_verylong_2Y_z
  eqs <- list(returnEq = eq_return, vulnEq = eq_vulnerability)
  inst_formula <- ~ Z_precovid_2Y_z + Z_early_nonoverlap_2Y_z + Z_verylong_2Y_z +
    MktBeta_z + SMBBeta_z + HMLBeta_z + RMWBeta_z + CMABeta_z + HistVol_z + PreMomentum_z
  fit_sys_2sls <- tryCatch(systemfit::systemfit(eqs, method = "2SLS", inst = inst_formula, data = sem_data), error = identity)
  if (!inherits(fit_sys_2sls, "error")) {
    writeLines(capture.output(summary(fit_sys_2sls)), system_summary_path)
    b <- coef(fit_sys_2sls)
    v <- tryCatch(vcov(fit_sys_2sls), error = function(e) NULL)
    se <- if (!is.null(v)) sqrt(diag(v)) else rep(NA_real_, length(b))
    system_results <- rbind(system_results, data.frame(
      method = "systemfit_2SLS",
      term = names(b),
      estimate = unname(b),
      standard_error = unname(se),
      statistic = unname(b / se),
      p_value = 2 * pnorm(abs(b / se), lower.tail = FALSE),
      note = "System estimate; reverse equation remains exploratory.",
      row.names = NULL
    ))
    residual_matrix <- residuals(fit_sys_2sls)
    if (is.matrix(residual_matrix) || is.data.frame(residual_matrix)) {
      ct <- cor.test(residual_matrix[, 1], residual_matrix[, 2])
      system_residual_correlation <- rbind(system_residual_correlation, data.frame(
        method = "systemfit_2SLS",
        residual_correlation = unname(ct$estimate),
        p_value = ct$p.value,
        note = "Correlation of structural residual estimates; 3SLS is mainly an efficiency refinement when errors are correlated.",
        row.names = NULL
      ))
    }
  } else {
    writeLines(paste("systemfit 2SLS failed:", fit_sys_2sls$message), system_summary_path)
  }
  fit_sys_3sls <- tryCatch(systemfit::systemfit(eqs, method = "3SLS", inst = inst_formula, data = sem_data), error = identity)
  if (!inherits(fit_sys_3sls, "error")) {
    writeLines(capture.output(summary(fit_sys_3sls)), system3_summary_path)
    b <- coef(fit_sys_3sls)
    v <- tryCatch(vcov(fit_sys_3sls), error = function(e) NULL)
    se <- if (!is.null(v)) sqrt(diag(v)) else rep(NA_real_, length(b))
    system_results <- rbind(system_results, data.frame(
      method = "systemfit_3SLS",
      term = names(b),
      estimate = unname(b),
      standard_error = unname(se),
      statistic = unname(b / se),
      p_value = 2 * pnorm(abs(b / se), lower.tail = FALSE),
      note = "Optional system estimator; not used as main evidence.",
      row.names = NULL
    ))
  } else {
    writeLines(paste("systemfit 3SLS failed:", fit_sys_3sls$message), system3_summary_path)
  }
} else {
  writeLines("systemfit package unavailable; system 2SLS/3SLS not run.", system_summary_path)
  writeLines("systemfit package unavailable; system 2SLS/3SLS not run.", system3_summary_path)
}

if (nrow(system_results) == 0) {
  system_results <- data.frame(
    method = character(),
    term = character(),
    estimate = numeric(),
    standard_error = numeric(),
    statistic = numeric(),
    p_value = numeric(),
    note = character()
  )
}
write.csv(system_results, file.path(out_dir, "step6_system_results.csv"), row.names = FALSE)
write.csv(system_results[system_results$method == "systemfit_2SLS", ], file.path(out_dir, "step6_system_2sls_results.csv"), row.names = FALSE)
write.csv(system_results[system_results$method == "systemfit_3SLS", ], file.path(out_dir, "step6_system_3sls_results.csv"), row.names = FALSE)

if (nrow(system_residual_correlation) == 0) {
  ols_resid_ct <- cor.test(resid(ols_return), resid(ols_vulnerability))
  system_residual_correlation <- data.frame(
    method = "OLS_equation_residuals_fallback",
    residual_correlation = unname(ols_resid_ct$estimate),
    p_value = ols_resid_ct$p.value,
    note = "Fallback residual correlation because systemfit residual matrix was unavailable.",
    row.names = NULL
  )
}
write.csv(system_residual_correlation, file.path(out_dir, "step6_system_residual_correlation.csv"), row.names = FALSE)

system_key_rows <- data.frame()
if (nrow(system_results) > 0) {
  key_terms <- c("returnEq_RateVulnerability_z", "vulnEq_CAR_3day")
  available_terms <- intersect(key_terms, system_results$term)
  if (length(available_terms) > 0) {
    system_key_rows <- data.frame(
      equation = ifelse(
        system_results$term[system_results$term %in% available_terms] == "returnEq_RateVulnerability_z",
        "return_equation",
        "vulnerability_equation_exploratory"
      ),
      method = system_results$method[system_results$term %in% available_terms],
      key_term = ifelse(
        system_results$term[system_results$term %in% available_terms] == "returnEq_RateVulnerability_z",
        "RateVulnerability_z",
        "CAR_3day"
      ),
      estimate = system_results$estimate[system_results$term %in% available_terms],
      robust_se = system_results$standard_error[system_results$term %in% available_terms],
      statistic = system_results$statistic[system_results$term %in% available_terms],
      p_value = system_results$p_value[system_results$term %in% available_terms],
      N = nrow(sem_data),
      note = "Systemfit conventional SE; not HC3.",
      row.names = NULL
    )
  }
}
model_comparison <- rbind(ols_vs_2sls_return, ols_vs_2sls_vulnerability, system_key_rows)
model_comparison$sign <- ifelse(model_comparison$estimate > 0, "positive", ifelse(model_comparison$estimate < 0, "negative", "zero"))
model_comparison$conclusion <- ifelse(
  grepl("vulnerability_equation", model_comparison$equation),
  "Exploratory only; reverse structure is not credible enough for main identification.",
  ifelse(model_comparison$p_value < 0.05, "Statistically significant at 5%.", "Not statistically significant at 5%.")
)
write.csv(model_comparison, file.path(out_dir, "step6_model_comparison.csv"), row.names = FALSE)

recursive_ct <- cor.test(resid(ols_return), resid(ols_vulnerability))
recursive_check <- data.frame(
  condition = c(
    "one_way_time_order_RV_before_CAR",
    "contemporaneous_reverse_CAR_determines_historical_RV",
    "structural_error_independence_required_for_recursive_OLS",
    "residual_correlation_proxy"
  ),
  result = c(TRUE, FALSE, recursive_ct$p.value >= 0.05, abs(unname(recursive_ct$estimate))),
  note = c(
    "The main empirical design constructs RV from pre-event historical windows.",
    "Same-event CAR cannot literally determine a historical RV estimate.",
    "A recursive system needs one-way ordering and uncorrelated structural errors.",
    paste0("OLS residual correlation = ", round(unname(recursive_ct$estimate), 3), ", p = ", round(recursive_ct$p.value, 3), ".")
  ),
  row.names = NULL
)
write.csv(recursive_check, file.path(out_dir, "step6_recursive_check.csv"), row.names = FALSE)

return_primary_fs <- first_stage_relevance[first_stage_relevance$equation == "return_equation_primary_IV", ][1, ]
return_system_fs <- first_stage_relevance[first_stage_relevance$equation == "return_equation_system_IV_set", ][1, ]
reverse_fs <- first_stage_relevance[first_stage_relevance$equation == "vulnerability_equation_exploratory", ][1, ]
return_overid_rejected <- any(overid_tests$equation == "return_equation" & overid_tests$reject_at_5pct %in% TRUE, na.rm = TRUE)
reverse_overid_rejected <- any(overid_tests$equation == "vulnerability_equation_exploratory" & overid_tests$reject_at_5pct %in% TRUE, na.rm = TRUE)
overid_rejected <- any(overid_tests$reject_at_5pct %in% TRUE, na.rm = TRUE)
system_resid_cor_p <- system_residual_correlation$p_value[1]
use_sem_as_main <- FALSE

sem_summary_flags <- data.frame(
  item = c(
    "return equation order condition satisfied",
    "vulnerability equation order condition satisfied",
    "vulnerability equation economic exclusion plausible",
    "primary return first stage strong F>=10",
    "system return first stage strong F>=10",
    "exploratory reverse first stage strong F>=10",
    "return equation overidentification rejected at 5pct",
    "exploratory reverse overidentification rejected at 5pct",
    "overidentification rejected anywhere at 5pct",
    "system residuals correlated at 5pct",
    "use SEM as main result"
  ),
  result = c(
    identification_conditions$order_condition_satisfied[identification_conditions$equation == "return_equation"],
    identification_conditions$order_condition_satisfied[identification_conditions$equation == "vulnerability_equation_exploratory"],
    identification_conditions$economic_exclusion_plausible[identification_conditions$equation == "vulnerability_equation_exploratory"],
    return_primary_fs$joint_first_stage_F >= 10,
    return_system_fs$joint_first_stage_F >= 10,
    reverse_fs$joint_first_stage_F >= 10,
    return_overid_rejected,
    reverse_overid_rejected,
    overid_rejected,
    system_resid_cor_p < 0.05,
    use_sem_as_main
  ),
  value = c(
    "",
    "",
    "",
    round(return_primary_fs$joint_first_stage_F, 3),
    round(return_system_fs$joint_first_stage_F, 3),
    round(reverse_fs$joint_first_stage_F, 3),
    "",
    "",
    "",
    round(system_resid_cor_p, 3),
    ""
  ),
  interpretation = c(
    "Return equation has enough excluded predetermined variables by the order condition.",
    "Reverse equation is mechanically identifiable by the order condition.",
    "False: factor betas and controls plausibly affect RV directly, so they are weak excluded IVs for CAR.",
    "The strongest predetermined IV is relevant enough for auxiliary IV evidence.",
    "The larger historical IV set is only moderate because the instruments are correlated.",
    "Even if relevant, the reverse equation fails the economic timing/exclusion logic.",
    "Return-equation overidentification tests do not reject in the tested IV sets.",
    "The exploratory reverse equation has at least one rejection, reinforcing that it should not be main evidence.",
    "Non-rejection is not proof of validity; rejection is a warning.",
    "A low p-value would motivate 3SLS as an efficiency check, not a fix for identification.",
    "False: the complete two-way SEM is exploratory; the main result should be the return equation plus robustness checks."
  ),
  row.names = NULL
)
write.csv(sem_summary_flags, file.path(out_dir, "step6_sem_summary_flags.csv"), row.names = FALSE)

fit_warning <- data.frame(
  item = c("2SLS_R2", "2SLS_ordinary_F", "recommended_inference"),
  lecture6_point = c(
    "The ordinary R-squared from a 2SLS structural equation can be negative.",
    "The usual OLS-style overall F statistic is not generally valid for 2SLS.",
    "Use coefficient tests, IV first-stage diagnostics, DWH, overidentification tests, and Wald-style restrictions."
  ),
  implementation = c(
    "Step 6 does not use ordinary 2SLS R-squared for conclusions.",
    "Step 6 does not use ordinary OLS F statistics for 2SLS model fit.",
    "Model comparison reports key coefficients and diagnostic files."
  ),
  row.names = NULL
)
write.csv(fit_warning, file.path(out_dir, "step6_2sls_fit_warning.csv"), row.names = FALSE)

fmt <- function(x, digits = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
get_row <- function(tab, equation, method_pattern) {
  tab[tab$equation == equation & grepl(method_pattern, tab$method), ][1, ]
}
ols_ret_row <- get_row(model_comparison, "return_equation", "^OLS")
iv_ret_primary_row <- get_row(model_comparison, "return_equation", "primary_Z_precovid")
iv_ret_system_row <- get_row(model_comparison, "return_equation", "system_IV_set")
ols_rev_row <- get_row(model_comparison, "vulnerability_equation_exploratory", "^OLS")
iv_rev_row <- get_row(model_comparison, "vulnerability_equation_exploratory", "2SLS_manual")

coverage_lines <- c(
  "# Step 6 lecture coverage check",
  "",
  "This module maps Lecture 6, simultaneous equations models, into the empirical project while keeping the economic identification limits explicit.",
  "",
  "Covered topics:",
  "1. Simultaneous-equation motivation: `CAR_3day` and `RateVulnerability_z` are placed in a two-equation exploratory system.",
  "2. Structural equations: return equation and exploratory vulnerability equation are both estimated.",
  "3. Structural coefficients: key coefficients are saved in `step6_model_comparison.csv`.",
  "4. Structural errors: residual correlation is saved in `step6_system_residual_correlation.csv`.",
  "5. Endogenous variables: `CAR_3day` and `RateVulnerability_z` are listed in `step6_sem_variable_roles.csv`.",
  "6. Exogenous and predetermined variables: historical rate-exposure IVs and factor controls are classified in the variable-role table.",
  "7. Predetermined variables: pre-event windows such as `Z_precovid_2Y_z` are treated as predetermined external shifters.",
  "8. OLS inconsistency under simultaneity: OLS is reported only as a baseline comparison.",
  "9. Recursive-system logic: `step6_recursive_check.csv` notes that a valid recursive system requires one-way ordering and uncorrelated errors.",
  "10. ILS: `step6_ils_reduced_form_demo.csv` gives a just-identified reduced-form/Wald-ratio teaching example.",
  "11. 2SLS: manual HC3 2SLS results are saved for both equations.",
  "12. System estimation: `systemfit` 2SLS and optional 3SLS outputs are saved when available.",
  "13. Order condition: `step6_identification_conditions.csv` classifies under/just/over identification by excluded exogenous variables.",
  "14. Rank-condition intuition: `step6_first_stage_relevance.csv` reports first-stage F and partial R-squared as empirical relevance checks.",
  "15. Overidentification: Sargan and Hansen-style J tests are saved in `step6_overid_tests.csv`.",
  "16. 2SLS fit-statistic warning: `step6_2sls_fit_warning.csv` records the Lecture 6 warning that ordinary R-squared/F are not reliable for 2SLS.",
  "",
  "Important limitation:",
  "The reverse equation is not used as main evidence because same-event CAR cannot literally determine a historical pre-event RV measure, and the excluded controls for CAR plausibly affect RV directly."
)
writeLines(coverage_lines, file.path(out_dir, "step6_lecture_coverage_check.md"))

result_lines <- c(
  "# 第六讲联立方程模块结果",
  "",
  "## 总结判断",
  "",
  "第六讲要求关注变量相互决定、联立导致右侧内生变量、单方程 OLS 通常不一致。本文数据可以做一个完整的探索性两方程系统，但不能把它作为主识别结果。更稳妥的结论是：主方程仍应是 `RateVulnerability_z -> CAR_3day` 的收益反应方程；完整 `CAR_3day <-> RateVulnerability_z` 系统只作为第六讲概念演示和稳健性探索。",
  "",
  "原因是 `RateVulnerability_z` 来自事件前历史窗口，而 `CAR_3day` 是 2022-11-02 FOMC 事件窗口的市场反应。同一事件窗口的 CAR 不能反过来决定历史估计得到的 RV。因此，反向方程虽然机械上可估计，但经济解释和排除限制都不够强。",
  "",
  "## 识别条件",
  "",
  paste0("- 收益方程：阶条件满足，`RateVulnerability_z` 是右侧内生变量，排除的前定 IV 包括 `", paste(historical_ivs, collapse = "`, `"), "`。"),
  paste0("- 反向脆弱度方程：阶条件机械上满足，但 `", paste(reverse_excluded, collapse = "`, `"), "` 很可能直接影响 `RateVulnerability_z`，所以排除限制不可信。"),
  "",
  "## 第一阶段相关性",
  "",
  paste0("- 收益方程主 IV `Z_precovid_2Y_z`：第一阶段 F = ", fmt(return_primary_fs$joint_first_stage_F), "，partial R2 = ", fmt(return_primary_fs$partial_R2), "。"),
  paste0("- 收益方程系统 IV 组：第一阶段 F = ", fmt(return_system_fs$joint_first_stage_F), "，partial R2 = ", fmt(return_system_fs$partial_R2), "。"),
  paste0("- 探索性反向方程：第一阶段 F = ", fmt(reverse_fs$joint_first_stage_F), "，partial R2 = ", fmt(reverse_fs$partial_R2), "；即使相关性存在，经济排除限制仍然不足。"),
  "",
  "## OLS 与 2SLS 对比",
  "",
  paste0("- 收益方程 OLS-HC3：`RateVulnerability_z` 系数 = ", fmt(ols_ret_row$estimate), "，p = ", fmt(ols_ret_row$p_value), "。"),
  paste0("- 收益方程 2SLS，主 IV `Z_precovid_2Y_z`：系数 = ", fmt(iv_ret_primary_row$estimate), "，p = ", fmt(iv_ret_primary_row$p_value), "。"),
  paste0("- 收益方程 2SLS，系统 IV 组：系数 = ", fmt(iv_ret_system_row$estimate), "，p = ", fmt(iv_ret_system_row$p_value), "。"),
  paste0("- 反向方程 OLS-HC3：`CAR_3day` 系数 = ", fmt(ols_rev_row$estimate), "，p = ", fmt(ols_rev_row$p_value), "。"),
  paste0("- 反向方程 2SLS：`CAR_3day` 系数 = ", fmt(iv_rev_row$estimate), "，p = ", fmt(iv_rev_row$p_value), "；该结果只作探索性展示。"),
  "",
  "## 诊断结论",
  "",
  paste0("- 收益方程过度识别检验是否在 5% 水平拒绝：", return_overid_rejected, "。不拒绝不能证明工具变量外生，只能说明没有发现明显反证。"),
  paste0("- 探索性反向方程过度识别检验是否在 5% 水平拒绝：", reverse_overid_rejected, "。这个警告进一步说明反向方程不应作为主结果。"),
  paste0("- 系统残差相关检验 p = ", fmt(system_resid_cor_p), "。3SLS 只作为效率改进检查，不解决反向方程的经济识别问题。"),
  "",
  "## 最终写法",
  "",
  "论文中不要写“完整双向 SEM 证明了 CAR 和 RV 相互决定”。建议写：",
  "",
  "> 按第六讲联立方程模型的逻辑，本文进一步构造了收益反应方程和脆弱度方程组成的探索性两方程系统，并检查阶条件、工具变量相关性、内生性和过度识别限制。结果显示，收益方程具有较清晰的前定工具变量和时间顺序；但反向方程依赖的排除限制和经济含义较弱，因此完整 SEM 不作为主识别策略。本文主结论仍基于收益方程的 OLS-HC3、IV/2SLS 辅助证据和前几步稳健性检验。",
  "",
  "所有结果表保存在 `output/step6/`。"
)
writeLines(result_lines, file.path(out_dir, "step6_results_summary.md"))

message("Step 6 complete. Outputs written to ", out_dir)


write.csv(variable_roles, file.path(out_dir, "step6_sem_variable_roles.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Step 6: Lecture 6 concept definitions and empirical mapping
# ------------------------------------------------------------

concept_definitions <- data.frame(
  concept = c(
    "结构方程 structural equation",
    "结构系数 structural coefficient",
    "结构误差 structural error",
    "内生变量 endogenous variable",
    "外生变量 exogenous variable",
    "前定变量 predetermined variable",
    "右侧内生变量 endogenous RHS variable",
    "被排除外生变量 excluded exogenous variable",
    "约简型方程 reduced-form equation",
    "探索性联立系统 exploratory SEM"
  ),
  lecture6_definition = c(
    "联立方程模型中的单个方程，每个方程应具有其他条件不变下的合理经济解释。",
    "结构方程中的参数，反映系统内部变量之间的结构性关系。",
    "结构方程中的扰动项，表示结构方程未解释的部分；不同结构误差之间可能相关。",
    "由联立方程模型所描述的经济系统自身决定的变量；既可以作被解释变量，也可以作解释变量。",
    "由联立方程模型所描述的经济系统之外决定的变量；影响系统中的内生变量，但不受系统内变量影响。",
    "滞后的内生变量或滞后的外生变量；在联立方程识别中通常作为外生变量处理。",
    "出现在某个结构方程右侧、但可能与该方程结构误差相关的变量；若不处理，OLS 可能不一致。",
    "没有进入该结构方程、但用于识别右侧内生变量的外生变量；也就是工具变量的来源。",
    "用外生变量表示内生变量的方程；在 2SLS 第一阶段中体现为内生解释变量对工具变量和控制变量的回归。",
    "为了展示第六讲联立方程思想而构造的两方程系统；若经济解释或排除限制不足，则不作为主识别结果。"
  ),
  this_project_mapping = c(
    "主结构方程是收益方程：CAR_3day = alpha0 + alpha1 RateVulnerability_z + controls + u1。",
    "核心结构系数是 alpha1，即 RateVulnerability_z 对 CAR_3day 的影响。",
    "收益方程中的 u1 是无法观测的行业冲击、风险偏好变化、遗漏行业特征等。",
    "探索性 SEM 中，CAR_3day 与 RateVulnerability_z 被放入系统内讨论；但主模型中 CAR_3day 是结果变量，RateVulnerability_z 是潜在内生解释变量。",
    "Z_precovid_2Y_z、Z_early_nonoverlap_2Y_z、Z_verylong_2Y_z、Z_termspread_z 以及控制变量可被视为系统外给定变量；但其外生性需要经济论证。",
    "Z_precovid_2Y_z、Z_early_nonoverlap_2Y_z、Z_verylong_2Y_z 来自事件前历史窗口，因此更适合作为前定变量。",
    "主收益方程中的 RateVulnerability_z 是右侧内生变量，因为它可能与遗漏的行业久期、融资约束、避险属性等进入误差项的因素相关。",
    "主收益方程中，被排除外生变量是历史窗口利率暴露度，如 Z_precovid_2Y_z 和 Z_early_nonoverlap_2Y_z。",
    "第一阶段方程：RateVulnerability_z = pi0 + pi1 Z + controls + v，是主收益方程的约简型/第一阶段近似。",
    "CAR_3day <-> RateVulnerability_z 的完整双向系统只作为探索性模块；由于 CAR_3day 发生在历史 RV 之后，不能作为主识别。"
  ),
  use_in_main_text = c(
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE
  ),
  row.names = NULL
)

write.csv(
  concept_definitions,
  file.path(out_dir, "step6_concept_definitions.csv"),
  row.names = FALSE
)

equation_variable_map <- data.frame(
  equation = c(
    "main_return_equation",
    "main_return_equation",
    "main_return_equation",
    "main_return_equation",
    "exploratory_vulnerability_equation",
    "exploratory_vulnerability_equation",
    "exploratory_vulnerability_equation"
  ),
  variable_group = c(
    "dependent endogenous variable",
    "endogenous RHS variable",
    "included exogenous controls",
    "excluded exogenous instruments",
    "dependent endogenous variable",
    "endogenous RHS variable",
    "excluded exogenous instruments, questionable"
  ),
  variables = c(
    "CAR_3day",
    "RateVulnerability_z",
    paste(controls_return, collapse = "; "),
    paste(historical_ivs, collapse = "; "),
    "RateVulnerability_z",
    "CAR_3day",
    paste(reverse_excluded, collapse = "; ")
  ),
  econometric_meaning = c(
    "FOMC 事件窗口行业收益，是主结构方程的被解释变量。",
    "核心解释变量，可能与遗漏行业特征相关，因此作为潜在内生解释变量处理。",
    "进入收益方程的控制变量，用于控制行业资产定价特征、历史风险和事件前动量。",
    "不直接进入主收益方程、但解释 RateVulnerability_z 的事件前历史暴露变量。",
    "探索性反向方程的被解释变量；由于它来自历史窗口，该方程不作为主结论。",
    "探索性反向方程的右侧内生变量；时间顺序上不能解释为 CAR 决定历史 RV。",
    "这些变量机械上可以帮助识别反向方程，但经济上很可能直接影响 RV，因此排除限制较弱。"
  ),
  main_or_exploratory = c(
    "main", "main", "main", "main",
    "exploratory", "exploratory", "exploratory"
  ),
  row.names = NULL
)

write.csv(
  equation_variable_map,
  file.path(out_dir, "step6_equation_variable_map.csv"),
  row.names = FALSE
)

concept_md <- c(
  "# Step 6：第六讲概念定义与本文变量对应关系",
  "",
  "## 1. 为什么需要单独定义内生变量和外生变量？",
  "",
  "第六讲的联立方程模型不是普通多元回归。它关注的是系统内部变量相互影响、共同决定的问题。因此，在估计之前必须说明哪些变量是系统内生决定的，哪些变量是系统外部给定的，哪些变量虽然出现在右侧但可能与结构误差相关。",
  "",
  "## 2. 本文主收益结构方程",
  "",
  "`CAR_3day = alpha0 + alpha1 RateVulnerability_z + controls + u1`",
  "",
  "- `CAR_3day`：主结构方程的被解释变量，也是探索性 SEM 中的内生变量。",
  "- `RateVulnerability_z`：主收益方程的右侧内生变量，因为它可能与遗漏行业特征相关。",
  paste0("- included exogenous controls：", paste(controls_return, collapse = ", "), "。"),
  paste0("- excluded exogenous instruments：", paste(historical_ivs, collapse = ", "), "。"),
  "",
  "## 3. 探索性反向方程",
  "",
  "`RateVulnerability_z = gamma0 + gamma1 CAR_3day + historical_IVs + u2`",
  "",
  "这个方程只用于展示第六讲中双向或联立系统的概念，不作为主识别。原因是 `RateVulnerability_z` 是事件前历史窗口估计得到的行业利率脆弱度，而 `CAR_3day` 是 FOMC 事件窗口的收益反应；从时间顺序看，当期 CAR 不能反过来决定历史 RV。",
  "",
  "## 4. 结论",
  "",
  "本文主结论仍应基于收益结构方程的 OLS-HC3、2SLS 辅助检验和稳健性结果。完整 `CAR_3day <-> RateVulnerability_z` 联立系统只作为第六讲教学性和探索性模块。"
)

writeLines(
  concept_md,
  file.path(out_dir, "step6_concept_definitions.md")
)