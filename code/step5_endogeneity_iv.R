options(stringsAsFactors = FALSE)

this_file <- tryCatch(normalizePath(sys.frames()[[1]]$ofile), error = function(e) NA_character_)
if (is.na(this_file)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) this_file <- normalizePath(sub("^--file=", "", file_arg[1]))
}

project_root <- normalizePath(file.path(dirname(this_file), ".."))
raw_dir <- file.path(project_root, "data", "raw")
out_dir <- file.path(project_root, "output", "step5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

event_date <- as.Date("2022-11-02")
ff49_zip <- file.path(raw_dir, "49_Industry_Portfolios_daily_CSV.zip")
ff5_zip <- file.path(raw_dir, "F-F_Research_Data_5_Factors_2x3_daily_CSV.zip")
frb_h15_csv <- file.path(raw_dir, "FRB_H15.csv")
step3_file <- file.path(project_root, "output", "step3", "step3_cross_section_data.csv")

if (!file.exists(ff49_zip)) stop("Missing ", ff49_zip)
if (!file.exists(ff5_zip)) stop("Missing ", ff5_zip)
if (!file.exists(frb_h15_csv)) stop("Missing ", frb_h15_csv)
if (!file.exists(step3_file)) stop("Missing Step 3 data. Run code/step3_multivariate_cross_section.R first.")

parse_zip_csv_lines <- function(zip_path) {
  temp_dir <- tempfile("zipcsv_")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
  files <- unzip(zip_path, exdir = temp_dir)
  csv_file <- files[grepl("\\.csv$", files, ignore.case = TRUE)][1]
  readLines(csv_file, warn = FALSE)
}

parse_ff49_daily <- function(zip_path) {
  lines <- parse_zip_csv_lines(zip_path)
  header_line <- grep("^\\s*,", lines)[1]
  data_start <- header_line + 1
  data_end <- data_start
  while (data_end <= length(lines) && grepl("^\\s*[0-9]{8}\\s*,", lines[data_end])) data_end <- data_end + 1
  table_text <- paste(lines[header_line:(data_end - 1)], collapse = "\n")
  dat <- read.csv(text = table_text, check.names = FALSE)
  names(dat)[1] <- "date"
  dat$date <- as.Date(as.character(dat$date), format = "%Y%m%d")
  for (col in setdiff(names(dat), "date")) {
    dat[[col]] <- as.numeric(dat[[col]])
    dat[[col]][dat[[col]] <= -99] <- NA_real_
  }
  dat[order(dat$date), ]
}

parse_ff5_daily <- function(zip_path) {
  lines <- parse_zip_csv_lines(zip_path)
  header_line <- grep("^\\s*,\\s*Mkt-RF\\s*,\\s*SMB\\s*,\\s*HML\\s*,\\s*RMW\\s*,\\s*CMA\\s*,\\s*RF\\s*$", lines)[1]
  data_start <- header_line + 1
  data_end <- data_start
  while (data_end <= length(lines) && grepl("^\\s*[0-9]{8}\\s*,", lines[data_end])) data_end <- data_end + 1
  table_text <- paste(lines[header_line:(data_end - 1)], collapse = "\n")
  dat <- read.csv(text = table_text, check.names = FALSE)
  names(dat)[1] <- "date"
  names(dat)[names(dat) == "Mkt-RF"] <- "MktRF"
  dat$date <- as.Date(as.character(dat$date), format = "%Y%m%d")
  for (col in setdiff(names(dat), "date")) dat[[col]] <- as.numeric(dat[[col]])
  dat[order(dat$date), ]
}

parse_frb_h15_dgs2 <- function(csv_path) {
  dat <- read.csv(csv_path, skip = 5, check.names = FALSE, na.strings = c("", "ND", "NA"))
  out <- data.frame(
    date = as.Date(dat[["Time Period"]]),
    DGS2 = suppressWarnings(as.numeric(dat[["RIFLGFCY02_N.B"]])),
    DGS10 = suppressWarnings(as.numeric(dat[["RIFLGFCY10_N.B"]]))
  )
  out <- out[order(out$date), ]
  out2 <- out[!is.na(out$DGS2), c("date", "DGS2")]
  out2$d_DGS2 <- c(NA_real_, diff(out2$DGS2))
  out10 <- out[!is.na(out$DGS10), c("date", "DGS10")]
  out10$d_DGS10 <- c(NA_real_, diff(out10$DGS10))
  merged <- Reduce(function(x, y) merge(x, y, by = "date", all = TRUE), list(out2, out10))
  merged$TermSpread <- merged$DGS10 - merged$DGS2
  spread_nonmissing <- merged[!is.na(merged$TermSpread), c("date", "TermSpread")]
  spread_nonmissing$d_TermSpread <- c(NA_real_, diff(spread_nonmissing$TermSpread))
  merge(merged, spread_nonmissing[, c("date", "d_TermSpread")], by = "date", all = TRUE)
}

zscore <- function(x) as.numeric(scale(x))

weighted_vcov <- function(model, type = "HC3") {
  x <- model.matrix(model)
  weights_model <- weights(model)
  if (is.null(weights_model)) weights_model <- rep(1, nrow(x))
  residual <- residuals(model)
  bread <- solve(t(x) %*% (weights_model * x))
  if (type == "const") return(vcov(model))
  hat_values <- hatvalues(model)
  adjustment <- switch(
    type,
    HC0 = rep(1, length(residual)),
    HC1 = rep(nrow(x) / (nrow(x) - ncol(x)), length(residual)),
    HC2 = 1 / pmax(1 - hat_values, .Machine$double.eps),
    HC3 = 1 / pmax(1 - hat_values, .Machine$double.eps)^2,
    stop("Unsupported vcov type")
  )
  meat <- t(x) %*% ((weights_model^2 * residual^2 * adjustment) * x)
  bread %*% meat %*% bread
}

coef_table <- function(model, model_name, method, term = "RateVulnerability_z", vcov_type = "HC3", depvar = "CAR", core_x = "RateVulnerability_z", controls = "", proxy = "", iv = "") {
  vc <- weighted_vcov(model, vcov_type)
  estimates <- coef(model)
  se <- sqrt(diag(vc))
  t_value <- estimates / se
  p_value <- 2 * pt(abs(t_value), df = df.residual(model), lower.tail = FALSE)
  out <- data.frame(
    Model = model_name,
    Method = method,
    DepVar = depvar,
    Core_X = core_x,
    Controls = controls,
    Proxy = proxy,
    IV = iv,
    term = names(estimates),
    beta = unname(estimates),
    SE_HC3 = unname(se),
    t_value = unname(t_value),
    p_value = unname(p_value),
    N = nobs(model),
    row.names = NULL
  )
  out[out$term == term, ]
}

winsorize <- function(x, probs = c(0.05, 0.95)) {
  qs <- quantile(x, probs = probs, na.rm = TRUE)
  pmin(pmax(x, qs[1]), qs[2])
}

fit_rate_sensitivity <- function(returns, rates, ff5 = NULL, event_index, start_lag, end_lag, rate_col = "d_DGS2", include_factors = TRUE) {
  estimation_index <- (event_index - start_lag):(event_index - end_lag)
  industry_cols <- setdiff(names(returns), "date")
  estimation_returns <- returns[estimation_index, c("date", industry_cols)]
  estimation_data <- merge(estimation_returns, rates[, c("date", rate_col)], by = "date", all.x = TRUE)
  if (include_factors) estimation_data <- merge(estimation_data, ff5, by = "date", all.x = TRUE)
  out <- data.frame(industry = industry_cols, sensitivity = NA_real_, n = NA_integer_)
  for (industry in industry_cols) {
    if (include_factors) {
      reg_data <- estimation_data[, c(industry, "RF", "MktRF", "SMB", "HML", "RMW", "CMA", rate_col)]
      names(reg_data) <- c("ret", "RF", "MktRF", "SMB", "HML", "RMW", "CMA", "d_rate")
      reg_data$excess_ret <- reg_data$ret - reg_data$RF
      reg_data <- reg_data[complete.cases(reg_data), ]
      fit <- lm(excess_ret ~ MktRF + SMB + HML + RMW + CMA + d_rate, data = reg_data)
    } else {
      reg_data <- estimation_data[, c(industry, rate_col)]
      names(reg_data) <- c("ret", "d_rate")
      reg_data <- reg_data[complete.cases(reg_data), ]
      fit <- lm(ret ~ d_rate, data = reg_data)
    }
    out$sensitivity[out$industry == industry] <- -unname(coef(fit)["d_rate"])
    out$n[out$industry == industry] <- nrow(reg_data)
  }
  out
}

fit_rate_sensitivity_dates <- function(returns, rates, ff5 = NULL, start_date, end_date, rate_col = "d_DGS2", include_factors = TRUE) {
  industry_cols <- setdiff(names(returns), "date")
  estimation_returns <- returns[returns$date >= as.Date(start_date) & returns$date <= as.Date(end_date), c("date", industry_cols)]
  estimation_data <- merge(estimation_returns, rates[, c("date", rate_col)], by = "date", all.x = TRUE)
  if (include_factors) estimation_data <- merge(estimation_data, ff5, by = "date", all.x = TRUE)
  out <- data.frame(industry = industry_cols, sensitivity = NA_real_, n = NA_integer_)
  for (industry in industry_cols) {
    if (include_factors) {
      reg_data <- estimation_data[, c(industry, "RF", "MktRF", "SMB", "HML", "RMW", "CMA", rate_col)]
      names(reg_data) <- c("ret", "RF", "MktRF", "SMB", "HML", "RMW", "CMA", "d_rate")
      reg_data$excess_ret <- reg_data$ret - reg_data$RF
      reg_data <- reg_data[complete.cases(reg_data), ]
      fit <- lm(excess_ret ~ MktRF + SMB + HML + RMW + CMA + d_rate, data = reg_data)
    } else {
      reg_data <- estimation_data[, c(industry, rate_col)]
      names(reg_data) <- c("ret", "d_rate")
      reg_data <- reg_data[complete.cases(reg_data), ]
      fit <- lm(ret ~ d_rate, data = reg_data)
    }
    out$sensitivity[out$industry == industry] <- -unname(coef(fit)["d_rate"])
    out$n[out$industry == industry] <- nrow(reg_data)
  }
  out
}

manual_2sls <- function(data, y, x_endog, controls, instruments) {
  y_vec <- as.matrix(data[, y, drop = FALSE])
  x_matrix <- model.matrix(as.formula(paste("~", paste(c(x_endog, controls), collapse = " + "))), data = data)
  z_matrix <- model.matrix(as.formula(paste("~", paste(c(instruments, controls), collapse = " + "))), data = data)
  projection_z <- z_matrix %*% solve(crossprod(z_matrix)) %*% t(z_matrix)
  x_hat <- projection_z %*% x_matrix
  beta <- solve(t(x_hat) %*% x_matrix) %*% t(x_hat) %*% y_vec
  rownames(beta) <- colnames(x_matrix)
  fitted_values <- x_matrix %*% beta
  residual <- as.vector(y_vec - fitted_values)
  n <- nrow(x_matrix)
  k <- ncol(x_matrix)
  bread <- solve(t(x_hat) %*% x_matrix)
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

iv_coef_row <- function(iv_fit, model_name, iv_label) {
  idx <- match(iv_fit$x_endog, iv_fit$names)
  data.frame(
    Model = model_name,
    Method = "2SLS_manual_HC3",
    DepVar = iv_fit$y,
    Core_X = iv_fit$x_endog,
    Controls = paste(iv_fit$controls, collapse = "; "),
    Proxy = "",
    IV = iv_label,
    term = iv_fit$x_endog,
    beta = iv_fit$coefficients[idx],
    SE_HC3 = iv_fit$se_hc3[idx],
    t_value = iv_fit$t_value[idx],
    p_value = iv_fit$p_value[idx],
    N = iv_fit$n,
    row.names = NULL
  )
}

first_stage_table <- function(data, x_endog, controls, instruments, label) {
  full <- lm(as.formula(paste(x_endog, "~", paste(c(instruments, controls), collapse = " + "))), data = data)
  restricted <- lm(as.formula(paste(x_endog, "~", paste(controls, collapse = " + "))), data = data)
  a <- anova(restricted, full)
  sst_restricted <- sum(resid(restricted)^2)
  sst_full <- sum(resid(full)^2)
  partial_r2 <- (sst_restricted - sst_full) / sst_restricted
  coef_mat <- summary(full)$coefficients
  rows <- data.frame(
    IV_set = label,
    IV = instruments,
    coefficient = coef(full)[instruments],
    std_error = coef_mat[instruments, "Std. Error"],
    p_value = coef_mat[instruments, "Pr(>|t|)"],
    first_stage_F = a$F[2],
    first_stage_p = a$`Pr(>F)`[2],
    partial_R2 = partial_r2,
    N = nobs(full),
    row.names = NULL
  )
  rows
}

dwh_test <- function(data, y, x_endog, controls, instruments, label) {
  first <- lm(as.formula(paste(x_endog, "~", paste(c(instruments, controls), collapse = " + "))), data = data)
  data$first_stage_resid <- resid(first)
  augmented <- lm(as.formula(paste(y, "~", paste(c(x_endog, controls, "first_stage_resid"), collapse = " + "))), data = data)
  coef_mat <- summary(augmented)$coefficients
  data.frame(
    test = paste0("DWH residual inclusion: ", label),
    coefficient_resid = coef(augmented)["first_stage_resid"],
    p_value = coef_mat["first_stage_resid", "Pr(>|t|)"],
    N = nobs(augmented),
    row.names = NULL
  )
}

sargan_test <- function(iv_fit, label) {
  overid_df <- ncol(iv_fit$z_matrix) - ncol(iv_fit$x_matrix)
  if (overid_df <= 0) {
    return(data.frame(test = paste0("Sargan overidentification: ", label), statistic = NA_real_, df = overid_df, p_value = NA_real_, note = "Exactly identified; overidentification test not available."))
  }
  aux <- lm(iv_fit$residuals ~ iv_fit$z_matrix[, -1, drop = FALSE])
  statistic <- iv_fit$n * summary(aux)$r.squared
  data.frame(
    test = paste0("Sargan overidentification: ", label),
    statistic = statistic,
    df = overid_df,
    p_value = pchisq(statistic, df = overid_df, lower.tail = FALSE),
    note = "Tests joint validity of overidentifying restrictions under homoskedasticity.",
    row.names = NULL
  )
}

hansen_j_test <- function(iv_fit, label) {
  overid_df <- ncol(iv_fit$z_matrix) - ncol(iv_fit$x_matrix)
  if (overid_df <= 0) {
    return(data.frame(
      test = paste0("Hansen-style J overidentification: ", label),
      statistic = NA_real_,
      df = overid_df,
      p_value = NA_real_,
      note = "Exactly identified; robust overidentification test not available.",
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
  data.frame(
    test = paste0("Hansen-style J overidentification: ", label),
    statistic = statistic,
    df = overid_df,
    p_value = pchisq(statistic, df = overid_df, lower.tail = FALSE),
    note = "Robust moment-based J using 2SLS residuals; interpret cautiously in N = 49.",
    row.names = NULL
  )
}

ff49 <- parse_ff49_daily(ff49_zip)
ff5 <- parse_ff5_daily(ff5_zip)
rates <- parse_frb_h15_dgs2(frb_h15_csv)
step3 <- read.csv(step3_file)

event_index <- match(event_date, ff49$date)
industry_cols <- setdiff(names(ff49), "date")

event_window <- function(start_lag, end_lag) {
  idx <- (event_index + start_lag):(event_index + end_lag)
  out <- data.frame(industry = industry_cols, value = NA_real_)
  for (industry in industry_cols) out$value[out$industry == industry] <- sum(ff49[[industry]][idx], na.rm = TRUE)
  out
}

car_1 <- event_window(0, 0)
car_3 <- event_window(-1, 1)
car_5 <- event_window(-2, 2)
names(car_1)[2] <- "CAR_1"
names(car_3)[2] <- "CAR_3_rebuilt"
names(car_5)[2] <- "CAR_5"

rv_long <- fit_rate_sensitivity(ff49, rates, ff5, event_index, 520, 30, "d_DGS2", include_factors = TRUE)
names(rv_long)[names(rv_long) == "sensitivity"] <- "RateVulnerability_long"
rv_10y <- fit_rate_sensitivity(ff49, rates, ff5, event_index, 260, 30, "d_DGS10", include_factors = TRUE)
names(rv_10y)[names(rv_10y) == "sensitivity"] <- "RateVulnerability_10Y"
rv_simple_long <- fit_rate_sensitivity(ff49, rates, ff5, event_index, 760, 521, "d_DGS2", include_factors = FALSE)
names(rv_simple_long)[names(rv_simple_long) == "sensitivity"] <- "Z_presample_simple_2Y"
rv_presample_2y <- fit_rate_sensitivity(ff49, rates, ff5, event_index, 760, 521, "d_DGS2", include_factors = TRUE)
names(rv_presample_2y)[names(rv_presample_2y) == "sensitivity"] <- "Z_presample_2Y"
rv_presample_10y <- fit_rate_sensitivity(ff49, rates, ff5, event_index, 760, 521, "d_DGS10", include_factors = TRUE)
names(rv_presample_10y)[names(rv_presample_10y) == "sensitivity"] <- "Z_presample_10Y"
rv_recent_nonoverlap <- fit_rate_sensitivity(ff49, rates, ff5, event_index, 520, 261, "d_DGS2", include_factors = TRUE)
names(rv_recent_nonoverlap)[names(rv_recent_nonoverlap) == "sensitivity"] <- "Z_recent_nonoverlap_2Y"
rv_verylong <- fit_rate_sensitivity(ff49, rates, ff5, event_index, 1500, 521, "d_DGS2", include_factors = TRUE)
names(rv_verylong)[names(rv_verylong) == "sensitivity"] <- "Z_verylong_2Y"
rv_precovid <- fit_rate_sensitivity_dates(ff49, rates, ff5, "2017-01-03", "2019-12-31", "d_DGS2", include_factors = TRUE)
names(rv_precovid)[names(rv_precovid) == "sensitivity"] <- "Z_precovid_2Y"
rv_termspread <- fit_rate_sensitivity(ff49, rates, ff5, event_index, 520, 261, "d_TermSpread", include_factors = TRUE)
names(rv_termspread)[names(rv_termspread) == "sensitivity"] <- "Z_termspread"
rv_early_split <- fit_rate_sensitivity(ff49, rates, ff5, event_index, 520, 391, "d_DGS2", include_factors = TRUE)
names(rv_early_split)[names(rv_early_split) == "sensitivity"] <- "RV_early_2Y"
rv_late_split <- fit_rate_sensitivity(ff49, rates, ff5, event_index, 390, 261, "d_DGS2", include_factors = TRUE)
names(rv_late_split)[names(rv_late_split) == "sensitivity"] <- "RV_late_2Y"

df <- Reduce(function(x, y) merge(x, y, by = "industry", all.x = TRUE), list(
  step3, car_1, car_3, car_5,
  rv_long[, c("industry", "RateVulnerability_long")],
  rv_10y[, c("industry", "RateVulnerability_10Y")],
  rv_simple_long[, c("industry", "Z_presample_simple_2Y")],
  rv_presample_2y[, c("industry", "Z_presample_2Y")],
  rv_presample_10y[, c("industry", "Z_presample_10Y")],
  rv_recent_nonoverlap[, c("industry", "Z_recent_nonoverlap_2Y")],
  rv_verylong[, c("industry", "Z_verylong_2Y")],
  rv_precovid[, c("industry", "Z_precovid_2Y")],
  rv_termspread[, c("industry", "Z_termspread")],
  rv_early_split[, c("industry", "RV_early_2Y")],
  rv_late_split[, c("industry", "RV_late_2Y")]
))

df$CAR <- df$CAR
df$RateVulnerability_long_z <- zscore(df$RateVulnerability_long)
df$RateVulnerability_10Y_z <- zscore(df$RateVulnerability_10Y)
df$RateVulnerability_winsor_z <- zscore(winsorize(df$RateVulnerability))
df$Z_presample_simple_2Y_z <- zscore(df$Z_presample_simple_2Y)
df$Z_presample_2Y_z <- zscore(df$Z_presample_2Y)
df$Z_presample_10Y_z <- zscore(df$Z_presample_10Y)
df$Z_recent_nonoverlap_2Y_z <- zscore(df$Z_recent_nonoverlap_2Y)
df$Z_verylong_2Y_z <- zscore(df$Z_verylong_2Y)
df$Z_precovid_2Y_z <- zscore(df$Z_precovid_2Y)
df$Z_termspread_z <- zscore(df$Z_termspread)
df$RV_early_2Y_z <- zscore(df$RV_early_2Y)
df$RV_late_2Y_z <- zscore(df$RV_late_2Y)
df$GoldDummy <- as.integer(df$industry == "Gold")
df$FinanceDummy <- as.integer(df$Finance == 1)
df$DefensiveDummy <- as.integer(df$Defensive == 1)
df$HighGrowthProxy <- as.integer(df$HMLBeta_z < median(df$HMLBeta_z, na.rm = TRUE))
df$HighVolProxy <- as.integer(df$HistVol_z > median(df$HistVol_z, na.rm = TRUE))

controls <- c("MktBeta_z", "SMBBeta_z", "HMLBeta_z", "RMWBeta_z", "CMABeta_z", "HistVol_z", "PreMomentum_z")
controls_text <- paste(controls, collapse = "; ")

fit_lm <- function(y, x, extra = character(), data = df) {
  lm(as.formula(paste(y, "~", paste(c(x, controls, extra), collapse = " + "))), data = data)
}

proxy_models <- do.call(rbind, list(
  coef_table(lm(CAR ~ RateVulnerability_z, data = df), "M1_simple", "OLS_HC3", controls = "None"),
  coef_table(fit_lm("CAR", "RateVulnerability_z"), "M2_controls", "OLS_HC3", controls = controls_text),
  coef_table(fit_lm("CAR", "RateVulnerability_z", "GoldDummy"), "M3_proxy_Gold", "OLS_HC3", controls = controls_text, proxy = "GoldDummy"),
  coef_table(fit_lm("CAR", "RateVulnerability_z", "HighGrowthProxy"), "M3_proxy_HighGrowth", "OLS_HC3", controls = controls_text, proxy = "HighGrowthProxy"),
  coef_table(fit_lm("CAR", "RateVulnerability_z", "HighVolProxy"), "M3_proxy_HighVol", "OLS_HC3", controls = controls_text, proxy = "HighVolProxy"),
  coef_table(fit_lm("CAR", "RateVulnerability_z", "FinanceDummy"), "M3_proxy_Finance", "OLS_HC3", controls = controls_text, proxy = "FinanceDummy"),
  coef_table(fit_lm("CAR", "RateVulnerability_z", "DefensiveDummy"), "M3_proxy_Defensive", "OLS_HC3", controls = controls_text, proxy = "DefensiveDummy")
))

measurement_y <- do.call(rbind, list(
  coef_table(fit_lm("CAR_1", "RateVulnerability_z"), "Y_AR0", "OLS_HC3", depvar = "CAR_1", controls = controls_text),
  coef_table(fit_lm("CAR", "RateVulnerability_z"), "Y_CAR3", "OLS_HC3", depvar = "CAR_3", controls = controls_text),
  coef_table(fit_lm("CAR_5", "RateVulnerability_z"), "Y_CAR5", "OLS_HC3", depvar = "CAR_5", controls = controls_text),
  coef_table(fit_lm("CAR", "RateVulnerability_z", data = df[df$industry != "Gold", ]), "Y_CAR3_exGold", "OLS_HC3", depvar = "CAR_3_exGold", controls = controls_text)
))

measurement_x <- do.call(rbind, list(
  coef_table(fit_lm("CAR", "RateVulnerability_z"), "X_main_260d_2Y", "OLS_HC3", controls = controls_text),
  coef_table(fit_lm("CAR", "RateVulnerability_long_z"), "X_long_520d_2Y", "OLS_HC3", term = "RateVulnerability_long_z", core_x = "RateVulnerability_long_z", controls = controls_text),
  coef_table(fit_lm("CAR", "RateVulnerability_10Y_z"), "X_260d_10Y", "OLS_HC3", term = "RateVulnerability_10Y_z", core_x = "RateVulnerability_10Y_z", controls = controls_text),
  coef_table(fit_lm("CAR", "RateVulnerability_winsor_z"), "X_winsorized_2Y", "OLS_HC3", term = "RateVulnerability_winsor_z", core_x = "RateVulnerability_winsor_z", controls = controls_text)
))

mechanism_models <- do.call(rbind, list(
  coef_table(lm(CAR ~ RateVulnerability_z * HighGrowthProxy + MktBeta_z + SMBBeta_z + HMLBeta_z + RMWBeta_z + CMABeta_z + HistVol_z + PreMomentum_z, data = df), "Mechanism_HighGrowth_interaction_main", "OLS_HC3", controls = controls_text, proxy = "HighGrowthProxy"),
  coef_table(lm(CAR ~ RateVulnerability_z * HighVolProxy + MktBeta_z + SMBBeta_z + HMLBeta_z + RMWBeta_z + CMABeta_z + HistVol_z + PreMomentum_z, data = df), "Mechanism_HighVol_interaction_main", "OLS_HC3", controls = controls_text, proxy = "HighVolProxy")
))
mechanism_interactions <- do.call(rbind, list(
  coef_table(lm(CAR ~ RateVulnerability_z * HighGrowthProxy + MktBeta_z + SMBBeta_z + HMLBeta_z + RMWBeta_z + CMABeta_z + HistVol_z + PreMomentum_z, data = df), "Mechanism_HighGrowth_interaction", "OLS_HC3", term = "RateVulnerability_z:HighGrowthProxy", controls = controls_text, proxy = "HighGrowthProxy"),
  coef_table(lm(CAR ~ RateVulnerability_z * HighVolProxy + MktBeta_z + SMBBeta_z + HMLBeta_z + RMWBeta_z + CMABeta_z + HistVol_z + PreMomentum_z, data = df), "Mechanism_HighVol_interaction", "OLS_HC3", term = "RateVulnerability_z:HighVolProxy", controls = controls_text, proxy = "HighVolProxy")
))

iv1 <- manual_2sls(df, "CAR", "RateVulnerability_z", controls, "Z_presample_2Y_z")
iv2 <- manual_2sls(df, "CAR", "RateVulnerability_z", controls, c("Z_presample_2Y_z", "Z_presample_10Y_z"))
iv3 <- manual_2sls(df, "CAR", "RateVulnerability_z", controls, "Z_presample_simple_2Y_z")

iv_results <- do.call(rbind, list(
  iv_coef_row(iv1, "IV_one_presample_2Y", "Z_presample_2Y_z"),
  iv_coef_row(iv2, "IV_two_presample_2Y_10Y", "Z_presample_2Y_z; Z_presample_10Y_z"),
  iv_coef_row(iv3, "IV_simple_presample_2Y", "Z_presample_simple_2Y_z")
))

first_stage <- do.call(rbind, list(
  first_stage_table(df, "RateVulnerability_z", controls, "Z_presample_2Y_z", "IV_one_presample_2Y"),
  first_stage_table(df, "RateVulnerability_z", controls, c("Z_presample_2Y_z", "Z_presample_10Y_z"), "IV_two_presample_2Y_10Y"),
  first_stage_table(df, "RateVulnerability_z", controls, "Z_presample_simple_2Y_z", "IV_simple_presample_2Y")
))

endogeneity_diagnostics <- do.call(rbind, list(
  dwh_test(df, "CAR", "RateVulnerability_z", controls, "Z_presample_2Y_z", "IV_one_presample_2Y"),
  dwh_test(df, "CAR", "RateVulnerability_z", controls, c("Z_presample_2Y_z", "Z_presample_10Y_z"), "IV_two_presample_2Y_10Y"),
  dwh_test(df, "CAR", "RateVulnerability_z", controls, "Z_presample_simple_2Y_z", "IV_simple_presample_2Y")
))

overid_tests <- do.call(rbind, list(
  sargan_test(iv1, "IV_one_presample_2Y"),
  sargan_test(iv2, "IV_two_presample_2Y_10Y"),
  sargan_test(iv3, "IV_simple_presample_2Y")
))

iv_candidate_sets <- list(
  Z_recent_nonoverlap_2Y_z = "event -520 to -261 trading days, factor-adjusted 2Y sensitivity",
  Z_early_nonoverlap_2Y_z = "event -760 to -521 trading days, factor-adjusted 2Y sensitivity",
  Z_verylong_2Y_z = "event -1500 to -521 trading days, factor-adjusted 2Y sensitivity",
  Z_precovid_2Y_z = "calendar 2017-2019, factor-adjusted 2Y sensitivity",
  Z_termspread_z = "event -520 to -261 trading days, factor-adjusted 10Y-2Y term-spread sensitivity",
  RV_early_2Y_z = "split-half early window, event -520 to -391 trading days",
  RV_late_2Y_z = "split-half late window, event -390 to -261 trading days"
)

df$Z_early_nonoverlap_2Y_z <- df$Z_presample_2Y_z

iv_search_candidates <- data.frame(
  IV = names(iv_candidate_sets),
  construction = unname(unlist(iv_candidate_sets)),
  target_endogenous_variable = "RateVulnerability_z",
  intended_use = ifelse(
    names(iv_candidate_sets) %in% c("RV_early_2Y_z", "RV_late_2Y_z"),
    "split-sample measurement-error IV",
    "pre-determined non-overlapping historical rate-exposure IV"
  ),
  row.names = NULL
)

search_first_stage_list <- lapply(names(iv_candidate_sets), function(iv_name) {
  first_stage_table(df, "RateVulnerability_z", controls, iv_name, iv_name)
})
step5b_first_stage_rank <- do.call(rbind, search_first_stage_list)
step5b_first_stage_rank$strength_bucket <- cut(
  step5b_first_stage_rank$first_stage_F,
  breaks = c(-Inf, 5, 10, Inf),
  labels = c("weak_F_below_5", "moderate_F_5_to_10", "strong_F_above_10"),
  right = FALSE
)
step5b_first_stage_rank <- step5b_first_stage_rank[order(-step5b_first_stage_rank$first_stage_F), ]

search_2sls_list <- lapply(names(iv_candidate_sets), function(iv_name) {
  fit <- manual_2sls(df, "CAR", "RateVulnerability_z", controls, iv_name)
  row <- iv_coef_row(fit, paste0("IV_search_", iv_name), iv_name)
  row$first_stage_F <- step5b_first_stage_rank$first_stage_F[match(iv_name, step5b_first_stage_rank$IV)]
  row$partial_R2 <- step5b_first_stage_rank$partial_R2[match(iv_name, step5b_first_stage_rank$IV)]
  row
})
step5b_2sls_rank <- do.call(rbind, search_2sls_list)
step5b_2sls_rank <- step5b_2sls_rank[order(-step5b_2sls_rank$first_stage_F), ]

split_iv_models <- list(
  early_instruments_late = list(y = "CAR", x = "RV_late_2Y_z", z = "RV_early_2Y_z"),
  late_instruments_early = list(y = "CAR", x = "RV_early_2Y_z", z = "RV_late_2Y_z")
)
split_first_stage <- do.call(rbind, lapply(names(split_iv_models), function(name) {
  spec <- split_iv_models[[name]]
  first_stage_table(df, spec$x, controls, spec$z, name)
}))
split_2sls <- do.call(rbind, lapply(names(split_iv_models), function(name) {
  spec <- split_iv_models[[name]]
  fit <- manual_2sls(df, spec$y, spec$x, controls, spec$z)
  row <- iv_coef_row(fit, name, spec$z)
  row$target_endogenous_variable <- spec$x
  row
}))

step5b_dwh_strong_iv <- do.call(rbind, list(
  dwh_test(df, "CAR", "RateVulnerability_z", controls, "Z_precovid_2Y_z", "Step5B strongest IV: Z_precovid_2Y_z")
))

overid_iv_sets <- list(
  precovid_plus_verylong = c("Z_precovid_2Y_z", "Z_verylong_2Y_z"),
  precovid_plus_early_nonoverlap = c("Z_precovid_2Y_z", "Z_early_nonoverlap_2Y_z"),
  precovid_plus_termspread = c("Z_precovid_2Y_z", "Z_termspread_z")
)

step5b_overid_results <- do.call(rbind, lapply(names(overid_iv_sets), function(label) {
  instruments <- overid_iv_sets[[label]]
  fit <- manual_2sls(df, "CAR", "RateVulnerability_z", controls, instruments)
  rbind(
    transform(sargan_test(fit, label), test_type = "Sargan"),
    transform(hansen_j_test(fit, label), test_type = "Hansen_style_J")
  )
}))

step5b_overid_2sls <- do.call(rbind, lapply(names(overid_iv_sets), function(label) {
  instruments <- overid_iv_sets[[label]]
  fit <- manual_2sls(df, "CAR", "RateVulnerability_z", controls, instruments)
  iv_coef_row(fit, paste0("overid_", label), paste(instruments, collapse = "; "))
}))

main_comparison <- do.call(rbind, list(
  proxy_models,
  measurement_y,
  measurement_x,
  iv_results,
  mechanism_interactions
))

write.csv(df, file.path(out_dir, "step5_data.csv"), row.names = FALSE)
write.csv(proxy_models, file.path(out_dir, "step5_proxy_models.csv"), row.names = FALSE)
write.csv(measurement_y, file.path(out_dir, "step5_measurement_error_y.csv"), row.names = FALSE)
write.csv(measurement_x, file.path(out_dir, "step5_measurement_error_x.csv"), row.names = FALSE)
write.csv(first_stage, file.path(out_dir, "step5_first_stage_iv.csv"), row.names = FALSE)
write.csv(iv_results, file.path(out_dir, "step5_2sls_results.csv"), row.names = FALSE)
write.csv(endogeneity_diagnostics, file.path(out_dir, "step5_endogeneity_diagnostics.csv"), row.names = FALSE)
write.csv(overid_tests, file.path(out_dir, "step5_overidentification_tests.csv"), row.names = FALSE)
write.csv(mechanism_models, file.path(out_dir, "step5_mechanism_main_effects.csv"), row.names = FALSE)
write.csv(mechanism_interactions, file.path(out_dir, "step5_mechanism_interactions.csv"), row.names = FALSE)
write.csv(main_comparison, file.path(out_dir, "step5_main_comparison_table.csv"), row.names = FALSE)
write.csv(iv_search_candidates, file.path(out_dir, "step5b_iv_search_candidates.csv"), row.names = FALSE)
write.csv(step5b_first_stage_rank, file.path(out_dir, "step5b_first_stage_rank.csv"), row.names = FALSE)
write.csv(step5b_2sls_rank, file.path(out_dir, "step5b_2sls_rank.csv"), row.names = FALSE)
write.csv(split_first_stage, file.path(out_dir, "step5b_split_sample_first_stage.csv"), row.names = FALSE)
write.csv(split_2sls, file.path(out_dir, "step5b_split_sample_2sls.csv"), row.names = FALSE)
write.csv(step5b_dwh_strong_iv, file.path(out_dir, "step5b_dwh_strong_iv.csv"), row.names = FALSE)
write.csv(step5b_overid_results, file.path(out_dir, "step5b_overidentification_tests.csv"), row.names = FALSE)
write.csv(step5b_overid_2sls, file.path(out_dir, "step5b_overidentified_2sls.csv"), row.names = FALSE)

plot_data <- main_comparison
plot_data$label <- paste(plot_data$Model, plot_data$Core_X, sep = ": ")
plot_data$ci_low <- plot_data$beta - qt(0.975, df = pmax(plot_data$N - 10, 1)) * plot_data$SE_HC3
plot_data$ci_high <- plot_data$beta + qt(0.975, df = pmax(plot_data$N - 10, 1)) * plot_data$SE_HC3
plot_data <- plot_data[order(plot_data$beta), ]

png(file.path(out_dir, "step5_coefficient_plot.png"), width = 1500, height = 1100, res = 150)
op <- par(mar = c(5, 15, 4, 2) + 0.1)
ypos <- seq_len(nrow(plot_data))
plot(
  plot_data$beta,
  ypos,
  xlim = range(c(plot_data$ci_low, plot_data$ci_high), na.rm = TRUE),
  yaxt = "n",
  pch = 19,
  col = ifelse(plot_data$p_value < 0.05, "#d62728", "#1f77b4"),
  xlab = "Coefficient estimate with approximate 95% CI",
  ylab = "",
  main = "Step 5: core coefficient across endogeneity checks"
)
segments(plot_data$ci_low, ypos, plot_data$ci_high, ypos, col = "gray50")
abline(v = 0, lty = 2, col = "gray30")
axis(2, at = ypos, labels = plot_data$label, las = 2, cex.axis = 0.65)
grid(nx = NULL, ny = NA, col = "gray85")
par(op)
dev.off()

iv_assessment <- c(
  "Step 5 IV assessment:",
  "Candidate IVs use pre-sample industry rate sensitivity estimated before the main event-estimation window.",
  "Relevance is testable through the first stage. Exogeneity and exclusion are not mechanically testable and require economic judgment.",
  "The pre-sample IVs are temporally predetermined, so reverse causality from the 2022-11-02 CAR is unlikely.",
  "However, they may still capture persistent industry risk characteristics, duration, safe-haven exposure, or macro sensitivity, so the exclusion restriction is contestable.",
  "GoldDummy, HighGrowthProxy, HighVolProxy, FinanceDummy, and DefensiveDummy are proxy/control variables, not IVs, because they plausibly affect event-window CAR directly.",
  "Therefore IV/2SLS should be treated as exploratory robustness, not the paper's primary identification strategy.",
  "If first-stage F is weak or 2SLS estimates are imprecise, do not use IV as main conclusion; rely on OLS with controls, HC3, proxy checks, and measurement-error robustness."
)
writeLines(iv_assessment, file.path(out_dir, "step5_iv_assessment.txt"))

best_iv <- step5b_first_stage_rank[1, ]
step5b_notes <- c(
  "# Step 5B additional IV search",
  "",
  "Additional non-overlapping and split-sample historical rate-exposure IVs were tested.",
  paste0("Best first-stage candidate: ", best_iv$IV, ", F = ", round(best_iv$first_stage_F, 3), ", partial R2 = ", round(best_iv$partial_R2, 3), "."),
  "Decision rule: F > 10 can be discussed as stronger auxiliary IV evidence; 5 <= F < 10 is moderate and should be interpreted cautiously; F < 5 is weak/exploratory.",
  "Even when an IV is temporally predetermined, exclusion is still an economic assumption because historical rate exposure may proxy persistent industry risk characteristics.",
  "Split-sample IVs are most defensible for measurement-error concerns, but they can still be weak in a 49-industry cross-section."
)
writeLines(step5b_notes, file.path(out_dir, "step5b_iv_search_notes.md"))

coverage <- c(
  "# Fifth lecture coverage check",
  "",
  "The PDF cover renders as '第五讲 简单截面数据的经济计量方法（四）'. The text layer is not extractable with local tools, so coverage is checked against the provided Step 5 plan and the standard fifth-lecture topics on endogeneity.",
  "",
  "Covered topics:",
  "1. Endogenous explanatory variables and causal interpretation limits.",
  "2. Unobserved omitted variables and proxy-variable strategy.",
  "3. Measurement error in y and x, with alternative CAR windows and alternative vulnerability measures.",
  "4. Instrument relevance and exclusion restrictions.",
  "5. Manual 2SLS for one-IV and multi-IV cases.",
  "6. First-stage F and partial R-squared for weak-IV diagnosis.",
  "7. DWH residual-inclusion endogeneity diagnostics.",
  "8. Overidentification test only when more instruments than endogenous variables are available.",
  "9. Mechanism-style proxy interactions using high-growth and high-volatility proxies.",
  "10. Explicit caution that IV evidence is exploratory because exclusion is economically contestable.",
  "11. Additional Step 5B IV search over non-overlapping windows, term-spread exposure, and split-sample measurement-error IVs.",
  "12. Proxy-vs-IV distinction: industry dummies and growth/volatility proxies are included as controls/proxies, not instruments."
)
writeLines(coverage, file.path(out_dir, "step5_lecture_coverage_check.md"))

cat("Done.\n")
cat("Output directory:", out_dir, "\n\n")
cat("Proxy models:\n")
print(proxy_models, row.names = FALSE)
cat("\nMeasurement error in y:\n")
print(measurement_y, row.names = FALSE)
cat("\nMeasurement error in x:\n")
print(measurement_x, row.names = FALSE)
cat("\nFirst stage:\n")
print(first_stage, row.names = FALSE)
cat("\n2SLS:\n")
print(iv_results, row.names = FALSE)
cat("\nDWH diagnostics:\n")
print(endogeneity_diagnostics, row.names = FALSE)
cat("\nOveridentification:\n")
print(overid_tests, row.names = FALSE)
cat("\nStep 5B first-stage rank:\n")
print(step5b_first_stage_rank, row.names = FALSE)
cat("\nStep 5B 2SLS rank:\n")
print(step5b_2sls_rank, row.names = FALSE)
cat("\nStep 5B split-sample first stage:\n")
print(split_first_stage, row.names = FALSE)
cat("\nStep 5B split-sample 2SLS:\n")
print(split_2sls, row.names = FALSE)
cat("\nStep 5B DWH for strongest IVs:\n")
print(step5b_dwh_strong_iv, row.names = FALSE)
cat("\nStep 5B overidentification tests:\n")
print(step5b_overid_results, row.names = FALSE)
cat("\nStep 5B overidentified 2SLS:\n")
print(step5b_overid_2sls, row.names = FALSE)
