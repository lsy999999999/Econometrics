options(stringsAsFactors = FALSE)
options(timeout = 600)

this_file <- tryCatch(normalizePath(sys.frames()[[1]]$ofile), error = function(e) NA_character_)
if (is.na(this_file)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) this_file <- normalizePath(sub("^--file=", "", file_arg[1]))
}

project_root <- normalizePath(file.path(dirname(this_file), ".."))
raw_dir <- file.path(project_root, "data", "raw")
out_dir <- file.path(project_root, "output", "step3")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

event_date <- as.Date("2022-11-02")
estimation_start_lag <- 260
estimation_end_lag <- 30
pre_momentum_start_lag <- 60
pre_momentum_end_lag <- 2

ff49_zip <- file.path(raw_dir, "49_Industry_Portfolios_daily_CSV.zip")
ff49_zip_alt <- file.path(raw_dir, "49_Industry_Portfolios_Daily_CSV.zip")
if (!file.exists(ff49_zip) && file.exists(ff49_zip_alt)) ff49_zip <- ff49_zip_alt

frb_h15_csv <- file.path(raw_dir, "FRB_H15.csv")
dgs2_csv <- file.path(raw_dir, "DGS2.csv")

ff5_zip <- file.path(raw_dir, "F-F_Research_Data_5_Factors_2x3_daily_CSV.zip")
ff5_zip_alt <- file.path(raw_dir, "F-F_Research_Data_5_Factors_2x3_Daily_CSV.zip")
if (!file.exists(ff5_zip) && file.exists(ff5_zip_alt)) ff5_zip <- ff5_zip_alt

required_download_message <- paste(
  "Missing Fama-French 5 Factors daily data.",
  "Download the CSV zip from:",
  "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_5_Factors_2x3_daily_CSV.zip",
  "Save it as:",
  ff5_zip,
  sep = "\n"
)

stop_if_missing_inputs <- function() {
  if (!file.exists(ff49_zip)) stop("Missing 49 industry portfolio zip: ", ff49_zip)
  if (!file.exists(frb_h15_csv) && !file.exists(dgs2_csv)) {
    stop("Missing rate data. Expected FRB_H15.csv or DGS2.csv in: ", raw_dir)
  }
  if (!file.exists(ff5_zip)) stop(required_download_message)
}

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
  if (is.na(header_line)) stop("Cannot locate Fama-French 49 industry table header.")

  data_start <- header_line + 1
  data_end <- data_start
  while (data_end <= length(lines) && grepl("^\\s*[0-9]{8}\\s*,", lines[data_end])) {
    data_end <- data_end + 1
  }
  data_end <- data_end - 1
  if (data_end < data_start) stop("Cannot locate Fama-French 49 industry daily observations.")

  table_text <- paste(lines[header_line:data_end], collapse = "\n")
  dat <- read.csv(text = table_text, check.names = FALSE)
  names(dat)[1] <- "date"
  dat$date <- as.Date(as.character(dat$date), format = "%Y%m%d")
  industry_cols <- setdiff(names(dat), "date")
  for (col in industry_cols) {
    dat[[col]] <- as.numeric(dat[[col]])
    dat[[col]][dat[[col]] <= -99] <- NA_real_
  }
  dat[order(dat$date), ]
}

parse_ff5_daily <- function(zip_path) {
  lines <- parse_zip_csv_lines(zip_path)
  header_line <- grep("^\\s*,\\s*Mkt-RF\\s*,\\s*SMB\\s*,\\s*HML\\s*,\\s*RMW\\s*,\\s*CMA\\s*,\\s*RF\\s*$", lines)[1]
  if (is.na(header_line)) stop("Cannot locate Fama-French 5-factor daily table header.")

  data_start <- header_line + 1
  data_end <- data_start
  while (data_end <= length(lines) && grepl("^\\s*[0-9]{8}\\s*,", lines[data_end])) {
    data_end <- data_end + 1
  }
  data_end <- data_end - 1
  if (data_end < data_start) stop("Cannot locate Fama-French 5-factor daily observations.")

  table_text <- paste(lines[header_line:data_end], collapse = "\n")
  dat <- read.csv(text = table_text, check.names = FALSE)
  names(dat)[1] <- "date"
  names(dat)[names(dat) == "Mkt-RF"] <- "MktRF"
  dat$date <- as.Date(as.character(dat$date), format = "%Y%m%d")
  factor_cols <- setdiff(names(dat), "date")
  for (col in factor_cols) dat[[col]] <- as.numeric(dat[[col]])
  dat[order(dat$date), ]
}

parse_fred_dgs2 <- function(csv_path) {
  dat <- read.csv(csv_path, check.names = FALSE)
  names(dat) <- c("date", "DGS2")
  dat$date <- as.Date(dat$date)
  dat$DGS2 <- suppressWarnings(as.numeric(dat$DGS2))
  dat <- dat[order(dat$date), ]
  dat <- dat[!is.na(dat$DGS2), ]
  dat$d_DGS2 <- c(NA_real_, diff(dat$DGS2))
  dat
}

parse_frb_h15_dgs2 <- function(csv_path) {
  dat <- read.csv(
    csv_path,
    skip = 5,
    check.names = FALSE,
    na.strings = c("", "ND", "NA")
  )
  if (!("Time Period" %in% names(dat))) stop("Cannot locate date column in FRB H15 file.")
  if (!("RIFLGFCY02_N.B" %in% names(dat))) stop("Cannot locate 2-year Treasury column RIFLGFCY02_N.B.")
  out <- data.frame(
    date = as.Date(dat[["Time Period"]]),
    DGS2 = suppressWarnings(as.numeric(dat[["RIFLGFCY02_N.B"]]))
  )
  out <- out[order(out$date), ]
  out <- out[!is.na(out$DGS2), ]
  out$d_DGS2 <- c(NA_real_, diff(out$DGS2))
  out
}

read_rate_data <- function() {
  if (file.exists(frb_h15_csv)) {
    message("Using local FRB H15 file: ", frb_h15_csv)
    return(parse_frb_h15_dgs2(frb_h15_csv))
  }
  message("Using local FRED DGS2 file: ", dgs2_csv)
  parse_fred_dgs2(dgs2_csv)
}

zscore <- function(x) as.numeric(scale(x))

hc1_table <- function(model, model_name) {
  x <- model.matrix(model)
  residual <- residuals(model)
  n <- nrow(x)
  k <- ncol(x)
  bread <- solve(crossprod(x))
  meat <- t(x) %*% diag(as.numeric(residual^2), nrow = n) %*% x
  vcov_hc1 <- (n / (n - k)) * bread %*% meat %*% bread
  robust_se <- sqrt(diag(vcov_hc1))
  t_value <- coef(model) / robust_se
  p_value <- 2 * pt(abs(t_value), df = n - k, lower.tail = FALSE)
  data.frame(
    model = model_name,
    variable = names(coef(model)),
    coefficient = unname(coef(model)),
    robust_std_error = unname(robust_se),
    t_value = unname(t_value),
    p_value = unname(p_value),
    row.names = NULL
  )
}

hc1_vcov <- function(model) {
  x <- model.matrix(model)
  residual <- residuals(model)
  n <- nrow(x)
  k <- ncol(x)
  bread <- solve(crossprod(x))
  meat <- t(x) %*% diag(as.numeric(residual^2), nrow = n) %*% x
  (n / (n - k)) * bread %*% meat %*% bread
}

hc1_wald_test <- function(model, terms, test_name) {
  coefficients <- coef(model)
  vcov_hc1 <- hc1_vcov(model)
  missing_terms <- setdiff(terms, names(coefficients))
  if (length(missing_terms) > 0) {
    stop("Terms not found in model: ", paste(missing_terms, collapse = ", "))
  }
  r_matrix <- matrix(0, nrow = length(terms), ncol = length(coefficients))
  colnames(r_matrix) <- names(coefficients)
  rownames(r_matrix) <- terms
  for (i in seq_along(terms)) r_matrix[i, terms[i]] <- 1
  restrictions <- as.vector(r_matrix %*% coefficients)
  restricted_vcov <- r_matrix %*% vcov_hc1 %*% t(r_matrix)
  wald_chi2 <- as.numeric(t(restrictions) %*% solve(restricted_vcov) %*% restrictions)
  df <- length(terms)
  data.frame(
    test = test_name,
    terms = paste(terms, collapse = " ; "),
    statistic = wald_chi2,
    df = df,
    p_value = pchisq(wald_chi2, df = df, lower.tail = FALSE),
    covariance = "HC1",
    reference_distribution = "Chi-square",
    row.names = NULL
  )
}

ols_table <- function(model, model_name) {
  coef_mat <- summary(model)$coefficients
  data.frame(
    model = model_name,
    variable = rownames(coef_mat),
    coefficient = coef_mat[, "Estimate"],
    std_error = coef_mat[, "Std. Error"],
    t_value = coef_mat[, "t value"],
    p_value = coef_mat[, "Pr(>|t|)"],
    row.names = NULL
  )
}

model_stats <- function(model, model_name) {
  s <- summary(model)
  f <- s$fstatistic
  data.frame(
    model = model_name,
    N = nobs(model),
    R_squared = s$r.squared,
    Adj_R_squared = s$adj.r.squared,
    Residual_SE = s$sigma,
    F_statistic = unname(f["value"]),
    F_p_value = pf(f["value"], f["numdf"], f["dendf"], lower.tail = FALSE),
    row.names = NULL
  )
}

vif_table <- function(model) {
  model_frame <- model.frame(model)
  term_labels <- attr(terms(model), "term.labels")
  simple_terms <- term_labels[!grepl(":", term_labels) & !grepl("\\(", term_labels)]
  out <- lapply(simple_terms, function(term) {
    aux_formula <- as.formula(paste(term, "~", paste(setdiff(simple_terms, term), collapse = " + ")))
    aux <- lm(aux_formula, data = model_frame)
    r2 <- summary(aux)$r.squared
    data.frame(
      variable = term,
      tolerance = 1 - r2,
      VIF = 1 / (1 - r2),
      row.names = NULL
    )
  })
  do.call(rbind, out)
}

nested_f_test <- function(restricted, unrestricted, test_name) {
  a <- anova(restricted, unrestricted)
  data.frame(
    test = test_name,
    df_restricted = a$Res.Df[1],
    df_unrestricted = a$Res.Df[2],
    RSS_restricted = a$RSS[1],
    RSS_unrestricted = a$RSS[2],
    df_num = a$Df[2],
    F_statistic = a$F[2],
    p_value = a$`Pr(>F)`[2],
    row.names = NULL
  )
}

industry_group_data <- function(df) {
  tech <- c("BusSv", "Hardw", "Softw", "Chips", "LabEq", "Telcm", "ElcEq")
  finance <- c("Banks", "Insur", "RlEst", "Fin")
  defensive <- c("Food", "Soda", "Beer", "Smoke", "Hshld", "Hlth", "MedEq", "Drugs", "Util")
  df$TechGrowth <- as.integer(df$industry %in% tech)
  df$Finance <- as.integer(df$industry %in% finance)
  df$Defensive <- as.integer(df$industry %in% defensive)
  df
}

stop_if_missing_inputs()

ff49 <- parse_ff49_daily(ff49_zip)
ff5 <- parse_ff5_daily(ff5_zip)
dgs2 <- read_rate_data()

event_index <- match(event_date, ff49$date)
if (is.na(event_index)) stop("Event date is not in Fama-French 49 industry daily data: ", event_date)

event_window_index <- (event_index - 1):(event_index + 1)
estimation_index <- (event_index - estimation_start_lag):(event_index - estimation_end_lag)
pre_momentum_index <- (event_index - pre_momentum_start_lag):(event_index - pre_momentum_end_lag)

industry_cols <- setdiff(names(ff49), "date")
event_returns <- ff49[event_window_index, c("date", industry_cols)]
estimation_returns <- ff49[estimation_index, c("date", industry_cols)]
pre_returns <- ff49[pre_momentum_index, c("date", industry_cols)]

estimation_data <- merge(estimation_returns, ff5, by = "date", all.x = TRUE)
estimation_data <- merge(estimation_data, dgs2[, c("date", "d_DGS2")], by = "date", all.x = TRUE)
estimation_data <- estimation_data[complete.cases(estimation_data[, c("MktRF", "SMB", "HML", "RMW", "CMA", "RF", "d_DGS2")]), ]

results <- data.frame(
  industry = industry_cols,
  CAR = NA_real_,
  RateVulnerability = NA_real_,
  MktBeta = NA_real_,
  SMBBeta = NA_real_,
  HMLBeta = NA_real_,
  RMWBeta = NA_real_,
  CMABeta = NA_real_,
  HistVol = NA_real_,
  PreMomentum = NA_real_,
  n_estimation = NA_integer_
)

for (industry in industry_cols) {
  reg_data <- estimation_data[, c(industry, "RF", "MktRF", "SMB", "HML", "RMW", "CMA", "d_DGS2")]
  names(reg_data)[1] <- "ret"
  reg_data$excess_ret <- reg_data$ret - reg_data$RF
  reg_data <- reg_data[complete.cases(reg_data), ]
  fit <- lm(excess_ret ~ MktRF + SMB + HML + RMW + CMA + d_DGS2, data = reg_data)

  row <- match(industry, results$industry)
  results$CAR[row] <- sum(event_returns[[industry]], na.rm = TRUE)
  results$RateVulnerability[row] <- -unname(coef(fit)["d_DGS2"])
  results$MktBeta[row] <- unname(coef(fit)["MktRF"])
  results$SMBBeta[row] <- unname(coef(fit)["SMB"])
  results$HMLBeta[row] <- unname(coef(fit)["HML"])
  results$RMWBeta[row] <- unname(coef(fit)["RMW"])
  results$CMABeta[row] <- unname(coef(fit)["CMA"])
  results$HistVol[row] <- sd(reg_data$ret)
  results$PreMomentum[row] <- sum(pre_returns[[industry]], na.rm = TRUE)
  results$n_estimation[row] <- nrow(reg_data)
}

results <- industry_group_data(results)

vars_to_standardize <- c(
  "CAR", "RateVulnerability", "MktBeta", "SMBBeta", "HMLBeta",
  "RMWBeta", "CMABeta", "HistVol", "PreMomentum"
)
for (var in vars_to_standardize) results[[paste0(var, "_z")]] <- zscore(results[[var]])

m0 <- lm(CAR ~ RateVulnerability_z, data = results)
m1 <- lm(CAR ~ RateVulnerability_z + MktBeta_z + HistVol_z + PreMomentum_z, data = results)
m2 <- lm(
  CAR ~ RateVulnerability_z + MktBeta_z + SMBBeta_z + HMLBeta_z +
    RMWBeta_z + CMABeta_z + HistVol_z + PreMomentum_z,
  data = results
)
m_style_reduced <- lm(
  CAR ~ RateVulnerability_z + MktBeta_z + HistVol_z + PreMomentum_z,
  data = results
)
m_interact_reduced <- lm(
  CAR ~ RateVulnerability_z + TechGrowth + Finance + Defensive +
    MktBeta_z + HMLBeta_z + HistVol_z + PreMomentum_z,
  data = results
)
m_interact <- lm(
  CAR ~ RateVulnerability_z * TechGrowth +
    RateVulnerability_z * Finance +
    RateVulnerability_z * Defensive +
    MktBeta_z + HMLBeta_z + HistVol_z + PreMomentum_z,
  data = results
)
m3_full <- lm(
  CAR ~ RateVulnerability_z +
    MktBeta_z + SMBBeta_z + HMLBeta_z + RMWBeta_z + CMABeta_z +
    HistVol_z + PreMomentum_z +
    TechGrowth + Finance + Defensive +
    RateVulnerability_z:TechGrowth +
    RateVulnerability_z:Finance +
    RateVulnerability_z:Defensive,
  data = results
)
m3_full_reduced <- lm(
  CAR ~ RateVulnerability_z +
    MktBeta_z + SMBBeta_z + HMLBeta_z + RMWBeta_z + CMABeta_z +
    HistVol_z + PreMomentum_z +
    TechGrowth + Finance + Defensive,
  data = results
)
m_quad <- lm(
  CAR ~ RateVulnerability_z + I(RateVulnerability_z^2) + MktBeta_z + SMBBeta_z +
    HMLBeta_z + RMWBeta_z + CMABeta_z + HistVol_z + PreMomentum_z,
  data = results
)

fitted_m2 <- fitted(m2)
results$fitted_m2 <- fitted_m2
m_reset <- lm(
  CAR ~ RateVulnerability_z + MktBeta_z + SMBBeta_z + HMLBeta_z +
    RMWBeta_z + CMABeta_z + HistVol_z + PreMomentum_z +
    I(fitted_m2^2) + I(fitted_m2^3),
  data = results
)

m_std <- lm(
  CAR_z ~ RateVulnerability_z + MktBeta_z + SMBBeta_z + HMLBeta_z +
    RMWBeta_z + CMABeta_z + HistVol_z + PreMomentum_z,
  data = results
)

control_terms <- "MktBeta_z + SMBBeta_z + HMLBeta_z + RMWBeta_z + CMABeta_z + HistVol_z + PreMomentum_z"
rv_aux <- lm(as.formula(paste("RateVulnerability_z ~", control_terms)), data = results)
car_aux <- lm(as.formula(paste("CAR ~", control_terms)), data = results)
results$rv_resid <- resid(rv_aux)
results$car_resid <- resid(car_aux)
fwl_model <- lm(car_resid ~ 0 + rv_resid, data = results)
fwl_check <- data.frame(
  full_model_coefficient = unname(coef(m2)["RateVulnerability_z"]),
  fwl_coefficient = unname(coef(fwl_model)["rv_resid"]),
  absolute_difference = abs(unname(coef(m2)["RateVulnerability_z"]) - unname(coef(fwl_model)["rv_resid"]))
)

regression_table <- do.call(rbind, list(
  ols_table(m0, "M0_simple"),
  ols_table(m1, "M1_basic_controls"),
  ols_table(m2, "M2_full_controls"),
  ols_table(m_interact, "M3_interactions_exploratory"),
  ols_table(m3_full, "M3_full_interactions")
))
robust_table <- do.call(rbind, list(
  hc1_table(m0, "M0_simple"),
  hc1_table(m1, "M1_basic_controls"),
  hc1_table(m2, "M2_full_controls"),
  hc1_table(m_interact, "M3_interactions_exploratory"),
  hc1_table(m3_full, "M3_full_interactions")
))
stats_table <- do.call(rbind, list(
  model_stats(m0, "M0_simple"),
  model_stats(m1, "M1_basic_controls"),
  model_stats(m2, "M2_full_controls"),
  model_stats(m_interact, "M3_interactions_exploratory"),
  model_stats(m3_full, "M3_full_interactions")
))

f_tests <- do.call(rbind, list(
  nested_f_test(m0, m1, "M0 vs M1: basic controls jointly zero"),
  nested_f_test(m0, m2, "M0 vs M2: all controls jointly zero"),
  nested_f_test(m_style_reduced, m2, "FF style betas jointly zero"),
  nested_f_test(m_interact_reduced, m_interact, "Industry-group slope interactions jointly zero"),
  nested_f_test(m2, m3_full, "M2 vs M3_full: groups and interactions jointly zero"),
  nested_f_test(m3_full_reduced, m3_full, "M3_full: slope interactions jointly zero"),
  nested_f_test(m2, m_quad, "Rate vulnerability quadratic term zero"),
  nested_f_test(m2, m_reset, "RESET fitted^2 and fitted^3 jointly zero")
))

m2_loo <- do.call(
  rbind,
  lapply(results$industry, function(industry) {
    temp_data <- results[results$industry != industry, ]
    temp_model <- lm(
      CAR ~ RateVulnerability_z + MktBeta_z + SMBBeta_z + HMLBeta_z +
        RMWBeta_z + CMABeta_z + HistVol_z + PreMomentum_z,
      data = temp_data
    )
    temp_coef <- summary(temp_model)$coefficients
    data.frame(
      dropped_industry = industry,
      coefficient = unname(coef(temp_model)["RateVulnerability_z"]),
      std_error = temp_coef["RateVulnerability_z", "Std. Error"],
      t_value = temp_coef["RateVulnerability_z", "t value"],
      p_value = temp_coef["RateVulnerability_z", "Pr(>|t|)"],
      R_squared = summary(temp_model)$r.squared,
      Adj_R_squared = summary(temp_model)$adj.r.squared,
      row.names = NULL
    )
  })
)

coef_path <- data.frame(
  model = c("M0_simple", "M1_basic_controls", "M2_full_controls"),
  specification = c(
    "RateVulnerability only",
    "Add market beta, historical volatility, pre-momentum",
    "Add Fama-French style betas"
  ),
  coefficient = c(
    unname(coef(m0)["RateVulnerability_z"]),
    unname(coef(m1)["RateVulnerability_z"]),
    unname(coef(m2)["RateVulnerability_z"])
  ),
  p_value = c(
    summary(m0)$coefficients["RateVulnerability_z", "Pr(>|t|)"],
    summary(m1)$coefficients["RateVulnerability_z", "Pr(>|t|)"],
    summary(m2)$coefficients["RateVulnerability_z", "Pr(>|t|)"]
  ),
  R_squared = c(summary(m0)$r.squared, summary(m1)$r.squared, summary(m2)$r.squared)
)

interaction_terms <- c(
  "RateVulnerability_z:TechGrowth",
  "RateVulnerability_z:Finance",
  "RateVulnerability_z:Defensive"
)
robust_wald_interaction <- hc1_wald_test(
  m3_full,
  interaction_terms,
  "M3_full HC1 robust Wald: industry-group slope interactions jointly zero"
)

standardized_coefficients <- ols_table(m_std, "standardized_y_and_x")
standardized_coefficients <- standardized_coefficients[standardized_coefficients$variable != "(Intercept)", ]
standardized_coefficients$abs_coefficient <- abs(standardized_coefficients$coefficient)
standardized_coefficients <- standardized_coefficients[order(-standardized_coefficients$abs_coefficient), ]

write.csv(results, file.path(out_dir, "step3_cross_section_data.csv"), row.names = FALSE)
write.csv(regression_table, file.path(out_dir, "step3_regression_table.csv"), row.names = FALSE)
write.csv(robust_table, file.path(out_dir, "step3_regression_table_HC1.csv"), row.names = FALSE)
write.csv(ols_table(m3_full, "M3_full_interactions"), file.path(out_dir, "step3_regression_table_M3_full.csv"), row.names = FALSE)
write.csv(hc1_table(m3_full, "M3_full_interactions"), file.path(out_dir, "step3_regression_table_M3_full_HC1.csv"), row.names = FALSE)
write.csv(stats_table, file.path(out_dir, "step3_model_stats.csv"), row.names = FALSE)
write.csv(vif_table(m2), file.path(out_dir, "step3_vif.csv"), row.names = FALSE)
write.csv(f_tests, file.path(out_dir, "step3_f_tests.csv"), row.names = FALSE)
write.csv(standardized_coefficients, file.path(out_dir, "step3_standardized_coefficients.csv"), row.names = FALSE)
write.csv(fwl_check, file.path(out_dir, "step3_fwl_check.csv"), row.names = FALSE)
write.csv(m2_loo, file.path(out_dir, "step3_leave_one_out_M2.csv"), row.names = FALSE)
write.csv(coef_path, file.path(out_dir, "step3_coef_path.csv"), row.names = FALSE)
write.csv(robust_wald_interaction, file.path(out_dir, "step3_robust_wald_interaction.csv"), row.names = FALSE)

png(file.path(out_dir, "fig_scatter_partial_effect.png"), width = 1200, height = 850, res = 150)
plot(
  results$rv_resid,
  results$car_resid,
  pch = 19,
  col = "#1f77b4",
  xlab = "Residualized rate vulnerability",
  ylab = "Residualized CAR, percent",
  main = "FWL partial effect: rate vulnerability and CAR"
)
abline(fwl_model, col = "#d62728", lwd = 2)
grid(col = "gray85")
text(results$rv_resid, results$car_resid, labels = results$industry, pos = 3, cex = 0.55, col = "gray30")
dev.off()

png(file.path(out_dir, "fig_standardized_coefficients.png"), width = 1200, height = 850, res = 150)
bar_cols <- ifelse(standardized_coefficients$coefficient < 0, "#d62728", "#1f77b4")
coef_values <- standardized_coefficients$coefficient
y_padding <- 0.16 * diff(range(coef_values))
y_limits <- c(min(coef_values) - y_padding, max(coef_values) + y_padding)
label_offset <- 0.04 * diff(y_limits)
label_y <- ifelse(coef_values >= 0, coef_values + label_offset, coef_values - label_offset)
op <- par(mar = c(8, 5, 4, 2) + 0.1)
bp <- barplot(
  coef_values,
  names.arg = standardized_coefficients$variable,
  las = 2,
  col = bar_cols,
  ylim = y_limits,
  ylab = "Standardized coefficient",
  main = "Relative variable importance in the full model"
)
abline(h = 0, col = "gray30")
grid(nx = NA, ny = NULL, col = "gray85")
text(bp, label_y, labels = sprintf("%.2f", coef_values), cex = 0.75)
par(op)
dev.off()

interaction_grid <- expand.grid(
  RateVulnerability_z = seq(min(results$RateVulnerability_z), max(results$RateVulnerability_z), length.out = 100),
  group = c("Other", "TechGrowth", "Finance", "Defensive")
)
interaction_grid$TechGrowth <- as.integer(interaction_grid$group == "TechGrowth")
interaction_grid$Finance <- as.integer(interaction_grid$group == "Finance")
interaction_grid$Defensive <- as.integer(interaction_grid$group == "Defensive")
interaction_grid$MktBeta_z <- 0
interaction_grid$SMBBeta_z <- 0
interaction_grid$HMLBeta_z <- 0
interaction_grid$RMWBeta_z <- 0
interaction_grid$CMABeta_z <- 0
interaction_grid$HistVol_z <- 0
interaction_grid$PreMomentum_z <- 0
interaction_grid$predicted_CAR <- predict(m3_full, newdata = interaction_grid)

png(file.path(out_dir, "fig_interaction_effect.png"), width = 1200, height = 850, res = 150)
plot(
  range(interaction_grid$RateVulnerability_z),
  range(interaction_grid$predicted_CAR),
  type = "n",
  xlab = "Rate vulnerability (standardized)",
  ylab = "Predicted CAR, percent",
  main = "Industry-group interaction effects"
)
group_cols <- c(Other = "gray25", TechGrowth = "#1f77b4", Finance = "#2ca02c", Defensive = "#ff7f0e")
for (grp in names(group_cols)) {
  temp <- interaction_grid[interaction_grid$group == grp, ]
  lines(temp$RateVulnerability_z, temp$predicted_CAR, col = group_cols[grp], lwd = 2)
}
grid(col = "gray85")
legend("topright", legend = names(group_cols), col = group_cols, lwd = 2, bty = "n")
dev.off()

cat("Done.\n")
cat("Event date:", format(event_date), "\n")
cat("Event window dates:", paste(format(event_returns$date), collapse = ", "), "\n")
cat("Estimation rows after merge:", nrow(estimation_data), "\n")
cat("Output directory:", out_dir, "\n\n")
cat("Main coefficient, M2 full controls:\n")
print(regression_table[regression_table$model == "M2_full_controls" & regression_table$variable == "RateVulnerability_z", ], row.names = FALSE)
cat("\nHC1 robust coefficient, M2 full controls:\n")
print(robust_table[robust_table$model == "M2_full_controls" & robust_table$variable == "RateVulnerability_z", ], row.names = FALSE)
cat("\nFWL check:\n")
print(fwl_check, row.names = FALSE)
cat("\nF tests:\n")
print(f_tests, row.names = FALSE)
cat("\nM2 leave-one-out summary:\n")
print(
  data.frame(
    negative_coefficients = sum(m2_loo$coefficient < 0),
    significant_at_5pct = sum(m2_loo$p_value < 0.05),
    significant_at_10pct = sum(m2_loo$p_value < 0.10),
    max_p_value = max(m2_loo$p_value),
    weakest_drop = m2_loo$dropped_industry[which.max(m2_loo$p_value)],
    total = nrow(m2_loo)
  ),
  row.names = FALSE
)
cat("\nM3_full robust Wald interaction test:\n")
print(robust_wald_interaction, row.names = FALSE)
