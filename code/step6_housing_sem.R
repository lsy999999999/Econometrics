options(stringsAsFactors = FALSE)

this_file <- tryCatch(normalizePath(sys.frames()[[1]]$ofile), error = function(e) NA_character_)
if (is.na(this_file)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) this_file <- normalizePath(sub("^--file=", "", file_arg[1]))
}

project_root <- normalizePath(file.path(dirname(this_file), ".."))
raw_dir <- file.path(project_root, "data", "raw", "fred_housing")
out_dir <- file.path(project_root, "output", "step6_housing")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 180)

series_info <- data.frame(
  series_id = c(
    "HSN1F", "MSPNHSUS", "MORTGAGE30US", "DSPIC96",
    "UNRATE", "PERMIT", "MSACSR", "WPUSI012011"
  ),
  variable = c(
    "sales", "price", "mortgage_rate", "income",
    "unrate", "permit", "months_supply", "construction_cost"
  ),
  description = c(
    "New one-family houses sold, United States",
    "Median sales price for new houses sold, United States",
    "30-year fixed mortgage rate, weekly average",
    "Real disposable personal income",
    "Unemployment rate",
    "New privately owned housing units authorized by building permits",
    "Monthly supply of new houses",
    "PPI special index: construction materials"
  ),
  frequency = c("monthly", "monthly", "weekly_to_monthly_mean", "monthly", "monthly", "monthly", "monthly", "monthly"),
  transformation = c("log", "log", "monthly mean, level", "log", "level", "log", "level", "log"),
  equation_role = c(
    "endogenous quantity", "endogenous price", "demand shifter",
    "demand shifter", "demand shifter", "supply shifter",
    "supply shifter", "supply shifter"
  ),
  source_url = paste0("https://fred.stlouisfed.org/series/", c(
    "HSN1F", "MSPNHSUS", "MORTGAGE30US", "DSPIC96",
    "UNRATE", "PERMIT", "MSACSR", "WPUSI012011"
  )),
  row.names = NULL
)
write.csv(series_info, file.path(out_dir, "housing_variable_dictionary.csv"), row.names = FALSE)

download_fred <- function(series_id) {
  raw_file <- file.path(raw_dir, paste0(series_id, ".csv"))
  needs_download <- !file.exists(raw_file)
  if (!needs_download) {
    cached <- tryCatch(read.csv(raw_file, na.strings = c(".", "", "NA")), error = function(e) NULL)
    if (is.null(cached) || nrow(cached) < 100) {
      needs_download <- TRUE
    } else {
      cached_date <- as.Date(cached[[1]])
      needs_download <- min(cached_date, na.rm = TRUE) > as.Date("2000-01-01") ||
        max(cached_date, na.rm = TRUE) < as.Date("2024-12-01")
    }
  }
  if (needs_download) {
    url <- paste0(
      "https://fred.stlouisfed.org/graph/fredgraph.csv?cosd=2000-01-01&coed=2024-12-31&id=",
      series_id
    )
    download.file(url, raw_file, mode = "wb", quiet = TRUE, method = "libcurl")
  }
  dat <- read.csv(raw_file, na.strings = c(".", "", "NA"))
  names(dat)[1] <- "date"
  dat$date <- as.Date(dat$date)
  names(dat)[2] <- series_id
  dat[[series_id]] <- suppressWarnings(as.numeric(dat[[series_id]]))
  dat
}

to_month <- function(x) as.Date(paste0(format(x, "%Y-%m"), "-01"))

monthly_mean <- function(dat, value_col) {
  dat$month <- to_month(dat$date)
  aggregate(dat[[value_col]], list(date = dat$month), function(x) mean(x, na.rm = TRUE))
}

series_list <- lapply(series_info$series_id, download_fred)
names(series_list) <- series_info$series_id

monthly_list <- list()
for (sid in names(series_list)) {
  dat <- series_list[[sid]]
  if (sid == "MORTGAGE30US") {
    m <- monthly_mean(dat, sid)
    names(m)[2] <- sid
  } else {
    dat$date <- to_month(dat$date)
    m <- dat[, c("date", sid)]
  }
  monthly_list[[sid]] <- m
}

raw_monthly <- Reduce(function(x, y) merge(x, y, by = "date", all = TRUE), monthly_list)
for (i in seq_len(nrow(series_info))) {
  names(raw_monthly)[names(raw_monthly) == series_info$series_id[i]] <- series_info$variable[i]
}
raw_monthly <- raw_monthly[order(raw_monthly$date), ]
write.csv(raw_monthly, file.path(out_dir, "housing_sem_raw.csv"), row.names = FALSE)

df <- raw_monthly[raw_monthly$date >= as.Date("2000-01-01") & raw_monthly$date <= as.Date("2024-12-01"), ]
df$log_sales <- log(df$sales)
df$log_price <- log(df$price)
df$log_income <- log(df$income)
df$log_permit <- log(df$permit)
df$log_construction_cost <- log(df$construction_cost)
df$month <- factor(format(df$date, "%m"))
df <- df[complete.cases(df[, c(
  "log_sales", "log_price", "mortgage_rate", "log_income", "unrate",
  "log_permit", "months_supply", "log_construction_cost"
)]), ]
write.csv(df, file.path(out_dir, "housing_sem_clean.csv"), row.names = FALSE)

hc_vcov <- function(model, type = "HC3") {
  x <- model.matrix(model)
  residual <- residuals(model)
  bread <- MASS::ginv(crossprod(x))
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

coef_table <- function(model, equation, method, vcov_matrix = NULL) {
  if (is.null(vcov_matrix)) vcov_matrix <- hc_vcov(model, "HC3")
  est <- coef(model)
  se <- sqrt(diag(vcov_matrix))
  stat <- est / se
  p <- 2 * pt(abs(stat), df = max(df.residual(model), 1), lower.tail = FALSE)
  data.frame(
    equation = equation,
    method = method,
    term = names(est),
    estimate = unname(est),
    standard_error = unname(se),
    statistic = unname(stat),
    p_value = unname(p),
    N = nobs(model),
    row.names = NULL
  )
}

first_stage_relevance <- function(data, equation, endogenous_var, included_exog, excluded_exog) {
  full <- lm(as.formula(paste(endogenous_var, "~", paste(c(included_exog, excluded_exog), collapse = " + "))), data = data)
  restricted <- lm(as.formula(paste(endogenous_var, "~", paste(included_exog, collapse = " + "))), data = data)
  a <- anova(restricted, full)
  partial_r2 <- (sum(resid(restricted)^2) - sum(resid(full)^2)) / sum(resid(restricted)^2)
  data.frame(
    equation = equation,
    endogenous_rhs_variable = endogenous_var,
    excluded_IVs = paste(excluded_exog, collapse = "; "),
    first_stage_F = a$F[2],
    first_stage_p = a$`Pr(>F)`[2],
    partial_R2 = partial_r2,
    weak_IV_flag_F_below_10 = a$F[2] < 10,
    N = nobs(full),
    row.names = NULL
  )
}

manual_2sls <- function(data, y, x_endog, included_exog, excluded_iv) {
  rhs_x <- c(x_endog, included_exog)
  rhs_z <- c(included_exog, excluded_iv)
  y_vec <- as.matrix(data[, y, drop = FALSE])
  x_matrix <- model.matrix(as.formula(paste("~", paste(rhs_x, collapse = " + "))), data = data)
  z_matrix <- model.matrix(as.formula(paste("~", paste(rhs_z, collapse = " + "))), data = data)
  pz <- z_matrix %*% MASS::ginv(crossprod(z_matrix)) %*% t(z_matrix)
  x_hat <- pz %*% x_matrix
  bread <- MASS::ginv(t(x_hat) %*% x_matrix)
  beta <- bread %*% t(x_hat) %*% y_vec
  rownames(beta) <- colnames(x_matrix)
  residual <- as.vector(y_vec - x_matrix %*% beta)
  n <- nrow(x_matrix)
  k <- ncol(x_matrix)
  hat_diag <- diag(x_hat %*% bread %*% t(x_hat))
  adjustment <- 1 / pmax(1 - hat_diag, .Machine$double.eps)^2
  meat <- t(x_hat) %*% ((residual^2 * adjustment) * x_hat)
  vcov_hc3 <- bread %*% meat %*% t(bread)
  se <- sqrt(diag(vcov_hc3))
  stat <- as.vector(beta) / se
  p <- 2 * pt(abs(stat), df = max(n - k, 1), lower.tail = FALSE)
  list(
    coefficients = as.vector(beta),
    names = rownames(beta),
    se_hc3 = se,
    statistic = stat,
    p_value = p,
    residuals = residual,
    fitted = as.vector(x_matrix %*% beta),
    x_matrix = x_matrix,
    z_matrix = z_matrix,
    n = n,
    k = k,
    y = y,
    x_endog = x_endog,
    included_exog = included_exog,
    excluded_iv = excluded_iv
  )
}

iv_table <- function(fit, equation, method) {
  data.frame(
    equation = equation,
    method = method,
    term = fit$names,
    estimate = fit$coefficients,
    standard_error = fit$se_hc3,
    statistic = fit$statistic,
    p_value = fit$p_value,
    N = fit$n,
    row.names = NULL
  )
}

dwh_test <- function(data, equation, y, x_endog, included_exog, excluded_iv) {
  first <- lm(as.formula(paste(x_endog, "~", paste(c(included_exog, excluded_iv), collapse = " + "))), data = data)
  data$first_stage_resid <- resid(first)
  augmented <- lm(as.formula(paste(y, "~", paste(c(x_endog, included_exog, "first_stage_resid"), collapse = " + "))), data = data)
  ct <- coef_table(augmented, equation, "DWH_residual_inclusion")
  row <- ct[ct$term == "first_stage_resid", ]
  data.frame(
    equation = equation,
    test = "Durbin-Wu-Hausman residual inclusion",
    endogenous_rhs_variable = x_endog,
    residual_coefficient = row$estimate,
    p_value = row$p_value,
    reject_exogeneity_5pct = row$p_value < 0.05,
    N = nobs(augmented),
    row.names = NULL
  )
}

sargan_test <- function(fit, equation, label) {
  overid_df <- ncol(fit$z_matrix) - ncol(fit$x_matrix)
  if (overid_df <= 0) {
    return(data.frame(equation = equation, IV_set = label, test_type = "Sargan", statistic = NA_real_, df = overid_df, p_value = NA_real_, reject_at_5pct = NA, row.names = NULL))
  }
  aux <- lm(fit$residuals ~ fit$z_matrix[, -1, drop = FALSE])
  stat <- fit$n * summary(aux)$r.squared
  data.frame(equation = equation, IV_set = label, test_type = "Sargan", statistic = stat, df = overid_df, p_value = pchisq(stat, overid_df, lower.tail = FALSE), reject_at_5pct = pchisq(stat, overid_df, lower.tail = FALSE) < 0.05, row.names = NULL)
}

hansen_j_test <- function(fit, equation, label) {
  overid_df <- ncol(fit$z_matrix) - ncol(fit$x_matrix)
  if (overid_df <= 0) {
    return(data.frame(equation = equation, IV_set = label, test_type = "Hansen_style_J", statistic = NA_real_, df = overid_df, p_value = NA_real_, reject_at_5pct = NA, row.names = NULL))
  }
  z <- fit$z_matrix
  u <- fit$residuals
  n <- nrow(z)
  moments <- z * as.numeric(u)
  gbar <- colMeans(moments)
  s_matrix <- crossprod(scale(moments, center = FALSE, scale = FALSE)) / n
  stat <- as.numeric(n * t(gbar) %*% MASS::ginv(s_matrix) %*% gbar)
  data.frame(equation = equation, IV_set = label, test_type = "Hansen_style_J", statistic = stat, df = overid_df, p_value = pchisq(stat, overid_df, lower.tail = FALSE), reject_at_5pct = pchisq(stat, overid_df, lower.tail = FALSE) < 0.05, row.names = NULL)
}

descriptive_stats <- data.frame(
  variable = names(df)[sapply(df, is.numeric)],
  N = sapply(df[sapply(df, is.numeric)], function(x) sum(!is.na(x))),
  mean = sapply(df[sapply(df, is.numeric)], mean, na.rm = TRUE),
  sd = sapply(df[sapply(df, is.numeric)], sd, na.rm = TRUE),
  min = sapply(df[sapply(df, is.numeric)], min, na.rm = TRUE),
  p25 = sapply(df[sapply(df, is.numeric)], quantile, probs = 0.25, na.rm = TRUE),
  median = sapply(df[sapply(df, is.numeric)], median, na.rm = TRUE),
  p75 = sapply(df[sapply(df, is.numeric)], quantile, probs = 0.75, na.rm = TRUE),
  max = sapply(df[sapply(df, is.numeric)], max, na.rm = TRUE),
  row.names = NULL
)
write.csv(descriptive_stats, file.path(out_dir, "housing_descriptive_stats.csv"), row.names = FALSE)

cor_vars <- c("log_sales", "log_price", "mortgage_rate", "log_income", "unrate", "log_permit", "months_supply", "log_construction_cost")
cor_mat <- round(cor(df[, cor_vars], use = "pairwise.complete.obs"), 4)
write.csv(cor_mat, file.path(out_dir, "housing_correlation_matrix.csv"))

png(file.path(out_dir, "housing_time_series_price_quantity.png"), width = 1400, height = 900, res = 140)
par(mar = c(5, 5, 4, 5))
plot(df$date, df$log_price, type = "l", col = "#1f77b4", lwd = 2, xlab = "Date", ylab = "log price", main = "New home price and sales")
par(new = TRUE)
plot(df$date, df$log_sales, type = "l", col = "#d62728", lwd = 2, axes = FALSE, xlab = "", ylab = "")
axis(4)
mtext("log sales", side = 4, line = 3)
legend("topleft", legend = c("log price", "log sales"), col = c("#1f77b4", "#d62728"), lwd = 2, bty = "n")
dev.off()

png(file.path(out_dir, "housing_time_series_rate_sales.png"), width = 1400, height = 900, res = 140)
par(mar = c(5, 5, 4, 5))
plot(df$date, df$mortgage_rate, type = "l", col = "#2ca02c", lwd = 2, xlab = "Date", ylab = "mortgage rate", main = "Mortgage rate and new home sales")
par(new = TRUE)
plot(df$date, df$log_sales, type = "l", col = "#d62728", lwd = 2, axes = FALSE, xlab = "", ylab = "")
axis(4)
mtext("log sales", side = 4, line = 3)
legend("topright", legend = c("mortgage rate", "log sales"), col = c("#2ca02c", "#d62728"), lwd = 2, bty = "n")
dev.off()

demand_included <- c("mortgage_rate", "log_income", "unrate")
demand_excluded <- c("log_permit", "months_supply", "log_construction_cost")
supply_included <- c("log_permit", "months_supply", "log_construction_cost")
supply_excluded <- c("mortgage_rate", "log_income", "unrate")

demand_ols <- lm(log_sales ~ log_price + mortgage_rate + log_income + unrate, data = df)
supply_ols <- lm(log_sales ~ log_price + log_permit + months_supply + log_construction_cost, data = df)
write.csv(coef_table(demand_ols, "demand", "OLS_HC3"), file.path(out_dir, "housing_ols_demand_HC3.csv"), row.names = FALSE)
write.csv(coef_table(supply_ols, "supply", "OLS_HC3"), file.path(out_dir, "housing_ols_supply_HC3.csv"), row.names = FALSE)

first_stage <- rbind(
  first_stage_relevance(df, "demand", "log_price", demand_included, demand_excluded),
  first_stage_relevance(df, "supply", "log_price", supply_included, supply_excluded)
)
write.csv(first_stage, file.path(out_dir, "housing_first_stage_relevance.csv"), row.names = FALSE)

demand_2sls <- manual_2sls(df, "log_sales", "log_price", demand_included, demand_excluded)
supply_2sls <- manual_2sls(df, "log_sales", "log_price", supply_included, supply_excluded)
demand_2sls_tab <- iv_table(demand_2sls, "demand", "2SLS_manual_HC3")
supply_2sls_tab <- iv_table(supply_2sls, "supply", "2SLS_manual_HC3")
write.csv(demand_2sls_tab, file.path(out_dir, "housing_2sls_demand_HC3.csv"), row.names = FALSE)
write.csv(supply_2sls_tab, file.path(out_dir, "housing_2sls_supply_HC3.csv"), row.names = FALSE)

overid_tests <- rbind(
  sargan_test(demand_2sls, "demand", "supply_shifters"),
  hansen_j_test(demand_2sls, "demand", "supply_shifters"),
  sargan_test(supply_2sls, "supply", "demand_shifters"),
  hansen_j_test(supply_2sls, "supply", "demand_shifters")
)
write.csv(overid_tests, file.path(out_dir, "housing_overid_tests.csv"), row.names = FALSE)

endogeneity_tests <- rbind(
  dwh_test(df, "demand", "log_sales", "log_price", demand_included, demand_excluded),
  dwh_test(df, "supply", "log_sales", "log_price", supply_included, supply_excluded)
)
write.csv(endogeneity_tests, file.path(out_dir, "housing_endogeneity_tests.csv"), row.names = FALSE)

diagnostics <- rbind(
  transform(first_stage, diagnostic = "first_stage_relevance")[, c("equation", "diagnostic", "first_stage_F", "partial_R2", "weak_IV_flag_F_below_10", "N")],
  data.frame(
    equation = endogeneity_tests$equation,
    diagnostic = "DWH_price_endogeneity",
    first_stage_F = NA_real_,
    partial_R2 = endogeneity_tests$p_value,
    weak_IV_flag_F_below_10 = endogeneity_tests$reject_exogeneity_5pct,
    N = endogeneity_tests$N
  )
)
names(diagnostics)[names(diagnostics) == "partial_R2"] <- "value"
write.csv(diagnostics, file.path(out_dir, "housing_2sls_diagnostics.csv"), row.names = FALSE)

system_results <- data.frame()
system_resid <- data.frame()
if (requireNamespace("systemfit", quietly = TRUE)) {
  eqs <- list(
    demand = log_sales ~ log_price + mortgage_rate + log_income + unrate,
    supply = log_sales ~ log_price + log_permit + months_supply + log_construction_cost
  )
  inst <- ~ mortgage_rate + log_income + unrate + log_permit + months_supply + log_construction_cost
  fit_sys_2sls <- tryCatch(systemfit::systemfit(eqs, method = "2SLS", inst = inst, data = df), error = identity)
  fit_sys_3sls <- tryCatch(systemfit::systemfit(eqs, method = "3SLS", inst = inst, data = df), error = identity)
  if (!inherits(fit_sys_2sls, "error")) {
    writeLines(capture.output(summary(fit_sys_2sls)), file.path(out_dir, "housing_system_2sls_summary.txt"))
    b <- coef(fit_sys_2sls)
    se <- sqrt(diag(vcov(fit_sys_2sls)))
    system_results <- rbind(system_results, data.frame(method = "systemfit_2SLS", term = names(b), estimate = unname(b), standard_error = unname(se), statistic = unname(b / se), p_value = 2 * pnorm(abs(b / se), lower.tail = FALSE), row.names = NULL))
    r <- residuals(fit_sys_2sls)
    ct <- cor.test(r[, 1], r[, 2])
    system_resid <- data.frame(method = "systemfit_2SLS", residual_correlation = unname(ct$estimate), p_value = ct$p.value, row.names = NULL)
  } else {
    writeLines(paste("systemfit 2SLS failed:", fit_sys_2sls$message), file.path(out_dir, "housing_system_2sls_summary.txt"))
  }
  if (!inherits(fit_sys_3sls, "error")) {
    writeLines(capture.output(summary(fit_sys_3sls)), file.path(out_dir, "housing_system_3sls_summary.txt"))
    b <- coef(fit_sys_3sls)
    se <- sqrt(diag(vcov(fit_sys_3sls)))
    system_results <- rbind(system_results, data.frame(method = "systemfit_3SLS", term = names(b), estimate = unname(b), standard_error = unname(se), statistic = unname(b / se), p_value = 2 * pnorm(abs(b / se), lower.tail = FALSE), row.names = NULL))
  } else {
    writeLines(paste("systemfit 3SLS failed:", fit_sys_3sls$message), file.path(out_dir, "housing_system_3sls_summary.txt"))
  }
}
write.csv(system_results, file.path(out_dir, "housing_system_results.csv"), row.names = FALSE)
write.csv(system_results[system_results$method == "systemfit_2SLS", ], file.path(out_dir, "housing_system_2sls_results.csv"), row.names = FALSE)
write.csv(system_results[system_results$method == "systemfit_3SLS", ], file.path(out_dir, "housing_system_3sls_results.csv"), row.names = FALSE)
write.csv(system_resid, file.path(out_dir, "housing_system_residual_correlation.csv"), row.names = FALSE)

newey_west_table <- function(model, equation) {
  vc <- sandwich::NeweyWest(model, lag = 6, prewhite = FALSE, adjust = TRUE)
  coef_table(model, equation, "OLS_NeweyWest_lag6", vc)
}

yoy <- function(x) c(rep(NA_real_, 12), diff(log(x), lag = 12))
df_yoy <- df
df_yoy$yoy_sales <- yoy(df_yoy$sales)
df_yoy$yoy_price <- yoy(df_yoy$price)
df_yoy$yoy_income <- yoy(df_yoy$income)
df_yoy$yoy_permit <- yoy(df_yoy$permit)
df_yoy$yoy_construction_cost <- yoy(df_yoy$construction_cost)
df_yoy <- df_yoy[complete.cases(df_yoy[, c("yoy_sales", "yoy_price", "mortgage_rate", "yoy_income", "unrate", "yoy_permit", "months_supply", "yoy_construction_cost")]), ]

robust_results <- list()
robust_fit <- function(data, check, demand_formula, supply_formula) {
  d <- lm(demand_formula, data = data)
  s <- lm(supply_formula, data = data)
  dtab <- coef_table(d, "demand", check)
  stab <- coef_table(s, "supply", check)
  data.frame(
    check = check,
    demand_price_coef = dtab$estimate[dtab$term %in% c("log_price", "yoy_price")][1],
    demand_price_p = dtab$p_value[dtab$term %in% c("log_price", "yoy_price")][1],
    supply_price_coef = stab$estimate[stab$term %in% c("log_price", "yoy_price")][1],
    supply_price_p = stab$p_value[stab$term %in% c("log_price", "yoy_price")][1],
    demand_sign_ok = dtab$estimate[dtab$term %in% c("log_price", "yoy_price")][1] < 0,
    supply_sign_ok = stab$estimate[stab$term %in% c("log_price", "yoy_price")][1] > 0,
    N = nobs(d),
    note = "OLS robustness check; price remains simultaneous, so interpret as sensitivity evidence.",
    row.names = NULL
  )
}
robust_results[[1]] <- robust_fit(df[df$date <= as.Date("2019-12-01"), ], "pre_covid_2000_2019", log_sales ~ log_price + mortgage_rate + log_income + unrate, log_sales ~ log_price + log_permit + months_supply + log_construction_cost)
robust_results[[2]] <- robust_fit(df[df$date >= as.Date("2010-01-01"), ], "post_gfc_2010_2024", log_sales ~ log_price + mortgage_rate + log_income + unrate, log_sales ~ log_price + log_permit + months_supply + log_construction_cost)
robust_results[[3]] <- robust_fit(df_yoy, "yoy_growth_rates", yoy_sales ~ yoy_price + mortgage_rate + yoy_income + unrate, yoy_sales ~ yoy_price + yoy_permit + months_supply + yoy_construction_cost)
robust_results[[4]] <- robust_fit(df, "month_fixed_effects", log_sales ~ log_price + mortgage_rate + log_income + unrate + month, log_sales ~ log_price + log_permit + months_supply + log_construction_cost + month)
nw_demand <- newey_west_table(demand_ols, "demand")
nw_supply <- newey_west_table(supply_ols, "supply")
robust_results[[5]] <- data.frame(
  check = "newey_west_lag6",
  demand_price_coef = nw_demand$estimate[nw_demand$term == "log_price"],
  demand_price_p = nw_demand$p_value[nw_demand$term == "log_price"],
  supply_price_coef = nw_supply$estimate[nw_supply$term == "log_price"],
  supply_price_p = nw_supply$p_value[nw_supply$term == "log_price"],
  demand_sign_ok = nw_demand$estimate[nw_demand$term == "log_price"] < 0,
  supply_sign_ok = nw_supply$estimate[nw_supply$term == "log_price"] > 0,
  N = nobs(demand_ols),
  note = "Same OLS coefficients with Newey-West lag 6 standard errors.",
  row.names = NULL
)
robustness_summary <- do.call(rbind, robust_results)
write.csv(robustness_summary, file.path(out_dir, "housing_robustness_summary.csv"), row.names = FALSE)

identification <- data.frame(
  equation = c("demand", "supply"),
  dependent_variable = c("log_sales", "log_sales"),
  endogenous_rhs_variable = c("log_price", "log_price"),
  included_exogenous = c(paste(demand_included, collapse = "; "), paste(supply_included, collapse = "; ")),
  excluded_exogenous = c(paste(demand_excluded, collapse = "; "), paste(supply_excluded, collapse = "; ")),
  rhs_endogenous_count_Gi = c(1, 1),
  excluded_exogenous_count = c(length(demand_excluded), length(supply_excluded)),
  order_condition = c("overidentified", "overidentified"),
  economic_logic = c(
    "Supply shifters move the supply curve and identify demand.",
    "Demand shifters move the demand curve and identify supply."
  ),
  row.names = NULL
)
write.csv(identification, file.path(out_dir, "housing_identification_conditions.csv"), row.names = FALSE)

lecture_coverage <- data.frame(
  lecture6_topic = c(
    "双向或联立关系", "结构方程", "结构系数", "结构误差", "内生变量", "外生变量",
    "前定变量", "OLS 不一致", "识别条件", "2SLS", "系统估计", "过度识别检验", "内生性检验", "2SLS 拟合统计量警告"
  ),
  housing_module_mapping = c(
    "房价 log_price 和成交量 log_sales 由需求和供给共同决定。",
    "需求方程和供给方程。",
    "价格弹性、利率系数、收入系数、建筑许可系数等。",
    "需求冲击 u_d 与供给冲击 u_s；系统估计中检查残差相关。",
    "log_sales 与 log_price。",
    "mortgage_rate, log_income, unrate, log_permit, months_supply, log_construction_cost。",
    "月度宏观变量和供给变量可进一步加入滞后项；本模块把需求/供给移动变量作为系统外生变量。",
    "log_price 由供需共同决定，直接 OLS 估计结构方程会有联立性偏误。",
    "每个方程右侧 1 个内生变量，并排除 3 个外生变量，两个方程均过度识别。",
    "手工 2SLS 估计需求方程和供给方程，并输出 HC3 标准误。",
    "systemfit 2SLS 和 3SLS。",
    "Sargan 与 Hansen-style J。",
    "Durbin-Wu-Hausman 残差纳入检验。",
    "不使用普通 2SLS R2/F 作为结论依据。"
  ),
  row.names = NULL
)
write.csv(lecture_coverage, file.path(out_dir, "housing_lecture6_coverage.csv"), row.names = FALSE)

key <- function(tab, term) tab[tab$term == term, ][1, ]
d_ols_price <- key(coef_table(demand_ols, "demand", "OLS_HC3"), "log_price")
s_ols_price <- key(coef_table(supply_ols, "supply", "OLS_HC3"), "log_price")
d_iv_price <- key(demand_2sls_tab, "log_price")
s_iv_price <- key(supply_2sls_tab, "log_price")
d_rate <- key(demand_2sls_tab, "mortgage_rate")
d_income <- key(demand_2sls_tab, "log_income")
d_unrate <- key(demand_2sls_tab, "unrate")
s_permit <- key(supply_2sls_tab, "log_permit")
s_months <- key(supply_2sls_tab, "months_supply")
s_cost <- key(supply_2sls_tab, "log_construction_cost")

fmt <- function(x, digits = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))

summary_lines <- c(
  "# Step 6 Housing SEM Results",
  "",
  "## Research Question",
  "",
  "This module replaces the exploratory FOMC industry SEM with a standard new-housing-market supply-demand simultaneous-equations model. The goal is to demonstrate Lecture 6 with an economically coherent system: new home prices and new home sales are jointly determined by demand and supply.",
  "",
  "## Data",
  "",
  paste0("Monthly U.S. data from FRED are used from ", format(min(df$date), "%Y-%m"), " to ", format(max(df$date), "%Y-%m"), ", with N = ", nrow(df), "."),
  "The key endogenous variables are `log_sales` and `log_price`. `sales` is new one-family houses sold and `price` is the median sales price for new houses sold. Demand shifters are the mortgage rate, real disposable income, and unemployment. Supply shifters are permits, monthly supply of new houses, and construction material costs.",
  "",
  "## Structural Equations",
  "",
  "Demand: `log_sales = alpha0 + alpha1 log_price + alpha2 mortgage_rate + alpha3 log_income + alpha4 unrate + u_d`.",
  "Supply: `log_sales = beta0 + beta1 log_price + beta2 log_permit + beta3 months_supply + beta4 log_construction_cost + u_s`.",
  "Equilibrium condition: observed sales are both demanded and supplied, so `Qd = Qs = Q`.",
  "",
  "## Identification",
  "",
  "Both structural equations contain one endogenous right-hand-side variable, `log_price`. The demand equation excludes three supply shifters; the supply equation excludes three demand shifters. Therefore both equations satisfy the order condition and are overidentified.",
  "",
  "## First Stage",
  "",
  paste0("Demand equation instruments: first-stage F = ", fmt(first_stage$first_stage_F[first_stage$equation == "demand"]), ", partial R2 = ", fmt(first_stage$partial_R2[first_stage$equation == "demand"]), "."),
  paste0("Supply equation instruments: first-stage F = ", fmt(first_stage$first_stage_F[first_stage$equation == "supply"]), ", partial R2 = ", fmt(first_stage$partial_R2[first_stage$equation == "supply"]), "."),
  "",
  "## OLS versus 2SLS",
  "",
  paste0("Demand OLS price coefficient = ", fmt(d_ols_price$estimate), ", p = ", fmt(d_ols_price$p_value), "."),
  paste0("Demand 2SLS price coefficient = ", fmt(d_iv_price$estimate), ", p = ", fmt(d_iv_price$p_value), ". Expected sign: negative."),
  paste0("Supply OLS price coefficient = ", fmt(s_ols_price$estimate), ", p = ", fmt(s_ols_price$p_value), "."),
  paste0("Supply 2SLS price coefficient = ", fmt(s_iv_price$estimate), ", p = ", fmt(s_iv_price$p_value), ". Expected sign: positive."),
  "",
  "## Structural Coefficients",
  "",
  paste0("Demand 2SLS mortgage-rate coefficient = ", fmt(d_rate$estimate), ", p = ", fmt(d_rate$p_value), ". Expected sign: negative."),
  paste0("Demand 2SLS income coefficient = ", fmt(d_income$estimate), ", p = ", fmt(d_income$p_value), ". Expected sign: positive."),
  paste0("Demand 2SLS unemployment coefficient = ", fmt(d_unrate$estimate), ", p = ", fmt(d_unrate$p_value), ". Expected sign: negative."),
  paste0("Supply 2SLS permit coefficient = ", fmt(s_permit$estimate), ", p = ", fmt(s_permit$p_value), ". Expected sign: positive."),
  paste0("Supply 2SLS months-supply coefficient = ", fmt(s_months$estimate), ", p = ", fmt(s_months$p_value), ". Expected sign: positive."),
  paste0("Supply 2SLS construction-cost coefficient = ", fmt(s_cost$estimate), ", p = ", fmt(s_cost$p_value), ". Expected sign: negative."),
  "",
  "## Diagnostics",
  "",
  paste0("DWH demand p-value = ", fmt(endogeneity_tests$p_value[endogeneity_tests$equation == "demand"]), "; DWH supply p-value = ", fmt(endogeneity_tests$p_value[endogeneity_tests$equation == "supply"]), "."),
  paste0("Overidentification rejected at 5% in demand equation: ", any(overid_tests$equation == "demand" & overid_tests$reject_at_5pct %in% TRUE), "."),
  paste0("Overidentification rejected at 5% in supply equation: ", any(overid_tests$equation == "supply" & overid_tests$reject_at_5pct %in% TRUE), "."),
  if (nrow(system_resid) > 0) paste0("System 2SLS residual correlation = ", fmt(system_resid$residual_correlation), ", p = ", fmt(system_resid$p_value), ".") else "System residual correlation was not available.",
  "",
  "## Interpretation",
  "",
  "This housing module is a better Lecture 6 application than the previous FOMC CAR-RV exploratory SEM because price and quantity have a standard simultaneous-equilibrium interpretation. The demand and supply equations each have clear ceteris-paribus meanings, and the excluded demand/supply shifters provide a transparent identification strategy.",
  "",
  "All tables and figures are saved in `output/step6_housing/`."
)
writeLines(summary_lines, file.path(out_dir, "housing_sem_results_summary.md"))

message("Step 6 housing SEM complete. Outputs written to ", out_dir)
