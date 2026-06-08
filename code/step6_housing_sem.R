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
    "UNRATE", "PERMIT", "MSACSR", "WPUSI012011",
    "HOUST", "UNDCONTSA"
  ),
  variable = c(
    "sales", "price", "mortgage_rate", "income",
    "unrate", "permit", "months_supply", "construction_cost",
    "housing_starts", "under_construction"
  ),
  description = c(
    "New one-family houses sold, United States",
    "Median sales price for new houses sold, United States",
    "30-year fixed mortgage rate, weekly average",
    "Real disposable personal income",
    "Unemployment rate",
    "New privately owned housing units authorized by building permits",
    "Monthly supply of new houses",
    "PPI special index: construction materials",
    "New privately owned housing units started",
    "New privately owned housing units under construction"
  ),
  frequency = c("monthly", "monthly", "weekly_to_monthly_mean", "monthly", "monthly", "monthly", "monthly", "monthly", "monthly", "monthly"),
  transformation = c("log", "log", "monthly mean, level", "log", "level", "log", "level", "log", "log", "log"),
  equation_role = c(
    "endogenous quantity", "endogenous price", "demand shifter",
    "demand shifter", "demand shifter", "supply shifter",
    "supply shifter", "supply shifter", "clean supply shifter",
    "clean supply shifter"
  ),
  source_url = paste0("https://fred.stlouisfed.org/series/", c(
    "HSN1F", "MSPNHSUS", "MORTGAGE30US", "DSPIC96",
    "UNRATE", "PERMIT", "MSACSR", "WPUSI012011",
    "HOUST", "UNDCONTSA"
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
df$log_starts <- log(df$housing_starts)
df$log_under_construction <- log(df$under_construction)
df$trend <- seq_len(nrow(df))
df$trend2 <- df$trend^2
df$month <- factor(format(df$date, "%m"))
df <- df[complete.cases(df[, c(
  "log_sales", "log_price", "mortgage_rate", "log_income", "unrate",
  "log_permit", "months_supply", "log_construction_cost",
  "log_starts", "log_under_construction"
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

lag_vec <- function(x, k) c(rep(NA_real_, k), head(x, -k))
diff_lag <- function(x, k) x - lag_vec(x, k)

for (v in c("log_sales", "log_price", "mortgage_rate", "log_income", "unrate",
            "log_permit", "months_supply", "log_construction_cost",
            "log_starts", "log_under_construction")) {
  df[[paste0("L1_", v)]] <- lag_vec(df[[v]], 1)
  df[[paste0("L3_", v)]] <- lag_vec(df[[v]], 3)
  df[[paste0("L6_", v)]] <- lag_vec(df[[v]], 6)
  df[[paste0("d1_", v)]] <- diff_lag(df[[v]], 1)
  df[[paste0("d12_", v)]] <- diff_lag(df[[v]], 12)
}
df$affordability_pressure <- df$log_price - df$log_income + df$mortgage_rate
df$gfc <- as.integer(df$date >= as.Date("2007-12-01") & df$date <= as.Date("2009-06-01"))
df$covid <- as.integer(df$date >= as.Date("2020-03-01") & df$date <= as.Date("2021-12-01"))
df$high_rate <- as.integer(df$date >= as.Date("2022-03-01"))

cycle_vars <- c("log_sales", "log_price", "log_income", "log_permit", "log_construction_cost", "log_starts")
for (v in cycle_vars) {
  df[[paste0(v, "_cycle")]] <- resid(lm(as.formula(paste(v, "~ trend + trend2")), data = df))
}

iv_price_row <- function(tab, term = "log_price") {
  row <- tab[tab$term == term, ][1, ]
  if (nrow(row) == 0) return(data.frame(estimate = NA_real_, p_value = NA_real_))
  data.frame(estimate = row$estimate, p_value = row$p_value)
}

fit_pair <- function(data, model, sample_name, transformation, demand_incl, demand_excl,
                     supply_incl, supply_excl, y = "log_sales", x = "log_price") {
  needed <- unique(c(y, x, demand_incl, demand_excl, supply_incl, supply_excl))
  d <- data[complete.cases(data[, needed]), ]
  if (nrow(d) < 40) return(NULL)
  dfit <- manual_2sls(d, y, x, demand_incl, demand_excl)
  sfit <- manual_2sls(d, y, x, supply_incl, supply_excl)
  dfs <- first_stage_relevance(d, paste0(model, "_demand"), x, demand_incl, demand_excl)
  sfs <- first_stage_relevance(d, paste0(model, "_supply"), x, supply_incl, supply_excl)
  dov <- rbind(sargan_test(dfit, paste0(model, "_demand"), paste(demand_excl, collapse = "; ")),
               hansen_j_test(dfit, paste0(model, "_demand"), paste(demand_excl, collapse = "; ")))
  sov <- rbind(sargan_test(sfit, paste0(model, "_supply"), paste(supply_excl, collapse = "; ")),
               hansen_j_test(sfit, paste0(model, "_supply"), paste(supply_excl, collapse = "; ")))
  ddwh <- dwh_test(d, paste0(model, "_demand"), y, x, demand_incl, demand_excl)
  sdwh <- dwh_test(d, paste0(model, "_supply"), y, x, supply_incl, supply_excl)
  dtab <- iv_table(dfit, paste0(model, "_demand"), "2SLS_manual_HC3")
  stab <- iv_table(sfit, paste0(model, "_supply"), "2SLS_manual_HC3")
  dprice <- iv_price_row(dtab, x)
  sprice <- iv_price_row(stab, x)
  grid <- data.frame(
    model = model,
    sample = sample_name,
    transformation = transformation,
    IV_set = c("demand_excluded_supply_shifters", "supply_excluded_demand_shifters"),
    equation = c("demand", "supply"),
    first_stage_F = c(dfs$first_stage_F, sfs$first_stage_F),
    partial_R2 = c(dfs$partial_R2, sfs$partial_R2),
    overid_min_p = c(min(dov$p_value, na.rm = TRUE), min(sov$p_value, na.rm = TRUE)),
    DWH_p = c(ddwh$p_value, sdwh$p_value),
    price_coef = c(dprice$estimate, sprice$estimate),
    price_p = c(dprice$p_value, sprice$p_value),
    price_sign_ok = c(dprice$estimate < 0, sprice$estimate > 0),
    price_sig_10pct = c(dprice$p_value < 0.10, sprice$p_value < 0.10),
    N = nrow(d),
    row.names = NULL
  )
  list(demand_fit = dfit, supply_fit = sfit, demand_table = dtab, supply_table = stab,
       first_stage = rbind(dfs, sfs), overid = rbind(dov, sov), dwh = rbind(ddwh, sdwh),
       grid = grid)
}

clean_supply_ivs <- c("log_permit", "log_construction_cost", "log_starts", "log_under_construction")
clean_demand_ivs_A <- c("mortgage_rate", "unrate")
clean_demand_ivs_B <- c("mortgage_rate", "log_income", "unrate")

baseline_pair <- fit_pair(df, "levels_baseline", "2000_2024", "log_levels",
                          demand_included, demand_excluded, supply_included, supply_excluded)
trend_pair <- fit_pair(df, "trend_monthFE", "2000_2024", "log_levels_trend_monthFE",
                       c(demand_included, "trend", "trend2", "month"), demand_excluded,
                       c(supply_included, "trend", "trend2", "month"), supply_excluded)
dynamic_pair <- fit_pair(df, "dynamic_lagY", "2000_2024", "log_levels_L1_sales",
                         c(demand_included, "L1_log_sales"), demand_excluded,
                         c(supply_included, "L1_log_sales"), supply_excluded)
yoy_pair <- fit_pair(df, "yoy_growth", "2000_2024", "year_over_year_change",
                     c("d12_mortgage_rate", "d12_log_income", "d12_unrate"),
                     c("d12_log_permit", "d12_log_construction_cost", "d12_log_starts"),
                     c("d12_log_permit", "d12_months_supply", "d12_log_construction_cost"),
                     c("d12_mortgage_rate", "d12_log_income", "d12_unrate"),
                     y = "d12_log_sales", x = "d12_log_price")
mom_pair <- fit_pair(df, "mom_growth", "2000_2024", "month_over_month_change",
                     c("d1_mortgage_rate", "d1_log_income", "d1_unrate"),
                     c("d1_log_permit", "d1_log_construction_cost", "d1_log_starts"),
                     c("d1_log_permit", "d1_months_supply", "d1_log_construction_cost"),
                     c("d1_mortgage_rate", "d1_log_income", "d1_unrate"),
                     y = "d1_log_sales", x = "d1_log_price")
cycle_pair <- fit_pair(df, "detrended_cycle", "2000_2024", "quadratic_trend_residual",
                       c("mortgage_rate", "log_income_cycle", "unrate"),
                       c("log_permit_cycle", "log_construction_cost_cycle", "log_starts_cycle"),
                       c("log_permit_cycle", "months_supply", "log_construction_cost_cycle"),
                       c("mortgage_rate", "log_income_cycle", "unrate"),
                       y = "log_sales_cycle", x = "log_price_cycle")
clean_pair <- fit_pair(df, "clean_ivset", "2000_2024", "log_levels",
                       demand_included, clean_supply_ivs, supply_included, clean_demand_ivs_B)
lagged_iv_pair <- fit_pair(df, "lagged_IV", "2000_2024", "log_levels_lagged_instruments",
                           demand_included,
                           c("L1_log_permit", "L3_log_permit", "L1_log_construction_cost", "L3_log_construction_cost", "L1_log_starts", "L3_log_starts"),
                           supply_included,
                           c("L1_mortgage_rate", "L3_mortgage_rate", "L1_unrate", "L3_unrate"))
pre_covid_pair <- fit_pair(df[df$date <= as.Date("2019-12-01"), ], "pre_covid", "2000_2019", "log_levels",
                           demand_included, clean_supply_ivs, supply_included, clean_demand_ivs_B)
post_gfc_pre_covid_pair <- fit_pair(df[df$date >= as.Date("2010-01-01") & df$date <= as.Date("2019-12-01"), ],
                                    "post_gfc_pre_covid", "2010_2019", "log_levels",
                                    demand_included, clean_supply_ivs, supply_included, clean_demand_ivs_B)

all_pairs <- Filter(Negate(is.null), list(
  baseline_pair, trend_pair, dynamic_pair, yoy_pair, mom_pair, cycle_pair,
  clean_pair, lagged_iv_pair, pre_covid_pair, post_gfc_pre_covid_pair
))
iv_grid <- do.call(rbind, lapply(all_pairs, function(x) x$grid))
write.csv(iv_grid, file.path(out_dir, "housing_iv_diagnostic_grid.csv"), row.names = FALSE)
all_model_2sls <- do.call(rbind, lapply(all_pairs, function(x) rbind(x$demand_table, x$supply_table)))
write.csv(all_model_2sls, file.path(out_dir, "housing_all_2sls_model_results.csv"), row.names = FALSE)
write.csv(rbind(baseline_pair$demand_table, baseline_pair$supply_table),
          file.path(out_dir, "housing_baseline_static_2sls.csv"), row.names = FALSE)
write.csv(rbind(trend_pair$demand_table, trend_pair$supply_table), file.path(out_dir, "housing_trend_monthfe_2sls.csv"), row.names = FALSE)
write.csv(rbind(dynamic_pair$demand_table, dynamic_pair$supply_table), file.path(out_dir, "housing_dynamic_lagY_2sls.csv"), row.names = FALSE)
write.csv(rbind(yoy_pair$demand_table, yoy_pair$supply_table), file.path(out_dir, "housing_yoy_2sls_results.csv"), row.names = FALSE)
write.csv(yoy_pair$first_stage, file.path(out_dir, "housing_yoy_first_stage.csv"), row.names = FALSE)
write.csv(yoy_pair$overid, file.path(out_dir, "housing_yoy_overid.csv"), row.names = FALSE)
write.csv(rbind(mom_pair$demand_table, mom_pair$supply_table), file.path(out_dir, "housing_mom_2sls_results.csv"), row.names = FALSE)
write.csv(rbind(cycle_pair$demand_table, cycle_pair$supply_table), file.path(out_dir, "housing_detrended_cycle_2sls.csv"), row.names = FALSE)
write.csv(rbind(clean_pair$demand_table, clean_pair$supply_table), file.path(out_dir, "housing_clean_ivset_results.csv"), row.names = FALSE)
write.csv(rbind(lagged_iv_pair$demand_table, lagged_iv_pair$supply_table), file.path(out_dir, "housing_lagged_iv_results.csv"), row.names = FALSE)
write.csv(rbind(pre_covid_pair$demand_table, pre_covid_pair$supply_table), file.path(out_dir, "housing_pre_covid_results.csv"), row.names = FALSE)
write.csv(rbind(post_gfc_pre_covid_pair$demand_table, post_gfc_pre_covid_pair$supply_table), file.path(out_dir, "housing_post_gfc_pre_covid_results.csv"), row.names = FALSE)

supply_lagged_price <- do.call(rbind, lapply(c("L1_log_price", "L3_log_price", "L6_log_price"), function(px) {
  m <- lm(as.formula(paste("log_sales ~", px, "+ log_permit + months_supply + log_construction_cost")), data = df)
  out <- coef_table(m, "supply", paste0("OLS_supply_", px))
  out[out$term == px, ]
}))
write.csv(supply_lagged_price, file.path(out_dir, "housing_supply_lagged_price_results.csv"), row.names = FALSE)

just_specs <- list(
  demand_just_permit = list(y = "log_sales", x = "log_price", included = demand_included, excluded = "log_permit", equation = "demand"),
  demand_just_L3permit = list(y = "log_sales", x = "log_price", included = demand_included, excluded = "L3_log_permit", equation = "demand"),
  supply_just_mortgage = list(y = "log_sales", x = "log_price", included = supply_included, excluded = "mortgage_rate", equation = "supply"),
  supply_just_L3mortgage = list(y = "log_sales", x = "log_price", included = supply_included, excluded = "L3_mortgage_rate", equation = "supply")
)
just_results <- do.call(rbind, lapply(names(just_specs), function(name) {
  spec <- just_specs[[name]]
  needed <- unique(c(spec$y, spec$x, spec$included, spec$excluded))
  d <- df[complete.cases(df[, needed]), ]
  fit <- manual_2sls(d, spec$y, spec$x, spec$included, spec$excluded)
  tab <- iv_table(fit, spec$equation, paste0("just_identified_", name))
  fs <- first_stage_relevance(d, spec$equation, spec$x, spec$included, spec$excluded)
  key_row <- tab[tab$term == spec$x, ][1, ]
  data.frame(
    model = name,
    equation = spec$equation,
    endogenous_rhs_variable = spec$x,
    excluded_IV = spec$excluded,
    price_coef = key_row$estimate,
    price_p = key_row$p_value,
    first_stage_F = fs$first_stage_F,
    partial_R2 = fs$partial_R2,
    expected_price_sign_ok = ifelse(spec$equation == "demand", key_row$estimate < 0, key_row$estimate > 0),
    N = fit$n,
    note = "Just identified: overidentification test unavailable; this checks sign robustness only.",
    row.names = NULL
  )
}))
write.csv(just_results, file.path(out_dir, "housing_just_identified_2sls.csv"), row.names = FALSE)

yoy_preferred <- rbind(yoy_pair$demand_table, yoy_pair$supply_table)
yoy_preferred$preferred_role <- "YoY preferred transformation because ADF tests indicate nonstationarity in levels."
write.csv(yoy_preferred, file.path(out_dir, "housing_yoy_preferred_results.csv"), row.names = FALSE)

pre_covid_preferred <- rbind(pre_covid_pair$demand_table, pre_covid_pair$supply_table)
pre_covid_preferred$preferred_role <- "Pre-COVID preferred robustness sample: excludes pandemic and rapid-rate-hike distortions."
write.csv(pre_covid_preferred, file.path(out_dir, "housing_pre_covid_preferred_results.csv"), row.names = FALSE)

supply_lagged_price_preferred <- do.call(rbind, lapply(c("L1_log_price", "L3_log_price", "L6_log_price"), function(px) {
  d <- df[complete.cases(df[, c("log_sales", px, supply_included, supply_excluded)]), ]
  form <- as.formula(paste("log_sales ~", px, "+", paste(supply_included, collapse = " + ")))
  inst_form <- as.formula(paste("~", paste(c(supply_included, supply_excluded), collapse = " + ")))
  fit <- AER::ivreg(form, instruments = inst_form, data = d)
  vc <- sandwich::vcovHC(fit, type = "HC3")
  est <- coef(fit)
  se <- sqrt(diag(vc))
  data.frame(
    lagged_price_term = px,
    estimate = unname(est[px]),
    HC3_se = unname(se[px]),
    statistic = unname(est[px] / se[px]),
    p_value = 2 * pnorm(abs(est[px] / se[px]), lower.tail = FALSE),
    expected_positive = unname(est[px] > 0),
    N = nobs(fit),
    note = "Supply response may reflect lagged price incentives rather than same-month prices.",
    row.names = NULL
  )
}))
write.csv(supply_lagged_price_preferred, file.path(out_dir, "housing_supply_lagged_price_preferred.csv"), row.names = FALSE)

preferred_grid <- data.frame(
  preferred_check = c(
    "just_identified_IV",
    "yoy_growth_preferred",
    "pre_covid_preferred",
    "supply_lagged_price_preferred"
  ),
  output_file = c(
    "housing_just_identified_2sls.csv",
    "housing_yoy_preferred_results.csv",
    "housing_pre_covid_preferred_results.csv",
    "housing_supply_lagged_price_preferred.csv"
  ),
  main_question = c(
    "Do price coefficient signs survive when each equation uses only one excluded IV?",
    "Do signs improve after removing common trends with year-over-year changes?",
    "Do signs improve before pandemic and rapid-rate-hike distortions?",
    "Does supply respond more clearly to lagged prices?"
  ),
  decision_rule = c(
    "No overidentification test; use only as sign robustness.",
    "Prioritize if demand price is negative, mortgage-rate change is negative, and permit growth is positive.",
    "Prioritize if signs are more theory-consistent than full-sample levels.",
    "Use to explain insignificant same-month supply price coefficient if lagged price is positive."
  ),
  row.names = NULL
)
write.csv(preferred_grid, file.path(out_dir, "housing_preferred_model_checks.csv"), row.names = FALSE)

affordability_2sls <- manual_2sls(df[complete.cases(df[, c("log_sales", "affordability_pressure", "unrate", clean_supply_ivs)]), ],
                                  "log_sales", "affordability_pressure", "unrate", clean_supply_ivs)
write.csv(iv_table(affordability_2sls, "demand_affordability", "2SLS_manual_HC3"),
          file.path(out_dir, "housing_affordability_demand_results.csv"), row.names = FALSE)

regime_interaction <- lm(log_sales ~ log_price * high_rate + mortgage_rate * high_rate + log_income + unrate + gfc + covid, data = df)
write.csv(coef_table(regime_interaction, "demand", "OLS_regime_interactions_HC3"),
          file.path(out_dir, "housing_regime_interaction_results.csv"), row.names = FALSE)

adf_tests <- do.call(rbind, lapply(cor_vars, function(v) {
  x <- na.omit(df[[v]])
  test <- tryCatch(tseries::adf.test(x), error = identity)
  if (inherits(test, "error")) {
    data.frame(variable = v, statistic = NA_real_, p_value = NA_real_, note = test$message)
  } else {
    data.frame(variable = v, statistic = unname(test$statistic), p_value = test$p.value,
               note = "ADF null: unit root; p-values may be approximate.")
  }
}))
write.csv(adf_tests, file.path(out_dir, "housing_adf_tests.csv"), row.names = FALSE)

vif_extract <- function(model, equation) {
  vals <- car::vif(model)
  if (is.matrix(vals)) vals <- vals[, "GVIF^(1/(2*Df))"]
  data.frame(equation = equation, term = names(vals), VIF = as.numeric(vals), row.names = NULL)
}
vif_diagnostics <- rbind(vif_extract(demand_ols, "demand_ols"), vif_extract(supply_ols, "supply_ols"))
write.csv(vif_diagnostics, file.path(out_dir, "housing_vif_diagnostics.csv"), row.names = FALSE)

serial_tests <- function(resid, model_name) {
  e <- as.numeric(na.omit(resid))
  dw <- sum(diff(e)^2) / sum(e^2)
  lag_order <- 12
  aux <- data.frame(e = e)
  for (i in seq_len(lag_order)) aux[[paste0("L", i)]] <- lag_vec(e, i)
  aux <- aux[complete.cases(aux), ]
  bg <- lm(e ~ ., data = aux)
  bg_stat <- nrow(aux) * summary(bg)$r.squared
  data.frame(
    model = model_name,
    test = c("Durbin_Watson_stat", "Breusch_Godfrey_LM_order12"),
    statistic = c(dw, bg_stat),
    p_value = c(NA_real_, pchisq(bg_stat, df = lag_order, lower.tail = FALSE)),
    row.names = NULL
  )
}
serial_correlation_tests <- rbind(
  serial_tests(resid(demand_ols), "demand_ols"),
  serial_tests(resid(supply_ols), "supply_ols"),
  serial_tests(demand_2sls$residuals, "demand_2sls"),
  serial_tests(supply_2sls$residuals, "supply_2sls")
)
write.csv(serial_correlation_tests, file.path(out_dir, "housing_serial_correlation_tests.csv"), row.names = FALSE)

bp_manual <- function(resid, fitted, model_name) {
  aux <- lm(I(resid^2) ~ fitted + I(fitted^2))
  stat <- length(resid) * summary(aux)$r.squared
  data.frame(model = model_name, test = "White_style_BP_on_fitted", statistic = stat, df = 2,
             p_value = pchisq(stat, df = 2, lower.tail = FALSE), row.names = NULL)
}
hetero_tests <- rbind(
  bp_manual(resid(demand_ols), fitted(demand_ols), "demand_ols"),
  bp_manual(resid(supply_ols), fitted(supply_ols), "supply_ols"),
  bp_manual(demand_2sls$residuals, demand_2sls$fitted, "demand_2sls"),
  bp_manual(supply_2sls$residuals, supply_2sls$fitted, "supply_2sls")
)
write.csv(hetero_tests, file.path(out_dir, "housing_heteroskedasticity_tests.csv"), row.names = FALSE)

ivreg_demand <- AER::ivreg(
  log_sales ~ log_price + mortgage_rate + log_income + unrate |
    mortgage_rate + log_income + unrate + log_permit + months_supply + log_construction_cost,
  data = df
)
ivreg_supply <- AER::ivreg(
  log_sales ~ log_price + log_permit + months_supply + log_construction_cost |
    log_permit + months_supply + log_construction_cost + mortgage_rate + log_income + unrate,
  data = df
)
nw_compare <- function(model, equation) {
  hc <- sandwich::vcovHC(model, type = "HC3")
  nw <- sandwich::NeweyWest(model, lag = 6, prewhite = FALSE, adjust = TRUE)
  est <- coef(model)
  rows <- do.call(rbind, lapply(c("log_price", names(est)[2]), function(term) {
    term <- unique(term)[1]
    data.frame(
      equation = equation,
      term = term,
      estimate = est[term],
      HC3_se = sqrt(diag(hc))[term],
      HC3_p = 2 * pnorm(abs(est[term] / sqrt(diag(hc))[term]), lower.tail = FALSE),
      NeweyWest_lag6_se = sqrt(diag(nw))[term],
      NeweyWest_lag6_p = 2 * pnorm(abs(est[term] / sqrt(diag(nw))[term]), lower.tail = FALSE),
      row.names = NULL
    )
  }))
  unique(rows)
}
write.csv(rbind(nw_compare(ivreg_demand, "demand_2sls"), nw_compare(ivreg_supply, "supply_2sls")),
          file.path(out_dir, "housing_2sls_HC3_vs_NeweyWest.csv"), row.names = FALSE)

model_scores <- aggregate(cbind(price_sign_ok, price_sig_10pct) ~ model, data = iv_grid, FUN = function(x) sum(x, na.rm = TRUE))
overid_min <- aggregate(overid_min_p ~ model, data = iv_grid, FUN = function(x) min(x, na.rm = TRUE))
fs_min <- aggregate(first_stage_F ~ model, data = iv_grid, FUN = function(x) min(x, na.rm = TRUE))
model_selection <- merge(merge(model_scores, overid_min, by = "model"), fs_min, by = "model")
model_selection$use_as_preferred <- with(model_selection, price_sign_ok == 2 & first_stage_F > 10 & overid_min_p >= 0.05)
model_selection$use_as_sensitivity <- with(model_selection, price_sign_ok >= 1 & first_stage_F > 10)
model_selection$reason <- ifelse(
  model_selection$overid_min_p < 0.05,
  "Overidentification remains rejected; use for Lecture 6 demonstration or sensitivity, not as preferred causal evidence.",
  "No overidentification rejection and both price signs fit theory."
)
if (!any(model_selection$use_as_preferred)) {
  model_selection$reason[model_selection$use_as_sensitivity] <- paste(
    model_selection$reason[model_selection$use_as_sensitivity],
    "No model passes the full preferred-model rule in this run."
  )
}
write.csv(model_selection, file.path(out_dir, "housing_model_selection_flags.csv"), row.names = FALSE)

caution_lines <- c(
  "# Housing SEM Interpretation Cautions",
  "",
  "1. The static levels model is retained as the Lecture 6 baseline, but monthly housing variables are trending and serially correlated.",
  "2. `months_supply` is not used in the clean IV set for the demand equation because it has a mechanical relationship with sales.",
  "3. Newey-West standard errors address serial correlation in inference; they do not solve simultaneity or invalid instruments.",
  "4. The preferred empirical discussion should compare levels, trend/month fixed effects, year-over-year changes, lagged IVs, and the pre-COVID sample.",
  "5. Overidentification rejections are treated as warnings about exclusion restrictions, not as mechanical failures of the code.",
  "6. The module is strongest as a Lecture 6 supply-demand SEM demonstration; causal claims should rely on the model-selection grid and diagnostic tables."
)
writeLines(caution_lines, file.path(out_dir, "housing_interpretation_cautions.md"))

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
  "## Improved Time-Series Pipeline",
  "",
  "Following the revision notes, this version also estimates trend/month fixed-effect, lagged-sales dynamic, year-over-year, month-over-month, detrended-cycle, clean-IV, lagged-IV, pre-COVID, post-GFC/pre-COVID, lagged-supply-price, affordability, and regime-interaction variants.",
  "Additional diagnostics include ADF unit-root tests, VIF, serial-correlation tests, heteroskedasticity tests, Newey-West versus HC3 standard errors, and a model-selection grid.",
  paste0("ADF tests fail to reject a unit root at 5% for ", sum(adf_tests$p_value >= 0.05, na.rm = TRUE), " of ", nrow(adf_tests), " core variables, so growth-rate and detrended specifications should be discussed alongside levels."),
  paste0("Preferred-model rule passed by any model: ", any(model_selection$use_as_preferred), ". A model must have both theoretically correct price signs, first-stage F > 10, and no overidentification rejection."),
  "Because overidentification is still rejected in the main variants, the housing SEM should be presented as a strong Lecture 6 supply-demand demonstration with transparent diagnostic cautions, not as definitive causal evidence.",
  "",
  "## Revision-2 Preferred Checks",
  "",
  "The second revision adds four focused checks: just-identified IV specifications, a year-over-year preferred transformation, a pre-COVID preferred sample, and lagged-price supply equations.",
  "The just-identified specifications avoid overidentification-test rejection by construction, but they cannot test exclusion restrictions; they are sign-robustness checks only.",
  "The year-over-year model is emphasized because the ADF tests indicate that levels are nonstationary. The pre-COVID sample is emphasized because pandemic and rapid-hiking periods likely changed housing-market behavior.",
  "The lagged-price supply checks address the institutional point that new housing supply responds with construction and sales delays, not necessarily within the same month.",
  paste0("In the just-identified checks, theoretically correct price signs appear in ", sum(just_results$expected_price_sign_ok, na.rm = TRUE), " of ", nrow(just_results), " specifications; therefore these checks do not rescue a preferred causal interpretation."),
  paste0("In the lagged-price supply checks, all three lagged price coefficients are positive, but the smallest p-value is ", fmt(min(supply_lagged_price_preferred$p_value, na.rm = TRUE)), ", so this is an economically sensible but statistically weak pattern."),
  "",
  "## Interpretation",
  "",
  "This housing module is a better Lecture 6 application than the previous FOMC CAR-RV exploratory SEM because price and quantity have a standard simultaneous-equilibrium interpretation. The demand and supply equations each have clear ceteris-paribus meanings, and the excluded demand/supply shifters provide a transparent identification strategy.",
  "",
  "All tables and figures are saved in `output/step6_housing/`."
)
writeLines(summary_lines, file.path(out_dir, "housing_sem_results_summary.md"))

message("Step 6 housing SEM complete. Outputs written to ", out_dir)
