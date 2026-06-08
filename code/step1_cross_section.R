options(stringsAsFactors = FALSE)
options(timeout = 600)

`%||%` <- function(left, right) {
  if (is.null(left)) right else left
}

this_file <- tryCatch(normalizePath(sys.frames()[[1]]$ofile), error = function(e) NA_character_)
if (is.na(this_file)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) this_file <- normalizePath(sub("^--file=", "", file_arg[1]))
}
project_root <- normalizePath(file.path(dirname(this_file), ".."))
raw_dir <- file.path(project_root, "data", "raw")
out_dir <- file.path(project_root, "output", "step1")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

event_date <- as.Date("2022-11-02")
estimation_start_lag <- 260
estimation_end_lag <- 30

ff_url <- "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/49_Industry_Portfolios_daily_CSV.zip"
fred_url <- "https://fred.stlouisfed.org/graph/fredgraph.csv?id=DGS2"

ff_zip <- file.path(raw_dir, "49_Industry_Portfolios_daily_CSV.zip")
if (!file.exists(ff_zip)) {
  ff_zip_alt <- file.path(raw_dir, "49_Industry_Portfolios_Daily_CSV.zip")
  if (file.exists(ff_zip_alt)) ff_zip <- ff_zip_alt
}
dgs2_csv <- file.path(raw_dir, "DGS2.csv")
frb_h15_csv <- file.path(raw_dir, "FRB_H15.csv")

download_if_missing <- function(url, dest) {
  needs_download <- !file.exists(dest) || file.info(dest)$size == 0
  if (!needs_download && grepl("\\.zip$", dest, ignore.case = TRUE)) {
    needs_download <- inherits(try(unzip(dest, list = TRUE), silent = TRUE), "try-error")
  }
  if (needs_download) {
    message("Downloading: ", url)
    temp_dest <- paste0(dest, ".download")
    if (file.exists(temp_dest)) unlink(temp_dest)
    download.file(url, temp_dest, mode = "wb", quiet = FALSE, method = "libcurl")
    file.rename(temp_dest, dest)
  }
}

parse_ff49_daily <- function(zip_path) {
  temp_dir <- tempfile("ff49_")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
  files <- unzip(zip_path, exdir = temp_dir)
  csv_file <- files[grepl("\\.csv$", files, ignore.case = TRUE)][1]
  lines <- readLines(csv_file, warn = FALSE)

  header_line <- grep("^\\s*,", lines)[1]
  if (is.na(header_line)) stop("Cannot locate Fama-French 49 industry value-weighted table header.")

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
  if (!("RIFLGFCY02_N.B" %in% names(dat))) stop("Cannot locate 2-year Treasury column RIFLGFCY02_N.B in FRB H15 file.")
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
  if (file.exists(frb_h15_csv) && file.info(frb_h15_csv)$size > 0) {
    message("Using local FRB H15 file: ", frb_h15_csv)
    return(parse_frb_h15_dgs2(frb_h15_csv))
  }
  if (file.exists(dgs2_csv) && file.info(dgs2_csv)$size > 0) {
    message("Using local FRED DGS2 file: ", dgs2_csv)
    return(parse_fred_dgs2(dgs2_csv))
  }
  download_if_missing(fred_url, dgs2_csv)
  parse_fred_dgs2(dgs2_csv)
}

write_model_summary <- function(model, file) {
  sink(file)
  cat("Step 1: Simple cross-sectional OLS\n")
  cat("Event date:", format(event_date), "\n")
  cat("Event window: [-1,+1] trading days around event date\n")
  cat("Estimation window:", estimation_start_lag, "to", estimation_end_lag, "trading days before event\n\n")
  print(summary(model))
  sink()
}

hc1_table <- function(model) {
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
    variable = names(coef(model)),
    coefficient = unname(coef(model)),
    robust_std_error = unname(robust_se),
    t_value = unname(t_value),
    p_value = unname(p_value),
    row.names = NULL
  )
}

download_if_missing(ff_url, ff_zip)

ff49 <- parse_ff49_daily(ff_zip)
dgs2 <- read_rate_data()

event_index <- match(event_date, ff49$date)
if (is.na(event_index)) stop("Event date is not in Fama-French daily return data: ", event_date)
if (event_index <= estimation_start_lag) stop("Not enough pre-event data for estimation window.")

event_window_index <- (event_index - 1):(event_index + 1)
estimation_index <- (event_index - estimation_start_lag):(event_index - estimation_end_lag)

industry_cols <- setdiff(names(ff49), "date")
event_returns <- ff49[event_window_index, c("date", industry_cols)]
estimation_returns <- ff49[estimation_index, c("date", industry_cols)]

estimation_data <- merge(estimation_returns, dgs2[, c("date", "d_DGS2")], by = "date", all.x = TRUE)
estimation_data <- estimation_data[!is.na(estimation_data$d_DGS2), ]

results <- data.frame(
  industry = industry_cols,
  CAR = NA_real_,
  rate_beta = NA_real_,
  RateVulnerability = NA_real_,
  n_estimation = NA_integer_
)

for (industry in industry_cols) {
  car <- sum(event_returns[[industry]], na.rm = TRUE)
  reg_data <- estimation_data[, c(industry, "d_DGS2")]
  names(reg_data) <- c("ret", "d_DGS2")
  reg_data <- reg_data[complete.cases(reg_data), ]
  fit <- lm(ret ~ d_DGS2, data = reg_data)
  beta <- unname(coef(fit)["d_DGS2"])

  row <- match(industry, results$industry)
  results$CAR[row] <- car
  results$rate_beta[row] <- beta
  results$RateVulnerability[row] <- -beta
  results$n_estimation[row] <- nrow(reg_data)
}

results$RateVulnerability_z <- as.numeric(scale(results$RateVulnerability))

model_raw <- lm(CAR ~ RateVulnerability, data = results)
model_z <- lm(CAR ~ RateVulnerability_z, data = results)

desc <- data.frame(
  variable = c("CAR", "RateVulnerability", "RateVulnerability_z"),
  n = c(sum(!is.na(results$CAR)), sum(!is.na(results$RateVulnerability)), sum(!is.na(results$RateVulnerability_z))),
  mean = c(mean(results$CAR), mean(results$RateVulnerability), mean(results$RateVulnerability_z)),
  sd = c(sd(results$CAR), sd(results$RateVulnerability), sd(results$RateVulnerability_z)),
  min = c(min(results$CAR), min(results$RateVulnerability), min(results$RateVulnerability_z)),
  max = c(max(results$CAR), max(results$RateVulnerability), max(results$RateVulnerability_z))
)

coef_z <- summary(model_z)$coefficients
reg_table <- data.frame(
  variable = rownames(coef_z),
  coefficient = coef_z[, "Estimate"],
  std_error = coef_z[, "Std. Error"],
  t_value = coef_z[, "t value"],
  p_value = coef_z[, "Pr(>|t|)"],
  row.names = NULL
)
model_stats <- data.frame(
  N = nobs(model_z),
  R_squared = summary(model_z)$r.squared,
  Adj_R_squared = summary(model_z)$adj.r.squared,
  Residual_SE = summary(model_z)$sigma,
  F_statistic = unname(summary(model_z)$fstatistic["value"]),
  F_p_value = pf(
    summary(model_z)$fstatistic["value"],
    summary(model_z)$fstatistic["numdf"],
    summary(model_z)$fstatistic["dendf"],
    lower.tail = FALSE
  )
)

robust_table <- hc1_table(model_z)

loo_results <- do.call(
  rbind,
  lapply(results$industry, function(industry) {
    temp_data <- results[results$industry != industry, ]
    temp_model <- lm(CAR ~ RateVulnerability_z, data = temp_data)
    temp_sum <- summary(temp_model)
    data.frame(
      dropped_industry = industry,
      coefficient = unname(coef(temp_model)["RateVulnerability_z"]),
      p_value = temp_sum$coefficients["RateVulnerability_z", "Pr(>|t|)"],
      R_squared = temp_sum$r.squared,
      row.names = NULL
    )
  })
)

write.csv(results, file.path(out_dir, "step1_cross_section_data.csv"), row.names = FALSE)
write.csv(desc, file.path(out_dir, "step1_descriptive_statistics.csv"), row.names = FALSE)
write.csv(reg_table, file.path(out_dir, "step1_ols_table_standardized.csv"), row.names = FALSE)
write.csv(robust_table, file.path(out_dir, "step1_ols_table_standardized_HC1.csv"), row.names = FALSE)
write.csv(loo_results, file.path(out_dir, "step1_leave_one_out.csv"), row.names = FALSE)
write.csv(model_stats, file.path(out_dir, "step1_model_stats_standardized.csv"), row.names = FALSE)
write_model_summary(model_z, file.path(out_dir, "step1_ols_summary_standardized.txt"))
write_model_summary(model_raw, file.path(out_dir, "step1_ols_summary_raw.txt"))

pdf(file.path(out_dir, "step1_scatter_ols.pdf"), width = 7, height = 5)
plot(
  results$RateVulnerability_z,
  results$CAR,
  pch = 19,
  col = "#1f77b4",
  xlab = "Rate vulnerability (standardized)",
  ylab = "CAR [-1,+1], percent",
  main = "Industry rate vulnerability and FOMC-window returns"
)
abline(model_z, col = "#d62728", lwd = 2)
grid(col = "gray85")
text(
  results$RateVulnerability_z,
  results$CAR,
  labels = results$industry,
  pos = 3,
  cex = 0.55,
  col = "gray30"
)
legend(
  "topright",
  legend = sprintf("OLS slope = %.3f, R^2 = %.3f", coef(model_z)[2], summary(model_z)$r.squared),
  bty = "n"
)
dev.off()

cat("Done.\n")
cat("Event date:", format(event_date), "\n")
cat("Event window dates:", paste(format(event_returns$date), collapse = ", "), "\n")
cat("Estimation rows after DGS2 merge:", nrow(estimation_data), "\n")
cat("Output directory:", out_dir, "\n\n")
cat("Standardized OLS table:\n")
print(reg_table, row.names = FALSE)
cat("\nModel stats:\n")
print(model_stats, row.names = FALSE)
cat("\nHC1 robust OLS table:\n")
print(robust_table, row.names = FALSE)
cat("\nLeave-one-out coefficient range:\n")
print(
  data.frame(
    min_coefficient = min(loo_results$coefficient),
    max_coefficient = max(loo_results$coefficient),
    max_p_value = max(loo_results$p_value),
    significant_at_5pct = sum(loo_results$p_value < 0.05),
    total = nrow(loo_results)
  ),
  row.names = FALSE
)
