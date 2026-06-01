library(OPTICS)

if (!exists("OPTICS.ms", mode = "function")) {
  stop("OPTICS.ms is not available in the loaded OPTICS package.")
}

simulate_mean_1d <- function(n_total, cps, a, dist = c("normal", "t10")) {
  dist <- match.arg(dist)
  cps_all <- c(0L, as.integer(cps), as.integer(n_total))
  mu <- numeric(n_total)
  for (seg in seq_len(length(cps_all) - 1L)) {
    idx <- (cps_all[seg] + 1L):cps_all[seg + 1L]
    mu[idx] <- if (seg %% 2L == 1L) a else -a
  }
  eps <- if (dist == "normal") stats::rnorm(n_total) else stats::rt(n_total, df = 10)
  matrix(mu + eps, nrow = 1L)
}

fit_with_diag <- function(expr) {
  out <- tryCatch(expr, error = function(e) e)
  if (inherits(out, "error")) {
    list(ok = FALSE, fit = NULL, msg = out$message)
  } else {
    list(ok = TRUE, fit = out, msg = "")
  }
}

run_one_rep_diag <- function(i, n_total, cps, a, dist, alpha, B, L, base_seed = 1234L) {
  true_k <- length(cps)
  k_max <- max(2L, ceiling(log(n_total / 2)))
  model_bs <- seq_len(k_max)
  model_sn <- model_bs + 1

  set.seed(base_seed + 100000L * as.integer(n_total) + 1000L * as.integer(100 * a) + 100L * i + if (dist == "t10") 1L else 0L)
  X <- simulate_mean_1d(n_total = n_total, cps = cps, a = a, dist = dist)

  base_bs <- fit_with_diag(OPTICS.mean(X, model = model_bs, method = "BinSeg", alpha = alpha, B = B))
  ms_bs <- fit_with_diag(OPTICS.ms(X, model = model_bs, engine = "mean", method = "BinSeg", alpha = alpha, B = B, L = L, trim_incomplete = TRUE))
  base_sn <- fit_with_diag(OPTICS.mean(X, model = model_sn, method = "SegNeigh", alpha = alpha, B = B))
  ms_sn <- fit_with_diag(OPTICS.ms(X, model = model_sn, engine = "mean", method = "SegNeigh", alpha = alpha, B = B, L = L, trim_incomplete = TRUE))

  data.frame(
    rep = i,
    n_total = n_total,
    a = a,
    dist = dist,
    base_bs_ok = base_bs$ok,
    base_bs_cov = if (base_bs$ok) as.numeric(true_k %in% base_bs$fit$A) else NA_real_,
    base_bs_len = if (base_bs$ok) length(base_bs$fit$A) else NA_real_,
    base_bs_err = base_bs$msg,
    ms_bs_ok = ms_bs$ok,
    ms_bs_cov = if (ms_bs$ok) as.numeric(true_k %in% ms_bs$fit$A_MS) else NA_real_,
    ms_bs_len = if (ms_bs$ok) length(ms_bs$fit$A_MS) else NA_real_,
    ms_bs_err = ms_bs$msg,
    base_sn_ok = base_sn$ok,
    base_sn_cov = if (base_sn$ok) as.numeric((true_k + 1) %in% base_sn$fit$A) else NA_real_,
    base_sn_len = if (base_sn$ok) length(base_sn$fit$A) else NA_real_,
    base_sn_err = base_sn$msg,
    ms_sn_ok = ms_sn$ok,
    ms_sn_cov = if (ms_sn$ok) as.numeric((true_k + 1) %in% ms_sn$fit$A_MS) else NA_real_,
    ms_sn_len = if (ms_sn$ok) length(ms_sn$fit$A_MS) else NA_real_,
    ms_sn_err = ms_sn$msg,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Smoke check configuration
# ---------------------------------------------------------------------------
n_total <- 1000
cps <- c(200, 400, 600, 800)
a_grid <- seq(0.5, 1, length.out = 5)
dist_grid <- c("normal", "t10")
reps <- 5
alpha <- 0.1
B <- 50
L <- 3

dir.create("./results", recursive = TRUE, showWarnings = FALSE)

message(sprintf("OPTICS version: %s", as.character(utils::packageVersion("OPTICS"))))
message(sprintf("has OPTICS.ms: %s", exists("OPTICS.ms", mode = "function")))
message(sprintf("n_total=%d, reps=%d, alpha=%.3f, B=%d, L=%d", n_total, reps, alpha, B, L))

detail <- list()
idx <- 1L
for (dist in dist_grid) {
  for (a in a_grid) {
    message(sprintf("[%s] dist=%s, a=%.3f", Sys.time(), dist, a))
    for (i in seq_len(reps)) {
      detail[[idx]] <- run_one_rep_diag(
        i = i, n_total = n_total, cps = cps, a = a, dist = dist,
        alpha = alpha, B = B, L = L
      )
      idx <- idx + 1L
    }
  }
}

detail_df <- do.call(rbind, detail)

summary_df <- stats::aggregate(
  cbind(
    base_bs_ok = as.numeric(base_bs_ok),
    ms_bs_ok = as.numeric(ms_bs_ok),
    base_sn_ok = as.numeric(base_sn_ok),
    ms_sn_ok = as.numeric(ms_sn_ok),
    base_bs_cov, ms_bs_cov, base_sn_cov, ms_sn_cov,
    base_bs_len, ms_bs_len, base_sn_len, ms_sn_len
  ) ~ dist + a,
  data = detail_df,
  FUN = function(x) mean(x, na.rm = TRUE)
)

write.csv(detail_df, "./results/check_ms_detail.csv", row.names = FALSE)
write.csv(summary_df, "./results/check_ms_summary.csv", row.names = FALSE)

print(summary_df)

# print first few non-empty error messages for debugging
err_cols <- c("base_bs_err", "ms_bs_err", "base_sn_err", "ms_sn_err")
for (col in err_cols) {
  msgs <- unique(detail_df[[col]])
  msgs <- msgs[nzchar(msgs)]
  if (length(msgs) > 0) {
    cat("\n", col, " unique errors:\n", sep = "")
    print(utils::head(msgs, 5))
  }
}
