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

timed_fit_with_diag <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  out <- fit_with_diag(expr)
  elapsed <- proc.time()[["elapsed"]] - t0
  list(ok = out$ok, fit = out$fit, msg = out$msg, elapsed = as.numeric(elapsed))
}

run_one_rep <- function(i, n_total, cps, a, dist, alpha, B, L, base_seed = 1234L) {
  true_k <- length(cps)
  k_max <- max(2L, ceiling(log(n_total / 2)))
  model_bs <- seq_len(k_max)
  model_sn <- model_bs + 1

  set.seed(base_seed + 100000L * as.integer(n_total) + 1000L * as.integer(100 * a) + 100L * i + if (dist == "t10") 1L else 0L)
  X <- simulate_mean_1d(n_total = n_total, cps = cps, a = a, dist = dist)

  base_bs <- timed_fit_with_diag(OPTICS.mean(X, model = model_bs, method = "BinSeg", alpha = alpha, B = B))
  ms_bs <- timed_fit_with_diag(OPTICS.ms(X, model = model_bs, engine = "mean", method = "BinSeg", alpha = alpha, B = B, L = L, trim_incomplete = TRUE))
  base_sn <- timed_fit_with_diag(OPTICS.mean(X, model = model_sn, method = "SegNeigh", alpha = alpha, B = B))
  ms_sn <- timed_fit_with_diag(OPTICS.ms(X, model = model_sn, engine = "mean", method = "SegNeigh", alpha = alpha, B = B, L = L, trim_incomplete = TRUE))

  data.frame(
    rep = i,
    n_total = n_total,
    a = a,
    dist = dist,
    base_bs_ok = base_bs$ok,
    base_bs_cov = if (base_bs$ok) as.numeric(true_k %in% base_bs$fit$A) else NA_real_,
    base_bs_len = if (base_bs$ok) length(base_bs$fit$A) else NA_real_,
    base_bs_time = base_bs$elapsed,
    ms_bs_ok = ms_bs$ok,
    ms_bs_cov = if (ms_bs$ok) as.numeric(true_k %in% ms_bs$fit$A_MS) else NA_real_,
    ms_bs_len = if (ms_bs$ok) length(ms_bs$fit$A_MS) else NA_real_,
    ms_bs_time = ms_bs$elapsed,
    base_sn_ok = base_sn$ok,
    base_sn_cov = if (base_sn$ok) as.numeric((true_k + 1) %in% base_sn$fit$A) else NA_real_,
    base_sn_len = if (base_sn$ok) length(base_sn$fit$A) else NA_real_,
    base_sn_time = base_sn$elapsed,
    ms_sn_ok = ms_sn$ok,
    ms_sn_cov = if (ms_sn$ok) as.numeric((true_k + 1) %in% ms_sn$fit$A_MS) else NA_real_,
    ms_sn_len = if (ms_sn$ok) length(ms_sn$fit$A_MS) else NA_real_,
    ms_sn_time = ms_sn$elapsed,
    stringsAsFactors = FALSE
  )
}

run_setting <- function(n_total, a, dist, reps, alpha, B, L) {
  cps <- round(n_total * c(0.2, 0.4, 0.6, 0.8))
  cps <- sort(unique(pmax(2L, pmin(n_total - 2L, as.integer(cps)))))

  out <- vector("list", reps)
  for (i in seq_len(reps)) {
    out[[i]] <- run_one_rep(
      i = i, n_total = n_total, cps = cps, a = a, dist = dist,
      alpha = alpha, B = B, L = L
    )
  }
  do.call(rbind, out)
}

summarize_results <- function(df) {
  stats::aggregate(
    cbind(
      base_bs_cov, ms_bs_cov, base_sn_cov, ms_sn_cov,
      base_bs_len, ms_bs_len, base_sn_len, ms_sn_len,
      base_bs_time, ms_bs_time, base_sn_time, ms_sn_time
    ) ~ n_total + a + dist,
    data = df,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
}

main <- function() {
  reps <- as.integer(Sys.getenv("REPS", "100"))
  L <- as.integer(Sys.getenv("L_SPLIT", "3"))
  alpha <- as.numeric(Sys.getenv("ALPHA", "0.1"))
  B <- as.integer(Sys.getenv("B_BOOT", "200"))
  n_vals <- as.integer(strsplit(Sys.getenv("N_TOTALS", "400,800,1000,1600"), ",")[[1]])
  a_vals <- as.numeric(strsplit(Sys.getenv("AMPLITUDES", "0.5,0.625,0.75,0.875,1.0"), ",")[[1]])
  dist_vals <- strsplit(Sys.getenv("DISTS", "normal,t10"), ",")[[1]]

  dir.create("./results", recursive = TRUE, showWarnings = FALSE)

  message(sprintf("OPTICS version: %s", as.character(utils::packageVersion("OPTICS"))))
  message(sprintf("has OPTICS.ms: %s", exists("OPTICS.ms", mode = "function")))
  message(sprintf("reps=%d, L=%d, alpha=%.3f, B=%d", reps, L, alpha, B))

  detail <- list()
  idx <- 1L
  total_jobs <- length(n_vals) * length(a_vals) * length(dist_vals)
  for (n_total in n_vals) {
    for (a in a_vals) {
      for (dist in dist_vals) {
        message(sprintf("[%s] job %d/%d: n_total=%d, a=%.3f, dist=%s",
                        Sys.time(), idx, total_jobs, n_total, a, dist))
        detail[[idx]] <- run_setting(
          n_total = n_total, a = a, dist = dist, reps = reps,
          alpha = alpha, B = B, L = L
        )
        idx <- idx + 1L
      }
    }
  }

  detail_df <- do.call(rbind, detail)
  summary_df <- summarize_results(detail_df)

  write.csv(summary_df, "./results/mean_multiple_compare_summary.csv", row.names = FALSE)
  print(summary_df)
}

main()
