library(OPTICS)
library(parallel)

if (!exists("OPTICS.ms", mode = "function")) {
  stop("OPTICS.ms is not available in the loaded OPTICS package. Install the updated OPTICS package first.")
}

# -----------------------------------------------------------------------------
# 1D mean-change data generator
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# One replicate: base OPTICS vs multiple-splitting OPTICS under BS/SN
# -----------------------------------------------------------------------------
run_one_rep <- function(i, n_total, cps, a, dist, model_bs, model_sn, true_k, alpha, B, L, base_seed = 1234L) {
  set.seed(base_seed + 100000L * as.integer(n_total) + 1000L * as.integer(100 * a) + 100L * i + if (dist == "t10") 1L else 0L)
  X <- simulate_mean_1d(n_total = n_total, cps = cps, a = a, dist = dist)

  fit_base_bs <- tryCatch(
    OPTICS.mean(X, model = model_bs, method = "BinSeg", alpha = alpha, B = B),
    error = function(e) NULL
  )
  fit_base_sn <- tryCatch(
    OPTICS.mean(X, model = model_sn, method = "SegNeigh", alpha = alpha, B = B),
    error = function(e) NULL
  )
  fit_ms_bs <- tryCatch(
    OPTICS.ms(X, model = model_bs, engine = "mean", method = "BinSeg", alpha = alpha, B = B, L = L, trim_incomplete = TRUE),
    error = function(e) NULL
  )
  fit_ms_sn <- tryCatch(
    OPTICS.ms(X, model = model_sn, engine = "mean", method = "SegNeigh", alpha = alpha, B = B, L = L, trim_incomplete = TRUE),
    error = function(e) NULL
  )

  c(
    base_bs_cov = if (is.null(fit_base_bs)) NA_real_ else as.numeric(true_k %in% fit_base_bs$A),
    base_bs_len = if (is.null(fit_base_bs)) NA_real_ else length(fit_base_bs$A),
    ms_bs_cov = if (is.null(fit_ms_bs)) NA_real_ else as.numeric(true_k %in% fit_ms_bs$A_MS),
    ms_bs_len = if (is.null(fit_ms_bs)) NA_real_ else length(fit_ms_bs$A_MS),
    base_sn_cov = if (is.null(fit_base_sn)) NA_real_ else as.numeric((true_k + 1) %in% fit_base_sn$A),
    base_sn_len = if (is.null(fit_base_sn)) NA_real_ else length(fit_base_sn$A),
    ms_sn_cov = if (is.null(fit_ms_sn)) NA_real_ else as.numeric((true_k + 1) %in% fit_ms_sn$A_MS),
    ms_sn_len = if (is.null(fit_ms_sn)) NA_real_ else length(fit_ms_sn$A_MS)
  )
}

# -----------------------------------------------------------------------------
# Parallel run for one distribution and one n_total over amplitude grid
# -----------------------------------------------------------------------------
run_setting_parallel <- function(n_total, a_grid, dist, simu_time, alpha, B, L, cl) {
  cps <- round(n_total * c(0.2, 0.4, 0.6, 0.8))
  cps <- sort(unique(pmax(2L, pmin(n_total - 2L, as.integer(cps)))))
  true_k <- length(cps)

  k_max <- max(2L, ceiling(log(n_total / 2)))
  model_bs <- seq_len(k_max)
  model_sn <- model_bs + 1

  res <- matrix(NA_real_, nrow = length(a_grid), ncol = 9L)
  colnames(res) <- c(
    "a",
    "base_BS_cov", "base_BS_len",
    "MS_BS_cov", "MS_BS_len",
    "base_SN_cov", "base_SN_len",
    "MS_SN_cov", "MS_SN_len"
  )

  for (j in seq_along(a_grid)) {
    a <- a_grid[j]
    message(sprintf("[%s] n_total=%d, dist=%s, a=%.3f (%d/%d)",
                    Sys.time(), n_total, dist, a, j, length(a_grid)))

    rep_list <- parLapply(
      cl = cl,
      X = seq_len(simu_time),
      fun = function(i) {
        run_one_rep(
          i = i, n_total = n_total, cps = cps, a = a, dist = dist,
          model_bs = model_bs, model_sn = model_sn, true_k = true_k, alpha = alpha, B = B, L = L
        )
      }
    )

    rep_mat <- do.call(rbind, rep_list)
    res[j, ] <- c(a, colMeans(rep_mat, na.rm = TRUE))
    print(res[j, ])
  }

  as.data.frame(res)
}

# -----------------------------------------------------------------------------
# Experiments
# -----------------------------------------------------------------------------
n_total <- 1000
cps <- c(200, 400, 600, 800)
a_grid <- seq(0.5, 1, length.out = 5)
simu_time <- 100
alpha <- 0.1
B <- 200
L <- 3

slurm_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
if (!is.na(slurm_cpus) && slurm_cpus > 0) {
  num_cores <- slurm_cpus
} else {
  num_cores <- max(1L, detectCores(logical = TRUE) - 1L)
}

cl <- makeCluster(num_cores)
on.exit(stopCluster(cl), add = TRUE)

clusterEvalQ(cl, {
  library(OPTICS)
  NULL
})

clusterExport(
  cl,
  varlist = c("simulate_mean_1d", "run_one_rep"),
  envir = environment()
)

clusterSetRNGStream(cl, iseed = 1234)

dir.create("./results", recursive = TRUE, showWarnings = FALSE)
message(sprintf("OPTICS version: %s", as.character(utils::packageVersion("OPTICS"))))
message(sprintf("n_total=%d, reps=%d, alpha=%.3f, B=%d, L=%d, cores=%d", n_total, simu_time, alpha, B, L, num_cores))

res_normal <- run_setting_parallel(
  n_total = n_total, a_grid = a_grid, dist = "normal",
  simu_time = simu_time, alpha = alpha, B = B, L = L, cl = cl
)
write.csv(res_normal, "./results/one_d_normal_multiple.csv", row.names = FALSE)

res_t10 <- run_setting_parallel(
  n_total = n_total, a_grid = a_grid, dist = "t10",
  simu_time = simu_time, alpha = alpha, B = B, L = L, cl = cl
)
write.csv(res_t10, "./results/one_d_t10_multiple.csv", row.names = FALSE)
