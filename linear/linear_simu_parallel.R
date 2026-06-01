library(OPTICS)
library(parallel)

# -----------------------------------------------------------------------------
# Parallel simulation for linear-model coefficient structural-breaks
# -----------------------------------------------------------------------------

Regression_model <- function(n, p, cps, a, dis = c("normal", "t")) {
  dis <- match.arg(dis)

  # x_i ~ N(0, I_d)
  X <- matrix(rnorm(n * p), ncol = p)

  # beta_k = (-1)^(k-1) * A * 1_d
  y <- numeric(n)
  cps_all <- c(0, cps, n)
  for (seg in seq_len(length(cps_all) - 1)) {
    start <- cps_all[seg] + 1
    end <- cps_all[seg + 1]
    beta <- rep(ifelse(seg %% 2 == 1, a, -a), p)
    y[start:end] <- X[start:end, , drop = FALSE] %*% beta
  }

  # eps_i ~ N(0, 1) or t(10)
  eps <- if (dis == "normal") rnorm(n) else rt(n, df = 10)
  y <- y + eps

  list(X = X, y = matrix(y, ncol = 1))
}

run_one_rep_linear <- function(n, p, cps, a, model, true_model, dis) {
  dat <- Regression_model(n = n, p = p, cps = cps, a = a, dis = dis)
  cv_bs <- OPTICS.linear(dat$X, dat$y, model = model, method = "BinSeg")
  cv_sn <- OPTICS.linear(dat$X, dat$y, model = model, method = "SegNeigh")

  c(
    OPTICS_BS_cov = as.numeric(true_model %in% cv_bs$A),
    OPTICS_BS_len = length(cv_bs$A),
    OPTICS_SN_cov = as.numeric(true_model %in% cv_sn$A),
    OPTICS_SN_len = length(cv_sn$A),
    A_COPSS_BS = as.numeric(true_model %in% c(cv_bs$model.min - 1, cv_bs$model.min, cv_bs$model.min + 1)),
    A_COPSS_SN = as.numeric(true_model %in% c(cv_sn$model.min - 1, cv_sn$model.min, cv_sn$model.min + 1)),
    COPSS_BS = as.numeric(true_model == cv_bs$model.min),
    COPSS_SN = as.numeric(true_model == cv_sn$model.min)
  )
}

run_setting_parallel_linear <- function(n, p, cps, a_grid, dis, simu_time = 100, cl) {
  model <- 1:ceiling(log(n))
  true_model <- length(cps)
  res <- matrix(NA_real_, nrow = length(a_grid), ncol = 9)

  for (iter in seq_along(a_grid)) {
    a <- a_grid[iter]
    message(sprintf("[%s] dist=%s, amplitude=%.3f (%d/%d)", Sys.time(), dis, a, iter, length(a_grid)))

    rep_list <- parLapply(
      cl = cl,
      X = seq_len(simu_time),
      fun = function(i, n, p, cps, a, model, true_model, dis) {
        run_one_rep_linear(n, p, cps, a, model, true_model, dis)
      },
      n = n,
      p = p,
      cps = cps,
      a = a,
      model = model,
      true_model = true_model,
      dis = dis
    )
    rep_mat <- do.call(rbind, rep_list)
    res[iter, ] <- c(a, colMeans(rep_mat))
    print(res[iter, ])
  }

  colnames(res) <- c(
    "a", "OPTICS(BS)_cov", "OPTICS(BS)_len", "OPTICS(SN)_cov", "OPTICS(SN)_len",
    "A_COPSS(BS)", "A_COPSS(SN)", "COPSS(BS)", "COPSS(SN)"
  )
  res
}

# -----------------------------------------------------------------------------
# Experiments
# -----------------------------------------------------------------------------

set.seed(1234)
n <- 1000
p <- 5
cps <- c(200, 400, 600, 800)
a_grid <- seq(0.2, 0.3, length.out = 5)
simu_time <- 100

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
clusterExport(cl, varlist = c("Regression_model", "run_one_rep_linear"), envir = environment())
clusterSetRNGStream(cl, iseed = 1234)

dir.create("./results", recursive = TRUE, showWarnings = FALSE)

res_normal <- run_setting_parallel_linear(
  n = n, p = p, cps = cps, a_grid = a_grid, dis = "normal", simu_time = simu_time, cl = cl
)
write.csv(res_normal, file = "./results/linear_normal.csv", row.names = FALSE)

res_t <- run_setting_parallel_linear(
  n = n, p = p, cps = cps, a_grid = a_grid, dis = "t", simu_time = simu_time, cl = cl
)
write.csv(res_t, file = "./results/linear_t.csv", row.names = FALSE)
