library(OPTICS)
library(MASS)
library(parallel)
library(changepoints)

# -----------------------------------------------------------------------------
# Parallel simulation for covariance change-point model
# z_i = Theta_k^{1/2} eta_i, eta_i ~ N(0, I_d)
# Theta_1 = I_d, Theta_2 = {A^{|u-v|}}_{u,v=1}^d
# -----------------------------------------------------------------------------

covariance_model <- function(n, p, cps, A) {
  theta_1 <- diag(p)
  theta_2 <- toeplitz(A^(0:(p - 1)))

  cps_all <- c(0, cps, n)
  Z <- matrix(0, nrow = p, ncol = n)

  for (seg in seq_len(length(cps_all) - 1)) {
    start <- cps_all[seg] + 1
    end <- cps_all[seg + 1]
    theta <- if (seg %% 2 == 1) theta_1 else theta_2
    Z[, start:end] <- t(MASS::mvrnorm(end - start + 1, mu = rep(0, p), Sigma = theta))
  }

  Z
}

run_one_rep_cov <- function(n, p, cps, A, model, true_model) {
  X <- covariance_model(n = n, p = p, cps = cps, A = A)
  cv <- OPTICS.covariance(X, model = model)

  c(
    OPTICS_cov = as.numeric(true_model %in% cv$A),
    OPTICS_len = length(cv$A),
    A_COPSS = as.numeric(true_model %in% c(cv$model.min - 1, cv$model.min, cv$model.min + 1)),
    COPSS = as.numeric(true_model == cv$model.min)
  )
}

run_setting_parallel_cov <- function(n, p, cps, A_grid, simu_time = 100, cl) {
  model <- 1:ceiling(log(n))
  true_model <- length(cps)
  res <- matrix(NA_real_, nrow = length(A_grid), ncol = 5)

  for (iter in seq_along(A_grid)) {
    A <- A_grid[iter]
    message(sprintf("[%s] amplitude=%.3f (%d/%d)", Sys.time(), A, iter, length(A_grid)))

    rep_list <- parLapply(
      cl = cl,
      X = seq_len(simu_time),
      fun = function(i, n, p, cps, A, model, true_model) {
        run_one_rep_cov(n = n, p = p, cps = cps, A = A, model = model, true_model = true_model)
      },
      n = n,
      p = p,
      cps = cps,
      A = A,
      model = model,
      true_model = true_model
    )

    rep_mat <- do.call(rbind, rep_list)
    res[iter, ] <- c(A, colMeans(rep_mat))
    print(res[iter, ])
  }

  colnames(res) <- c("a", "OPTICS_cov", "OPTICS_len", "A_COPSS", "COPSS")
  res
}

# -----------------------------------------------------------------------------
# Experiment setting from paper
# -----------------------------------------------------------------------------

set.seed(1234)
n <- 1000
p <- 5
cps <- c(200, 400, 600, 800)
A_grid <- seq(0.3, 0.5, length.out = 5)
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
  library(MASS)
  NULL
})

clusterExport(
  cl,
  varlist = c("covariance_model", "run_one_rep_cov"),
  envir = environment()
)

clusterSetRNGStream(cl, iseed = 1234)

dir.create("./results", recursive = TRUE, showWarnings = FALSE)
res_cov <- run_setting_parallel_cov(n = n, p = p, cps = cps, A_grid = A_grid, simu_time = simu_time, cl = cl)
write.csv(res_cov, "./results/covariance.csv", row.names = FALSE)
