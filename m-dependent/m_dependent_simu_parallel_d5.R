library(OPTICS)
library(parallel)

# -----------------------------------------------------------------------------
# Parallel simulation for m-dependent multivariate mean-change model (d = 5)
# y_i = mu_k + eps_i, eps_i = sum_{l=1}^M phi_l * eta_{i+l},
# eta_t ~ N(0, Sigma), phi_l = sqrt(1/M).
# -----------------------------------------------------------------------------

MCP_model_mdep <- function(n, p, cps, a, sparsity, M = 3L, rho = 0, dense = FALSE) {
  if (sparsity > p) warning("sparsity is larger than dimension")
  if (!is.numeric(M) || M <= 0 || M != as.integer(M)) stop("'M' must be a positive integer.")

  Sigma <- if (dense) {
    s <- matrix(rho, nrow = p, ncol = p)
    diag(s) <- 1
    s
  } else {
    toeplitz(rho^(0:(p - 1)))
  }

  cps_all <- c(0, cps, n)
  mean_mat <- matrix(0, nrow = p, ncol = n)
  for (seg in seq_len(length(cps_all) - 1)) {
    start <- cps_all[seg] + 1
    end <- cps_all[seg + 1]
    seg_mean <- matrix(0, nrow = p, ncol = 1)
    seg_mean[1:sparsity, ] <- ifelse(seg %% 2 == 1, a, -a)
    mean_mat[, start:end] <- seg_mean
  }

  eta <- mvtnorm::rmvnorm(n + M, sigma = Sigma)
  phi <- rep(sqrt(1 / M), M)

  eps <- matrix(0, nrow = n, ncol = p)
  for (i in seq_len(n)) {
    idx <- i + seq_len(M)
    eps[i, ] <- colSums(eta[idx, , drop = FALSE] * phi)
  }

  mean_mat + t(eps)
}

run_one_rep_mdep_d5 <- function(n, p, cps, a, model, true_model, M) {
  X <- MCP_model_mdep(n = n, p = p, cps = cps, a = a, sparsity = p, M = M)

  # m-dependent OPTICS
  m_cv_bs <- OPTICS.m_dependent(X, model = model, method = "BinSeg", M = M)
  m_cv_sn <- OPTICS.m_dependent(X, model = model + 1, method = "SegNeigh", M = M)

  # baseline (ignoring m-dependence): multivariate mean OPTICS
  cv_bs <- OPTICS.dmean(X, model = model, method = "BinSeg")
  cv_sn <- OPTICS.dmean(X, model = model + 1, method = "SegNeigh")

  c(
    M_OPTICS_BS_cov = as.numeric(true_model %in% m_cv_bs$A),
    M_OPTICS_BS_len = length(m_cv_bs$A),
    M_OPTICS_SN_cov = as.numeric((true_model + 1) %in% m_cv_sn$A),
    M_OPTICS_SN_len = length(m_cv_sn$A),
    OPTICS_BS_cov = as.numeric(true_model %in% cv_bs$A),
    OPTICS_BS_len = length(cv_bs$A),
    OPTICS_SN_cov = as.numeric((true_model + 1) %in% cv_sn$A),
    OPTICS_SN_len = length(cv_sn$A),
    A_COPSS_BS = as.numeric(true_model %in% c(cv_bs$model.min - 1, cv_bs$model.min, cv_bs$model.min + 1)),
    A_COPSS_SN = as.numeric((true_model + 1) %in% c(cv_sn$model.min - 1, cv_sn$model.min, cv_sn$model.min + 1)),
    COPSS_BS = as.numeric(true_model == cv_bs$model.min),
    COPSS_SN = as.numeric((true_model + 1) == cv_sn$model.min)
  )
}

run_setting_parallel_mdep_d5 <- function(n, p, cps, a_grid, M, simu_time = 100, cl) {
  model <- 1:ceiling(log(n))
  true_model <- length(cps)
  res <- matrix(NA_real_, nrow = length(a_grid), ncol = 13)

  for (iter in seq_along(a_grid)) {
    a <- a_grid[iter]
    message(sprintf("[%s] d=%d, M=%d, amplitude=%.3f (%d/%d)", Sys.time(), p, M, a, iter, length(a_grid)))

    rep_list <- parLapply(
      cl = cl,
      X = seq_len(simu_time),
      fun = function(i, n, p, cps, a, model, true_model, M) {
        run_one_rep_mdep_d5(n = n, p = p, cps = cps, a = a, model = model, true_model = true_model, M = M)
      },
      n = n,
      p = p,
      cps = cps,
      a = a,
      model = model,
      true_model = true_model,
      M = M
    )

    rep_mat <- do.call(rbind, rep_list)
    res[iter, ] <- c(a, colMeans(rep_mat))
    print(res[iter, ])
  }

  colnames(res) <- c(
    "a",
    "M_OPTICS(BS)_cov", "M_OPTICS(BS)_len", "M_OPTICS(SN)_cov", "M_OPTICS(SN)_len",
    "OPTICS(BS)_cov", "OPTICS(BS)_len", "OPTICS(SN)_cov", "OPTICS(SN)_len",
    "A_COPSS(BS)", "A_COPSS(SN)", "COPSS(BS)", "COPSS(SN)"
  )
  res
}

# -----------------------------------------------------------------------------
# Experiment setting
# -----------------------------------------------------------------------------

set.seed(1234)
n <- 1000
p <- 5
cps <- c(200, 400, 600, 800)
M <- 3
a_grid <- seq(0.5, 1, length.out = 5)
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

clusterExport(
  cl,
  varlist = c("MCP_model_mdep", "run_one_rep_mdep_d5"),
  envir = environment()
)

clusterSetRNGStream(cl, iseed = 1234)

dir.create("./results", recursive = TRUE, showWarnings = FALSE)
res_mdep_d5 <- run_setting_parallel_mdep_d5(
  n = n, p = p, cps = cps, a_grid = a_grid, M = M, simu_time = simu_time, cl = cl
)
write.csv(res_mdep_d5, file = "./results/m_dependent_d5.csv", row.names = FALSE)
