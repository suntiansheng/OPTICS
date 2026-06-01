library(OPTICS)
library(FDRSeg)
library(parallel)

# -----------------------------------------------------------------------------
# Parallel simulation for multiple mean change-point models
# -----------------------------------------------------------------------------

MCP_model <- function(n, p, cps, a, sparsity, dis = c("normal", "t"), rho = 0, dense = FALSE) {
  dis <- match.arg(dis)
  if (sparsity > p) warning("sparsity is larger than dimension")

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

  eps <- if (dis == "normal") {
    mvtnorm::rmvnorm(n, sigma = Sigma)
  } else {
    mvtnorm::rmvt(n, sigma = Sigma, df = 10)
  }

  mean_mat + t(eps)
}

run_one_rep <- function(n, p, cps, a, model, true_model, dis, include_fdr_smuce = TRUE) {
  X <- MCP_model(n = n, p = p, cps = cps, a = a, sparsity = p, dis = dis)

  cv_bs <- OPTICS.mean(X, model = model, method = "BinSeg")
  cv_sn <- OPTICS.mean(X, model = model + 1, method = "SegNeigh")

  out <- c(
    OPTICS_BS_cov = as.numeric(true_model %in% cv_bs$A),
    OPTICS_BS_len = length(cv_bs$A),
    OPTICS_SN_cov = as.numeric(true_model %in% (cv_sn$A - 1)),
    OPTICS_SN_len = length(cv_sn$A),
    A_COPSS_BS = as.numeric(true_model %in% c(cv_bs$model.min - 1, cv_bs$model.min, cv_bs$model.min + 1)),
    A_COPSS_SN = as.numeric(true_model %in% c(cv_sn$model.min - 2, cv_sn$model.min - 1, cv_sn$model.min)),
    COPSS_BS = as.numeric(true_model == cv_bs$model.min),
    COPSS_SN = as.numeric(true_model == cv_sn$model.min - 1)
  )

  if (include_fdr_smuce) {
    fdrseg_fit <- fdrseg(as.numeric(X))
    fdrseg_k <- length(fdrseg_fit$left) - 1

    smuce_fit <- smuce(as.numeric(X))
    smuce_k <- length(smuce_fit$left) - 1

    out <- c(
      out,
      A_FDRseg = as.numeric(true_model %in% c(fdrseg_k - 1, fdrseg_k, fdrseg_k + 1)),
      A_SMUCE = as.numeric(true_model %in% c(smuce_k - 1, smuce_k, smuce_k + 1)),
      fdr_seg = as.numeric(true_model == fdrseg_k),
      smuce = as.numeric(true_model == smuce_k)
    )
  }

  out
}

run_setting_parallel <- function(n, p, cps, a_grid, dis, simu_time = 100, cl) {
  k_max <- ceiling(log(n))
  model <- 1:k_max
  true_model <- length(cps)
  include_fdr_smuce <- (p == 1)

  res <- matrix(NA_real_, nrow = length(a_grid), ncol = if (include_fdr_smuce) 13 else 9)

  for (iter in seq_along(a_grid)) {
    a <- a_grid[iter]
    message(sprintf("[%s] p=%d, dist=%s, amplitude=%.3f (%d/%d)", Sys.time(), p, dis, a, iter, length(a_grid)))

    rep_list <- parLapply(
      cl = cl,
      X = seq_len(simu_time),
      fun = function(i) {
        run_one_rep(
          n = n, p = p, cps = cps, a = a, model = model,
          true_model = true_model, dis = dis, include_fdr_smuce = include_fdr_smuce
        )
      }
    )
    rep_mat <- do.call(rbind, rep_list)
    res[iter, ] <- c(a, colMeans(rep_mat))
    print(res[iter, ])
  }

  if (include_fdr_smuce) {
    colnames(res) <- c(
      "a", "OPTICS(BS)_cov", "OPTICS(BS)_len", "OPTICS(SN)_cov", "OPTICS(SN)_len",
      "A_COPSS(BS)", "A_COPSS(SN)", "COPSS(BS)", "COPSS(SN)",
      "A_FDRseg", "A_SMUCE", "fdr_seg", "smuce"
    )
  } else {
    colnames(res) <- c(
      "a", "OPTICS(BS)_cov", "OPTICS(BS)_len", "OPTICS(SN)_cov", "OPTICS(SN)_len",
      "A_COPSS(BS)", "A_COPSS(SN)", "COPSS(BS)", "COPSS(SN)"
    )
  }

  res
}

# -----------------------------------------------------------------------------
# Experiments
# -----------------------------------------------------------------------------

n <- 1000
cps <- c(200, 400, 600, 800)
simu_time <- 50

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
  library(FDRSeg)
  NULL
})

clusterExport(
  cl,
  varlist = c(
    "MCP_model", "run_one_rep"
  ),
  envir = environment()
)

# Reproducible RNG streams across workers.
clusterSetRNGStream(cl, iseed = 1234)

dir.create("./results", recursive = TRUE, showWarnings = FALSE)

# a_grid <- seq(0.25, 0.5, length.out = 5)

# res_1d_normal <- run_setting_parallel(
#   n = n, p = 1, cps = cps, a_grid = a_grid, dis = "normal", simu_time = simu_time, cl = cl
# )
# write.csv(res_1d_normal, file = "./results/one_d_normal.csv", row.names = FALSE)

# res_1d_t <- run_setting_parallel(
#   n = n, p = 1, cps = cps, a_grid = a_grid, dis = "t", simu_time = simu_time, cl = cl
# )
# write.csv(res_1d_t, file = "./results/one_d_t.csv", row.names = FALSE)

a_grid <- seq(1, 2, length.out = 5)
res_5d_normal <- run_setting_parallel(
  n = n, p = 5, cps = cps, a_grid = a_grid, dis = "normal", simu_time = simu_time, cl = cl
)
write.csv(res_5d_normal, file = "./results/multi_d_norm.csv", row.names = FALSE)

res_5d_t <- run_setting_parallel(
  n = n, p = 5, cps = cps, a_grid = a_grid, dis = "t", simu_time = simu_time, cl = cl
)
write.csv(res_5d_t, file = "./results/multi_d_t.csv", row.names = FALSE)
