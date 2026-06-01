library(OPTICS)
library(parallel)

# -----------------------------------------------------------------------------
# Parallel simulation for heavy-tail mean change-point model
# Same mean-change setup as Section 5.1, but errors are t(df = 1).
# -----------------------------------------------------------------------------

MCP_model_heavy <- function(n, p, cps, a, sparsity, rho = 0, dense = FALSE) {
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

  eps <- mvtnorm::rmvt(n, sigma = Sigma, df = 1)
  mean_mat + t(eps)
}

run_one_rep_heavy <- function(n, p, cps, a, model, true_model) {
  X <- MCP_model_heavy(n = n, p = p, cps = cps, a = a, sparsity = p)

  ro_cv_bs <- OPTICS.heavy_tail(X, model = model, method = "BinSeg")
  ro_cv_sn <- OPTICS.heavy_tail(X, model = model + 1, method = "SegNeigh")

  cv_bs <- OPTICS.mean(X, model = model, method = "BinSeg")
  cv_sn <- OPTICS.mean(X, model = model + 1, method = "SegNeigh")

  c(
    RO_OPTICS_BS_cov = as.numeric(true_model %in% ro_cv_bs$A),
    RO_OPTICS_BS_len = length(ro_cv_bs$A),
    RO_OPTICS_SN_cov = as.numeric(true_model %in% (ro_cv_sn$A - 1)),
    RO_OPTICS_SN_len = length(ro_cv_sn$A),
    OPTICS_BS_cov = as.numeric(true_model %in% cv_bs$A),
    OPTICS_BS_len = length(cv_bs$A),
    OPTICS_SN_cov = as.numeric(true_model %in% (cv_sn$A - 1)),
    OPTICS_SN_len = length(cv_sn$A),
    A_COPSS_BS = as.numeric(true_model %in% c(cv_bs$model.min - 1, cv_bs$model.min, cv_bs$model.min + 1)),
    A_COPSS_SN = as.numeric(true_model %in% c(cv_sn$model.min - 2, cv_sn$model.min - 1, cv_sn$model.min)),
    COPSS_BS = as.numeric(true_model == cv_bs$model.min),
    COPSS_SN = as.numeric(true_model == cv_sn$model.min - 1)
  )
}

run_setting_parallel_heavy <- function(n, p, cps, a_grid, simu_time = 100, cl) {
  model <- 1:ceiling(log(n))
  true_model <- length(cps)
  res <- matrix(NA_real_, nrow = length(a_grid), ncol = 13)

  for (iter in seq_along(a_grid)) {
    a <- a_grid[iter]
    message(sprintf("[%s] amplitude=%.3f (%d/%d)", Sys.time(), a, iter, length(a_grid)))

    rep_list <- parLapply(
      cl = cl,
      X = seq_len(simu_time),
      fun = function(i, n, p, cps, a, model, true_model) {
        run_one_rep_heavy(n = n, p = p, cps = cps, a = a, model = model, true_model = true_model)
      },
      n = n,
      p = p,
      cps = cps,
      a = a,
      model = model,
      true_model = true_model
    )

    rep_mat <- do.call(rbind, rep_list)
    res[iter, ] <- c(a, colMeans(rep_mat))
    print(res[iter, ])
  }

  colnames(res) <- c(
    "a", "RO_OPTICS(BS)_cov", "RO_OPTICS(BS)_len", "RO_OPTICS(SN)_cov", "RO_OPTICS(SN)_len",
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
p <- 1
cps <- c(200, 400, 600, 800)
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
  varlist = c("MCP_model_heavy", "run_one_rep_heavy"),
  envir = environment()
)

clusterSetRNGStream(cl, iseed = 1234)

dir.create("./results", recursive = TRUE, showWarnings = FALSE)
res_heavy <- run_setting_parallel_heavy(n = n, p = p, cps = cps, a_grid = a_grid, simu_time = simu_time, cl = cl)
write.csv(res_heavy, file = "./results/heavy_one_d.csv", row.names = FALSE)
