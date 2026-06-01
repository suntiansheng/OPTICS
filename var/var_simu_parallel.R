library(OPTICS)
library(FDRSeg)
library(parallel)

# -----------------------------------------------------------------------------
# Parallel simulation for variance change-point model
# y_i = sigma_k * eps_i, eps_i ~ N(0, 0.25), sigma_1 = 1
# sigma_{k+1} / sigma_k alternates between A and 1/A.
# -----------------------------------------------------------------------------

variance_model <- function(n, cps, a) {
  cps_all <- c(0, cps, n)
  sigma <- numeric(n)
  sigma_level <- 1

  for (seg in seq_len(length(cps_all) - 1)) {
    start <- cps_all[seg] + 1
    end <- cps_all[seg + 1]
    sigma[start:end] <- sigma_level

    if (seg <= length(cps)) {
      if (seg %% 2 == 1) {
        sigma_level <- sigma_level * a
      } else {
        sigma_level <- sigma_level / a
      }
    }
  }

  eps_star <- rnorm(n, mean = 0, sd = sqrt(0.25))
  y <- sigma * eps_star
  matrix(y, nrow = 1)
}

run_one_rep_var <- function(n, cps, a, model, true_model) {
  X <- variance_model(n = n, cps = cps, a = a)
  X <- log(X^2)

  cv_bs <- OPTICS.mean(X, model = model, method = "BinSeg")
  cv_sn <- OPTICS.mean(X, model = model + 1, method = "SegNeigh")

  fdrseg_fit <- fdrseg(as.numeric(X))
  fdrseg_k <- length(fdrseg_fit$left) - 1

  smuce_fit <- smuce(as.numeric(X))
  smuce_k <- length(smuce_fit$left) - 1

  c(
    OPTICS_BS_cov = as.numeric(true_model %in% cv_bs$A),
    OPTICS_BS_len = length(cv_bs$A),
    OPTICS_SN_cov = as.numeric(true_model %in% (cv_sn$A - 1)),
    OPTICS_SN_len = length(cv_sn$A),
    A_COPSS_BS = as.numeric(true_model %in% c(cv_bs$model.min - 1, cv_bs$model.min, cv_bs$model.min + 1)),
    A_COPSS_SN = as.numeric(true_model %in% c(cv_sn$model.min - 2, cv_sn$model.min - 1, cv_sn$model.min)),
    COPSS_BS = as.numeric(true_model == cv_bs$model.min),
    COPSS_SN = as.numeric(true_model == cv_sn$model.min - 1),
    A_FDRseg = as.numeric(true_model %in% c(fdrseg_k - 1, fdrseg_k, fdrseg_k + 1)),
    A_SMUCE = as.numeric(true_model %in% c(smuce_k - 1, smuce_k, smuce_k + 1)),
    fdr_seg = as.numeric(true_model == fdrseg_k),
    smuce = as.numeric(true_model == smuce_k)
  )
}

run_setting_parallel_var <- function(n, cps, a_grid, simu_time = 100, cl) {
  model <- 1:ceiling(log(n))
  true_model <- length(cps)
  res <- matrix(NA_real_, nrow = length(a_grid), ncol = 13)

  for (iter in seq_along(a_grid)) {
    a <- a_grid[iter]
    message(sprintf("[%s] amplitude=%.3f (%d/%d)", Sys.time(), a, iter, length(a_grid)))

    rep_list <- parLapply(
      cl = cl,
      X = seq_len(simu_time),
      fun = function(i) {
        run_one_rep_var(n = n, cps = cps, a = a, model = model, true_model = true_model)
      }
    )
    rep_mat <- do.call(rbind, rep_list)
    res[iter, ] <- c(a, colMeans(rep_mat))
    print(res[iter, ])
  }

  colnames(res) <- c(
    "a", "OPTICS(BS)_cov", "OPTICS(BS)_len", "OPTICS(SN)_cov", "OPTICS(SN)_len",
    "A_COPSS(BS)", "A_COPSS(SN)", "COPSS(BS)", "COPSS(SN)",
    "A_FDRseg", "A_SMUCE", "fdr_seg", "smuce"
  )
  res
}

# -----------------------------------------------------------------------------
# Experiment setting from paper
# -----------------------------------------------------------------------------

set.seed(1234)
n <- 1000
cps <- c(200, 400, 600, 800)
a_grid <- seq(2, 6, length.out = 5)
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
  library(FDRSeg)
  NULL
})

clusterExport(
  cl,
  varlist = c("variance_model", "run_one_rep_var"),
  envir = environment()
)

clusterSetRNGStream(cl, iseed = 1234)

dir.create("./results", recursive = TRUE, showWarnings = FALSE)
res_var <- run_setting_parallel_var(n = n, cps = cps, a_grid = a_grid, simu_time = simu_time, cl = cl)
write.csv(res_var, file = "./results/var_normal.csv", row.names = FALSE)
