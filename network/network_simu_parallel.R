library(OPTICS)
library(parallel)

# -----------------------------------------------------------------------------
# Parallel simulation for network change-point model
# Z_i = Theta_k + W_i, with Theta_k from 3-community SBM.
# Memberships are reshuffled at each segment.
# -----------------------------------------------------------------------------

network_model <- function(n, p, cps, A, block_num = 3L) {
  Q1 <- A * matrix(
    c(0.6, 1.0, 0.6,
      1.0, 0.6, 0.5,
      0.6, 0.5, 0.6),
    nrow = 3, byrow = TRUE
  )
  Q2 <- A * matrix(
    c(0.6, 0.5, 0.6,
      0.5, 0.6, 1.0,
      0.6, 1.0, 0.6),
    nrow = 3, byrow = TRUE
  )

  cps_all <- c(0, cps, n)
  edge_num <- p * (p - 1) / 2
  d_c <- matrix(0, nrow = edge_num, ncol = n)

  for (seg in seq_len(length(cps_all) - 1)) {
    start <- cps_all[seg] + 1
    end <- cps_all[seg + 1]
    seg_n <- end - start + 1
    Q <- if (seg %% 2 == 1) Q1 else Q2

    membership <- sample.int(block_num, size = p, replace = TRUE)

    prob_vec <- numeric(edge_num)
    idx <- 1L
    for (u in 2:p) {
      for (v in 1:(u - 1)) {
        prob_vec[idx] <- Q[membership[u], membership[v]]
        idx <- idx + 1L
      }
    }

    obs_mat <- t(vapply(
      prob_vec,
      FUN = function(prob) rbinom(seg_n, size = 1, prob = prob),
      FUN.VALUE = integer(seg_n)
    ))
    d_c[, start:end] <- obs_mat
  }

  d_c
}

run_one_rep_network <- function(n, p, cps, A, model, true_model) {
  X <- network_model(n = n, p = p, cps = cps, A = A)
  cv <- OPTICS.network(X, model = model)

  c(
    OPTICS_cov = as.numeric(true_model %in% cv$A),
    OPTICS_len = length(cv$A),
    A_COPSS = as.numeric(true_model %in% c(cv$model.min - 1, cv$model.min, cv$model.min + 1)),
    COPSS = as.numeric(true_model == cv$model.min)
  )
}

run_setting_parallel_network <- function(n, p, cps, A_grid, simu_time = 100, cl) {
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
        run_one_rep_network(n = n, p = p, cps = cps, A = A, model = model, true_model = true_model)
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
# Experiment setting
# -----------------------------------------------------------------------------

set.seed(1234)
n <- 1000
p <- 5
cps <- c(200, 400, 600, 800)
A_grid <- seq(0.5, 1.0, by = 0.1)
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
  varlist = c("network_model", "run_one_rep_network"),
  envir = environment()
)

clusterSetRNGStream(cl, iseed = 1234)

dir.create("./results", recursive = TRUE, showWarnings = FALSE)
res_network <- run_setting_parallel_network(n = n, p = p, cps = cps, A_grid = A_grid, simu_time = simu_time, cl = cl)
write.csv(res_network, "./results/network.csv", row.names = FALSE)
