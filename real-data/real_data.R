library(OPTICS)
library(changepoint)
library(ecp)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(pracma)

dir.create("./figure", recursive = TRUE, showWarnings = FALSE)
dir.create("./result", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Shared plot style (JCGS-friendly)
# -----------------------------------------------------------------------------

col_signal <- "#4D4D4D"
col_cp <- "#B2182B"
col_optics <- "#1f78b4"
col_random <- "#e66101"

paper_theme <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "grey55", fill = NA, linewidth = 0.4),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey96", color = "grey70"),
      strip.text = element_text(face = "bold"),
      legend.title = element_blank()
    )
}

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------

data(ACGH, package = "ecp")
X <- ACGH$data
X <- t(as.matrix(X))
X <- X[1:10, , drop = FALSE]
colnames(X) <- NULL
rownames(X) <- paste0("ID", seq_len(nrow(X)))

K_max <- ceiling(log(ncol(X)))
model <- seq(1, 3 * K_max, by = 3)

# -----------------------------------------------------------------------------
# 1) Moment bound plot (keep orders 2:100)
# -----------------------------------------------------------------------------

orders <- 2:100
n_sets <- nrow(X)
X_c <- X - rowMeans(X)

moment_m <- vapply(
  orders,
  FUN = function(k) rowMeans(abs(X_c)^k),
  FUN.VALUE = numeric(n_sets)
)
moment_m <- t(moment_m)
upper <- orders^(orders / 2)

df_mom <- as.data.frame(log10(moment_m)) |>
  setNames(rownames(X)) |>
  mutate(order = orders) |>
  pivot_longer(
    cols = -order,
    names_to = "dataset",
    values_to = "log_moment"
  ) |>
  drop_na(log_moment)

df_bound <- tibble(order = orders, log_upper = log10(upper))

plt_moment <- ggplot() +
  geom_line(
    data = df_mom,
    aes(x = order, y = log_moment, group = dataset),
    color = col_signal, alpha = 0.45, linewidth = 0.35
  ) +
  geom_line(
    data = df_bound,
    aes(x = order, y = log_upper),
    color = "black", linewidth = 0.9
  ) +
  labs(
    x = "Moment order k",
    y = expression(log[10] * " moment")
  ) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  paper_theme(12)

ggsave("./figure/moment_upper_bound.png", plt_moment, width = 8, height = 5, dpi = 320)
ggsave("./figure/moment_upper_bound.pdf", plt_moment, width = 8, height = 5)

# -----------------------------------------------------------------------------
# 2) Single-individual OPTICS sets and Hausdorff comparison
# -----------------------------------------------------------------------------

# Load precomputed confidence sets.
A_ls <- readRDS("./result/A_ls.RData")
if (length(A_ls) > nrow(X)) A_ls <- A_ls[seq_len(nrow(X))]
if (is.null(names(A_ls))) names(A_ls) <- rownames(X)

pairs_idx <- combn(seq_along(A_ls), 2)
A_dist <- apply(pairs_idx, 2, function(idx) hausdorff_dist(A_ls[[idx[1]]], A_ls[[idx[2]]]))

set.seed(1234)
random_ls <- lapply(A_ls, function(a_set) sample(model, length(a_set), replace = FALSE))
random_dist <- apply(pairs_idx, 2, function(idx) hausdorff_dist(random_ls[[idx[1]]], random_ls[[idx[2]]]))

df_dist <- bind_rows(
  data.frame(group = "OPTICS sets", value = A_dist),
  data.frame(group = "Random sets", value = random_dist)
)

plt_haus <- ggplot(df_dist, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.58, alpha = 0.7, outlier.shape = NA, color = "grey30") +
  geom_jitter(width = 0.10, size = 1.2, alpha = 0.20, color = "grey20") +
  scale_fill_manual(values = c("OPTICS sets" = col_optics, "Random sets" = col_random), guide = "none") +
  labs(y = "Hausdorff distance", x = NULL) +
  paper_theme(12)

ggsave("./figure/combined_boxplot.png", plt_haus, width = 7, height = 5, dpi = 320)
ggsave("./figure/combined_boxplot.pdf", plt_haus, width = 7, height = 5)

# -----------------------------------------------------------------------------
# 3) Separate plots (per individual), CPs selected using min(A_i)
# -----------------------------------------------------------------------------

cpt_ls <- lapply(seq_len(nrow(X)), function(i) {
  fit <- changepoint::cpt.mean(
    as.numeric(X[i, ]),
    penalty = "None",
    method = "BinSeg",
    Q = min(A_ls[[i]]),
    class = FALSE,
    minseglen = 20
  )
  fit[-length(fit)]
})
names(cpt_ls) <- rownames(X)
saveRDS(cpt_ls, file = "./result/cpt_ls.RData")

df_long <- as.data.frame(X) |>
  rownames_to_column("individual") |>
  pivot_longer(cols = -individual, names_to = "position_chr", values_to = "value") |>
  mutate(position = as.numeric(str_extract(position_chr, "\\d+"))) |>
  drop_na(position)

df_cp <- tibble(
  individual = rep(names(cpt_ls), lengths(cpt_ls)),
  cp = unlist(cpt_ls),
  khat = rep(vapply(A_ls, min, numeric(1)), lengths(cpt_ls))
)

id_levels <- rownames(X)
df_long <- df_long |>
  mutate(individual = factor(individual, levels = id_levels))
df_cp <- df_cp |>
  mutate(individual = factor(individual, levels = id_levels))

df_khat <- tibble(
  individual = factor(names(A_ls), levels = id_levels),
  khat = vapply(A_ls, min, numeric(1)),
  khat_label = paste0("hat(K)==", khat)
)

plt_sep <- ggplot(df_long, aes(x = position, y = value, group = individual)) +
  geom_line(color = col_signal, linewidth = 0.28, alpha = 0.9, na.rm = TRUE) +
  geom_vline(
    data = df_cp,
    aes(xintercept = cp),
    color = col_cp, linetype = "22", linewidth = 0.55, alpha = 0.9
  ) +
  geom_text(
    data = df_khat,
    aes(x = -Inf, y = Inf, label = khat_label),
    parse = TRUE, hjust = -0.05, vjust = 1.25, size = 2.8, inherit.aes = FALSE
  ) +
  facet_wrap(~ individual, ncol = 1, scales = "free_y") +
  labs(
    x = "Micro-array index",
    y = "Log-intensity ratio"
  ) +
  paper_theme(10) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_blank(),
    panel.spacing.y = unit(0.25, "lines")
  )

ggsave("./figure/SNP_line_changepoints1_10.pdf", plt_sep, width = 8.5, height = 12.5)
ggsave("./figure/SNP_line_changepoints1_10.png", plt_sep, width = 8.5, height = 12.5, dpi = 320)

# -----------------------------------------------------------------------------
# 4) Joint analysis (OPTICS.dmean with SegNeigh), unified style with separate
# -----------------------------------------------------------------------------

dfit <- OPTICS.dmean(X, model,delta=100)
khat_joint <- if (length(dfit$A) > 0) min(dfit$A) else dfit$model.min
print(dfit)
cat("khat_joint =", khat_joint, "\n")
saveRDS(dfit$A, file = "./result/A_joint_OPTICS.RData")

set.seed(1234)
wbs_fun <- getFromNamespace("network.detection.fun", "OPTICS")
cp_wbs <- sort(unique(as.integer(wbs_fun(X, num = khat_joint, delta = 20))))
cp_wbs_full <- sort(unique(as.integer(wbs_fun(X, num = max(model), delta = 20))))
cp_other <- setdiff(cp_wbs_full, cp_wbs)
saveRDS(cp_wbs, file = "./result/cp_joint_WBS.RData")

offset_map <- setNames(seq(0, by = 2, length.out = nrow(X)), rownames(X))
df_joint <- as.data.frame(X) |>
  rownames_to_column("individual") |>
  pivot_longer(cols = -individual, names_to = "position_chr", values_to = "value") |>
  mutate(
    position = as.numeric(str_extract(position_chr, "\\d+")),
    offset = offset_map[individual],
    value_off = value + offset
  ) |>
  drop_na(position)

cp_ticks <- sort(unique(c(cp_wbs, cp_other)))
base_ticks <- pretty(range(df_joint$position), n = 8)
joint_breaks <- base_ticks
joint_labels <- as.character(joint_breaks)

plt_joint <- ggplot(df_joint, aes(x = position, y = value_off, group = individual)) +
  geom_line(color = col_signal, linewidth = 0.28, alpha = 0.95) +
  geom_vline(
    xintercept = cp_other,
    color = "#f4a3a3", linetype = "22", linewidth = 0.45, alpha = 0.75
  ) +
  geom_vline(
    xintercept = cp_wbs,
    color = col_cp, linetype = "22", linewidth = 1.15, alpha = 0.98
  ) +
  geom_rug(
    data = data.frame(position = cp_other),
    aes(x = position),
    inherit.aes = FALSE,
    sides = "b",
    linewidth = 0.30,
    color = "#f4a3a3",
    alpha = 0.65
  ) +
  geom_rug(
    data = data.frame(position = cp_wbs),
    aes(x = position),
    inherit.aes = FALSE,
    sides = "b",
    linewidth = 0.55,
    color = col_cp,
    alpha = 0.90
  ) +
  labs(
    x = "Micro-array index",
    y = "Offset log-intensity ratio"
  ) +
  scale_x_continuous(breaks = joint_breaks, labels = joint_labels) +
  paper_theme(11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("./figure/joint.png", plt_joint, width = 12, height = 7, dpi = 320)
ggsave("./figure/joint.pdf", plt_joint, width = 12, height = 7)
