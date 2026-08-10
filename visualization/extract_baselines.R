library(dplyr)

setwd(".")

mean_ci <- function(x, confidence = 0.95) {
  x <- na.omit(x)
  n <- length(x)
  mean_x <- mean(x)
  std_err <- sd(x) / sqrt(n)
  error_margin <- qt((1 + confidence) / 2, df = n - 1) * std_err
  return(data.frame(mean = mean_x, ymin = mean_x - error_margin, ymax = mean_x + error_margin))
}

# ── IHC (Pathology) — 24 mpi ─────────────────────────────────────────────────
path_group_map <- list(
  "hipp_ant" = list(injection = "HIP", direction = "antero"),
  "hipp_ret" = list(injection = "HIP", direction = "retro"),
  "cp_ant"   = list(injection = "CP",  direction = "antero"),
  "cp_ret"   = list(injection = "CP",  direction = "retro")
)

ihc_time <- "24"

ihc_results <- lapply(names(path_group_map), function(group_name) {
  m  <- path_group_map[[group_name]]
  fp <- paste0("../sir_result_csvs/baselines/Total_Pathology_", ihc_time,
               "_MPI", m$injection, "_injection_rebuilt_rvm_", m$direction, ".csv")
  data <- read.csv(fp, header = FALSE)[, -1]
  colnames(data) <- c("Spearman_r", "v", "spread_rate", "injection_amount", "peak_timestep")
  baseline <- mean_ci(data$Spearman_r)$mean
  data.frame(pathology = paste0(group_name, "_ihc"), baseline = baseline)
})

# ── CA1 retro IHC — 24 mpi ───────────────────────────────────────────────────
# ca1_ret uses HIP injection, retro direction (same path pattern as hipp_ret)
ca1_fp <- paste0("../sir_result_csvs/baselines/Total_Pathology_", ihc_time,
                 "_MPIHIP_injection_rebuilt_rvm_retro.csv")

# Try a dedicated CA1 file if it exists, otherwise fall back
ca1_fp_alt <- paste0("../sir_result_csvs/baselines/Total_Pathology_CA1_", ihc_time,
                     "_MPIHIP_injection_rebuilt_rvm_retro.csv")

if (file.exists(ca1_fp_alt)) {
  ca1_data <- read.csv(ca1_fp_alt, header = FALSE)[, -1]
} else {
  ca1_data <- read.csv(ca1_fp, header = FALSE)[, -1]
}
colnames(ca1_data) <- c("Spearman_r", "v", "spread_rate", "injection_amount", "peak_timestep")
ca1_result <- data.frame(pathology = "ca1_ret_ihc", baseline = mean_ci(ca1_data$Spearman_r)$mean)

# ── ATROPHY (redo40) — 90 dpi = time index 3 ─────────────────────────────────
atr_group_map <- list(
  "hipp_ant_redo40" = list(region = "hipp", pff = "HuPff", direction = "antero"),
  "hipp_ret_redo40" = list(region = "hipp", pff = "HuPff", direction = "retro"),
  "cp_ant_redo40"   = list(region = "cp",   pff = "MsPff", direction = "antero"),
  "cp_ret_redo40"   = list(region = "cp",   pff = "MsPff", direction = "retro")
)

atr_time_idx <- 3  # 90 dpi

build_atr_path <- function(time_idx, region, pff, direction) {
  if (region == "hipp") {
    paste0("../sir_result_csvs/baselines/Time_", time_idx,
           "hipp_rgn_t_stats_hemi", pff,
           "_hemiPBS_rebuilt_rvm_", direction, ".csv")
  } else {
    paste0("../sir_result_csvs/baselines/Time_", time_idx,
           "cp_rgn_t_stats_full_hemi", pff,
           "_hemiPBS_rebuilt_rvm_", direction, ".csv")
  }
}

atr_results <- lapply(names(atr_group_map), function(group_name) {
  m  <- atr_group_map[[group_name]]
  fp <- build_atr_path(atr_time_idx, m$region, m$pff, m$direction)
  data <- read.csv(fp, header = FALSE)[, -1]
  colnames(data) <- c("Spearman_r", "v", "spread_rate", "injection_amount",
                      "k1_atrophy", "k2_atrophy", "peak_timestep")
  baseline <- mean_ci(data$Spearman_r)$mean
  data.frame(pathology = group_name, baseline = baseline)
})

# ── Combine and write ─────────────────────────────────────────────────────────
all_results <- do.call(rbind, c(ihc_results, list(ca1_result), atr_results))
rownames(all_results) <- NULL

print(all_results)

out_dir <- "../tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(all_results, file = file.path(out_dir, "baselines.csv"), row.names = TRUE)
message("Saved: ", file.path(out_dir, "baselines.csv"))
