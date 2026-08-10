library(ggplot2)
library(dplyr)

setwd(".")
fig_dir <- "../figures/"
PLOT_WIDTH  <- 6.5
PLOT_HEIGHT <- 7.0
PLOT_DPI    <- 150

mean_ci <- function(x, confidence = 0.95) {
  x <- na.omit(x)
  n <- length(x)
  mean_x <- mean(x)
  std_err <- sd(x) / sqrt(n)
  error_margin <- qt((1 + confidence) / 2, df = n-1) * std_err
  return(data.frame(mean = mean_x, ymin = mean_x - error_margin, ymax = mean_x + error_margin))
}

# ── PATHOLOGY groups ──────────────────────────────────────────────────────────
path_group_map <- list(
  "hipp_ant" = list(injection = "HIP", direction = "antero"),
  "hipp_ret" = list(injection = "HIP", direction = "retro"),
  "cp_ant"   = list(injection = "CP",  direction = "antero"),
  "cp_ret"   = list(injection = "CP",  direction = "retro")
)

path_col_order <- c("Spearman_r", "v", "spread_rate", "injection_amount", "peak_timestep")

path_facet_labels <- c(
  "Spearman_r"       = "Spearman r",
  "v"                = "Velocity (v)",
  "spread_rate"      = "Spread Rate",
  "injection_amount" = "Injection Amount",
  "peak_timestep"    = "Peak Timestep"
)

# Paired groups share axis limits
path_pairs <- list(
  c("hipp_ant", "hipp_ret"),
  c("cp_ant",   "cp_ret")
)

# ── ATROPHY groups ────────────────────────────────────────────────────────────
atr_group_map <- list(
  "hipp_ant_HuPff" = list(region = "hipp", pff = "HuPff", direction = "antero"),
  "hipp_ret_HuPff" = list(region = "hipp", pff = "HuPff", direction = "retro"),
  "cp_ant_HuPff"   = list(region = "cp",   pff = "HuPff", direction = "antero"),
  "cp_ret_HuPff"   = list(region = "cp",   pff = "HuPff", direction = "retro"),
  "cp_ant_MsPff"   = list(region = "cp",   pff = "MsPff", direction = "antero"),
  "cp_ret_MsPff"   = list(region = "cp",   pff = "MsPff", direction = "retro")
)

atr_col_order <- c("Spearman_r", "v", "spread_rate", "k1_atrophy", "k2_atrophy",
                   "injection_amount", "peak_timestep")

atr_facet_labels <- c(
  "Spearman_r"       = "Spearman r",
  "v"                = "Velocity (v)",
  "spread_rate"      = "Spread Rate",
  "k1_atrophy"       = "k1 Atrophy",
  "k2_atrophy"       = "k2 Atrophy",
  "injection_amount" = "Injection Amount",
  "peak_timestep"    = "Peak Timestep"
)

atr_pairs <- list(
  c("hipp_ant_HuPff", "hipp_ret_HuPff"),
  c("cp_ant_HuPff",   "cp_ret_HuPff"),
  c("cp_ant_MsPff",   "cp_ret_MsPff")
)

dpi_labels <- c("1" = "-7 dpi", "2" = "30 dpi", "3" = "90 dpi", "4" = "120 dpi")

# ── Shared helpers ────────────────────────────────────────────────────────────
process_wide <- function(data, time_val, group_name, col_order) {
  combined <- data.frame()
  for (col in col_order) {
    stats <- mean_ci(data[[col]])
    combined <- rbind(combined, data.frame(
      values = data[[col]],
      time   = as.character(time_val),
      group  = group_name,
      column = col,
      mean   = stats$mean,
      ymin   = stats$ymin,
      ymax   = stats$ymax
    ))
  }
  return(combined)
}

# Compute per-column y limits across a list of already-loaded data frames
compute_shared_limits <- function(data_list) {
  combined <- do.call(rbind, data_list)
  combined %>%
    group_by(column) %>%
    summarise(
      lo = min(c(values, ymin), na.rm = TRUE),
      hi = max(c(values, ymax), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    { setNames(mapply(function(lo, hi) c(lo, hi), .$lo, .$hi, SIMPLIFY = FALSE), .$column) }
}

build_plot <- function(all_data, group_name, x_levels, x_labels,
                       col_order, facet_labels, y_limits, x_title = NULL) {
  
  print(group_name)
  for(t in unique(all_data$time)) {
      print(paste0("Time ", as.character(t)))
      print(mean(na.omit(all_data[all_data$column == "Spearman_r" & all_data$time == t, "mean"])))

  }
  all_data$time   <- factor(all_data$time,   levels = x_levels)
  all_data$column <- factor(all_data$column, levels = col_order)
  
  time_colors <- setNames(
    colorRampPalette(c("#a8d1f5", "#08519c"))(length(x_levels)),
    x_levels
  )
  
  # Build a scales lookup so each facet row gets its shared limits
  scales_list <- setNames(
    lapply(col_order, function(col) {
      lims <- y_limits[[col]]
      scale_y_continuous(limits = lims)
    }),
    col_order
  )
  
  ggplot(all_data, aes(x = time, y = values, color = time)) +
    geom_jitter(width = 0.15, alpha = 0.35, size = 1.2) +
    geom_pointrange(
      aes(y = mean, ymin = ymin, ymax = ymax),
      color = "black", linewidth = 0.7, size = 0.5, fatten = 3
    ) +
    facet_grid(column ~ ., scales = "free_y",
               labeller = labeller(column = facet_labels)) +
    scale_color_manual(values = time_colors) +
    scale_x_discrete(labels = x_labels) +
    ggh4x::facetted_pos_scales(y = scales_list) +
    labs(
      title    = group_name,
      subtitle = "Points = individual values; black = mean ± 95% CI",
      x        = x_title,
      y        = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 13, margin = margin(b = 3)),
      plot.subtitle      = element_text(size = 9, color = "grey50", margin = margin(b = 6)),
      strip.text         = element_text(face = "bold", size = 7),
      strip.background   = element_rect(fill = "grey92", color = NA),
      axis.text.x        = element_text(size = 9, angle = 30, hjust = 1),
      axis.text.y        = element_text(size = 8),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.spacing      = unit(0.6, "lines"),
      panel.border       = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      legend.position    = "none"
    )
}

# ── PATHOLOGY pipeline ────────────────────────────────────────────────────────
read_pathology <- function(group_name, times) {
  combined_data <- data.frame()
  m <- path_group_map[[group_name]]
  for (time in times) {
    fp   <- paste0("../sir_result_csvs/baselines/Total_Pathology_", time,
                   "_MPI", m$injection, "_injection_rebuilt_rvm_", m$direction, ".csv")
    data <- read.csv(fp, header = FALSE)[, -1]
    colnames(data) <- c("Spearman_r", "v", "spread_rate", "injection_amount", "peak_timestep")
    combined_data  <- rbind(combined_data, process_wide(data, time, group_name, path_col_order))
  }
  return(combined_data)
}

plot_pathology <- function(pairs, times) {
  mpi_labels <- setNames(paste0(times, " mpi"), times)
  for (pair in pairs) {
    # Load both groups, compute shared limits, then plot each
    data_list  <- lapply(pair, read_pathology, times = times)
    y_limits   <- compute_shared_limits(data_list)
    for (i in seq_along(pair)) {
      g <- pair[i]
      p <- build_plot(data_list[[i]], g,
                      x_levels     = times,
                      x_labels     = mpi_labels,
                      col_order    = path_col_order,
                      facet_labels = path_facet_labels,
                      y_limits     = y_limits,
                      x_title      = "MPI Timepoint")
      ggsave(paste0(fig_dir, g, "_pathology_plot.png"), p,
             width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
      message("Saved pathology plot: ", g)
    }
  }
}

# ── ATROPHY pipeline ──────────────────────────────────────────────────────────
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

read_atrophy <- function(group_name, time_indices) {
  combined_data <- data.frame()
  m <- atr_group_map[[group_name]]
  for (ti in time_indices) {
    fp   <- build_atr_path(ti, m$region, m$pff, m$direction)
    data <- read.csv(fp, header = FALSE)[, -1]
    colnames(data) <- c("Spearman_r", "v", "spread_rate", "injection_amount",
                        "k1_atrophy", "k2_atrophy", "peak_timestep")
    combined_data  <- rbind(combined_data, process_wide(data, ti, group_name, atr_col_order))
  }
  return(combined_data)
}

plot_atrophy <- function(pairs, time_indices) {
  dpi_x_labels <- dpi_labels[as.character(time_indices)]
  for (pair in pairs) {
    data_list <- lapply(pair, read_atrophy, time_indices = time_indices)
    y_limits  <- compute_shared_limits(data_list)
    for (i in seq_along(pair)) {
      g <- pair[i]
      p <- build_plot(data_list[[i]], g,
                      x_levels     = as.character(time_indices),
                      x_labels     = dpi_x_labels,
                      col_order    = atr_col_order,
                      facet_labels = atr_facet_labels,
                      y_limits     = y_limits,
                      x_title      = "Timepoint (dpi)")
      ggsave(paste0(fig_dir, g, "_atrophy_plot.png"), p,
             width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
      message("Saved atrophy plot: ", g)
    }
  }
}

# ── Run ───────────────────────────────────────────────────────────────────────
plot_pathology(
  pairs = path_pairs,
  times = c("0.5", "1", "3", "6", "12", "18", "24")
)

plot_atrophy(
  pairs        = atr_pairs,
  time_indices = c(1, 2, 3, 4)
)