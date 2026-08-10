library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(Cairo)

OUTPUT_DIR <- "../figures/"

FILE_GROUPS_IHC <- list(
  list(
    directory       = "../optuna_csvs/35/",
    broad_pattern   = "^retFalse_Total_Pathology_24_MPICP_injection_.+_rebuilt_rvm\\.csv$",
    exclude_pattern = "^retFalse_Total_Pathology_24_MPICP_injection_[0-9]+_rebuilt_rvm\\.csv$",
    default_label   = "cp, ant, top 40 clearance",
    overrides       = list(
      "^retFalse_Total_Pathology_24_MPICP_injection_Zmat4_\\d+_rebuilt_rvm\\.csv$" = "cp, Zmat4, ant"
    )
  ),
  list(
    directory       = "../optuna_csvs/35/",
    broad_pattern   = "^retFalse_Total_Pathology_24_MPICP_injection_[0-9]+_rebuilt_rvm\\.csv$",
    exclude_pattern = NULL,
    default_label   = "cp, ant, baseline",
    overrides       = list()
  ),
  list(
    directory       = "../optuna_csvs/24/",
    broad_pattern   = "^retTrue_Total_Pathology_24_MPIHIP_injection_.+_rebuilt_rvm\\.csv$",
    exclude_pattern = "^retTrue_Total_Pathology_24_MPIHIP_injection_[0-9]+_rebuilt_rvm\\.csv$",
    default_label   = "hipp, ret, top 40 clearance",
    overrides       = list(
      "^retTrue_Total_Pathology_24_MPIHIP_injection_Ctnnbip1_\\d+_rebuilt_rvm\\.csv$" = "hipp, Ctnnbip1, ret"
    )
  ),
  list(
    directory       = "../optuna_csvs/24/",
    broad_pattern   = "^retTrue_Total_Pathology_24_MPIHIP_injection_[0-9]+_rebuilt_rvm\\.csv$",
    exclude_pattern = NULL,
    default_label   = "hipp, ret, baseline",
    overrides       = list()
  )
)

FILE_GROUPS_ATROPHY <- list(
  list(
    directory       = "../optuna_csvs/35/",
    broad_pattern   = "^retFalse_Time_3cp_rgn_t_stats_full_hemiMsPff_hemiPBS_.+_rebuilt_rvm\\.csv$",
    exclude_pattern = "^retFalse_Time_3cp_rgn_t_stats_full_hemiMsPff_hemiPBS_[0-9]+_rebuilt_rvm\\.csv$",
    default_label   = "cp, ant, top 40 clearance",
    overrides       = list(
      "^retFalse_Time_3cp_rgn_t_stats_full_hemiMsPff_hemiPBS_Dnajc17_\\d+_rebuilt_rvm\\.csv$" = "cp, Dnajc17, ant"
    )
  ),
  list(
    directory       = "../optuna_csvs/35/",
    broad_pattern   = "^retFalse_Time_3cp_rgn_t_stats_full_hemiMsPff_hemiPBS_[0-9]+_rebuilt_rvm\\.csv$",
    exclude_pattern = NULL,
    default_label   = "cp, ant, baseline",
    overrides       = list()
  ),
  list(
    directory       = "../optuna_csvs/40/",
    broad_pattern   = "^retTrue_Time_3hipp_rgn_t_stats_hemiHuPff_hemiPBS_.+_rebuilt_rvm\\.csv$",
    exclude_pattern = "^retTrue_Time_3hipp_rgn_t_stats_hemiHuPff_hemiPBS_[0-9]+_rebuilt_rvm\\.csv$",
    default_label   = "hipp, ret, top 40 clearance",
    overrides       = list(
      "^retTrue_Time_3hipp_rgn_t_stats_hemiHuPff_hemiPBS_Twist2_\\d+_rebuilt_rvm\\.csv$" = "hipp, Twist2, ret"
    )
  ),
  list(
    directory       = "../optuna_csvs/40/",
    broad_pattern   = "^retTrue_Time_3hipp_rgn_t_stats_hemiHuPff_hemiPBS_[0-9]+_rebuilt_rvm\\.csv$",
    exclude_pattern = NULL,
    default_label   = "hipp, ret, baseline",
    overrides       = list()
  )
)

# ── Color palettes: warm = CP/anterior, cool = hippocampal ────────────────────

CUSTOM_COLORS_IHC <- c(
  "cp, ant, baseline"           = "#E8A000",
  "cp, ant, top 40 clearance"   = "#8B0000",
  "cp, Zmat4, ant"              = "#C44B00",
  "hipp, ret, baseline"         = "#22D4E8",
  "hipp, ret, top 40 clearance" = "#0C3D6B",
  "hipp, Ctnnbip1, ret"         = "#3A84C9"
)

CUSTOM_COLORS_ATROPHY <- c(
  "cp, ant, baseline"           = "#E8A000",
  "cp, Dnajc17, ant"            = "#C44B00",
  "cp, ant, top 40 clearance"   = "#8B0000",
  "hipp, ret, baseline"         = "#22D4E8",
  "hipp, Twist2, ret"           = "#3A84C9",
  "hipp, ret, top 40 clearance" = "#0C3D6B"
)

# ── Progressive reveal palette (same colors, same order) ──────────────────────

CUSTOM_COLORS_ATROPHY_PROGRESSIVE <- CUSTOM_COLORS_ATROPHY

REVEAL_ORDER_ATROPHY <- c(
  "cp, ant, baseline",
  "cp, Dnajc17, ant",
  "cp, ant, top 40 clearance",
  "hipp, ret, baseline",
  "hipp, Twist2, ret",
  "hipp, ret, top 40 clearance"
)

# ── Trim data to the dense x-range for each simulation group ──────────────────

trim_to_dense <- function(data, x_col, trim_frac = 0.05) {
  data %>%
    group_by(simulation) %>%
    filter(
      .data[[x_col]] >= quantile(.data[[x_col]], trim_frac,  na.rm = TRUE),
      .data[[x_col]] <= quantile(.data[[x_col]], 1 - trim_frac, na.rm = TRUE)
    ) %>%
    ungroup()
}

# ── Generic param plot function ────────────────────────────────────────────────

plot_param_curves <- function(combined_data, custom_colors, output_dir,
                              base_name, title, x_col, x_label,
                              ylim = c(0.4, 0.75),
                              x_log = FALSE, xlim = NULL,
                              trim_frac = 0.05,
                              panel_width  = unit(8, "in"),
                              panel_height = unit(5, "in")) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  smooth_data <- trim_to_dense(combined_data, x_col, trim_frac)

  p <- ggplot(mapping = aes(x = .data[[x_col]], y = value,
                             colour = simulation, fill = simulation)) +
    geom_point(
      data  = combined_data,
      alpha = 0.08, size = 1.2, stroke = 0, shape = 16
    ) +
    geom_smooth(
      data     = smooth_data,
      method   = "loess", span = 0.1, se = TRUE,
      linetype = "dashed", alpha = 0.2
    ) +
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values  = custom_colors) +
    ylim(ylim[1], ylim[2]) +
    labs(title = title, x = x_label, y = "Spearman r") +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 26, face = "bold"),
      axis.title.x    = element_text(size = 26, face = "bold"),
      axis.title.y    = element_text(size = 26, face = "bold"),
      legend.title    = element_text(size = 26, face = "bold"),
      legend.text     = element_text(size = 24),
      legend.key.size = unit(1.5, "cm"),
      axis.text.x     = element_text(size = 24),
      axis.text.y     = element_text(size = 24)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 7)))

  if (x_log) {
    p <- p + scale_x_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    )
  } else if (!is.null(xlim)) {
    p <- p + xlim(xlim[1], xlim[2])
  }

  # Fix the panel (plot area) size; total image expands to accommodate legends
  gt <- ggplotGrob(p)
  panels <- which(gt$layout$name == "panel")
  gt$widths[unique(gt$layout$l[panels])]  <- panel_width
  gt$heights[unique(gt$layout$t[panels])] <- panel_height

  png(file.path(output_dir, paste0(base_name, ".png")),
      width = 14, height = 9, units = "in", res = 150)
  grid::grid.draw(gt)
  dev.off()

  invisible(p)
}

# ── Helper: load all file groups ──────────────────────────────────────────────

resolve_label <- function(filename, overrides, default_label) {
  for (pattern in names(overrides)) {
    if (grepl(pattern, filename)) return(overrides[[pattern]])
  }
  default_label
}

load_file_groups <- function(file_groups) {
  data_list <- list()

  for (group in file_groups) {
    csv_files <- list.files(group$directory,
                            pattern    = group$broad_pattern,
                            full.names = TRUE)
    if (!is.null(group$exclude_pattern)) {
      csv_files <- csv_files[!grepl(group$exclude_pattern, basename(csv_files), perl = TRUE)]
    }
    print(paste0("Group:", group$broad_pattern))
    print(paste0("Num files:", length(csv_files)))

    for (file in csv_files) {
      df <- read_csv(file, show_col_types = FALSE)
      df[["simulation"]] <- resolve_label(basename(file),
                                          group$overrides,
                                          group$default_label)
      data_list[[file]] <- df
    }
  }
  bind_rows(data_list, .id = "source")
}

# ── Main ───────────────────────────────────────────────────────────────────────

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── IHC plots ─────────────────────────────────────────────────────────────────

message("Loading IHC file groups …")
combined_ihc <- load_file_groups(FILE_GROUPS_IHC)

plot_param_curves(
  combined_data = combined_ihc,
  custom_colors = CUSTOM_COLORS_IHC,
  output_dir    = OUTPUT_DIR,
  base_name     = "ihc_velocity_curves",
  title         = "24 MPI pSyn Correlation vs aSyn Velocity",
  x_col         = "params_v",
  x_label       = "aSyn Velocity",
  ylim          = c(0.4, 0.75),
  x_log         = TRUE,
  panel_width  = unit(7, "in"),
  panel_height = unit(5, "in")
)

plot_param_curves(
  combined_data = combined_ihc,
  custom_colors = CUSTOM_COLORS_IHC,
  output_dir    = OUTPUT_DIR,
  base_name     = "ihc_spread_rate_curves",
  title         = "24 MPI pSyn Correlation vs aSyn Spreading Rate",
  x_col         = "params_spread_rate",
  x_label       = "aSyn Spreading Rate",
  ylim          = c(0.4, 0.75),
  panel_width  = unit(7, "in"),
  panel_height = unit(5, "in")
)

# ── Atrophy plots ──────────────────────────────────────────────────────────────

message("Loading atrophy file groups …")
combined_atrophy <- load_file_groups(FILE_GROUPS_ATROPHY)

plot_param_curves(
  combined_data = combined_atrophy,
  custom_colors = CUSTOM_COLORS_ATROPHY,
  output_dir    = OUTPUT_DIR,
  base_name     = "atrophy_velocity_curves",
  title         = "Atrophy Correlation vs aSyn Velocity",
  x_col         = "params_v",
  x_label       = "aSyn Velocity",
  ylim          = c(0.1, 0.8),
  x_log         = TRUE,
  panel_width  = unit(7, "in"),
  panel_height = unit(5, "in")
)

plot_param_curves(
  combined_data = combined_atrophy,
  custom_colors = CUSTOM_COLORS_ATROPHY,
  output_dir    = OUTPUT_DIR,
  base_name     = "atrophy_spread_rate_curves",
  title         = "Atrophy Correlation vs aSyn Spreading Rate",
  x_col         = "params_spread_rate",
  x_label       = "aSyn Spreading Rate",
  ylim          = c(0.1, 0.8),
  panel_width  = unit(7, "in"),
  panel_height = unit(5, "in")
)

plot_param_curves(
  combined_data = combined_atrophy,
  custom_colors = CUSTOM_COLORS_ATROPHY,
  output_dir    = OUTPUT_DIR,
  base_name     = "atrophy_k1_curves",
  title         = "Atrophy Correlation vs k1",
  x_col         = "params_k1_atrophy",
  x_label       = "k1",
  ylim          = c(0.1, 0.8),
  xlim          = c(0, 1),
  panel_width  = unit(7, "in"),
  panel_height = unit(5, "in")
)

plot_param_curves(
  combined_data = combined_atrophy,
  custom_colors = CUSTOM_COLORS_ATROPHY,
  output_dir    = OUTPUT_DIR,
  base_name     = "atrophy_k2_curves",
  title         = "Atrophy Correlation vs k2",
  x_col         = "params_k2_atrophy",
  x_label       = "k2",
  ylim          = c(0.1, 0.8),
  xlim          = c(0, 1),
  panel_width  = unit(7, "in"),
  panel_height = unit(5, "in")
)

# ── Progressive reveal: atrophy spread rate ────────────────────────────────────

for (n in seq_along(REVEAL_ORDER_ATROPHY)) {
  active_groups <- REVEAL_ORDER_ATROPHY[1:n]

  subset_data <- combined_atrophy %>%
    filter(simulation %in% active_groups) %>%
    mutate(simulation = factor(simulation, levels = active_groups))

  active_colors <- CUSTOM_COLORS_ATROPHY_PROGRESSIVE[active_groups]

  plot_param_curves(
    combined_data = subset_data,
    custom_colors = active_colors,
    output_dir    = OUTPUT_DIR,
    base_name     = paste0("atrophy_spread_rate_curves_reveal_", n),
    title         = "Atrophy Correlation vs aSyn Spreading Rate",
    x_col         = "params_spread_rate",
    x_label       = "aSyn Spreading Rate",
    ylim          = c(0.1, 0.8),
    panel_width  = unit(7, "in"),
    panel_height = unit(5, "in")
  )

  message(sprintf("Saved progressive plot %d/%d (%s)",
                  n, length(REVEAL_ORDER_ATROPHY),
                  paste(active_groups, collapse = " + ")))
}

message("Done. Figures saved to ", OUTPUT_DIR)