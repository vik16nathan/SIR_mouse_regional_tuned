library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(reticulate)

scipy_stats <- import("scipy.stats")

spearmanr <- function(x, y) {
  scipy_stats$spearmanr(as.numeric(x), as.numeric(y))$statistic
}

# ── Configuration ─────────────────────────────────────────────────────────────

WORKING_DIR        <- "."
ABA_NAMES_FILE     <- "../../derivatives/SIR_inputs/csvs/yohan_source_full.csv"
N_REGIONS          <- 209
OUTPUT_DIR         <- "../figures/"

EMPIRICAL_FILES    <- c(
  "../../preprocessed/steph_janice_regional_jacobians/rgn_t_stats_full_hemiMsPff_hemiPBS.csv",
  "../../preprocessed/steph_janice_regional_jacobians/hipp_rgn_t_stats_hemiHuPff_hemiPBS.csv",
  "../../preprocessed/shady_ihc/CP_injection.csv",
  "../../preprocessed/shady_ihc/HIP_injection.csv"
)
EMPIRICAL_COLS     <- c("Time.3", "Time.3", "Total.Pathology_24.MPI", "Total.Pathology_24.MPI")
CLEARANCE_GENES    <- c("Dnajc17", "Twist2", "Zmat4", "Ctnnbip1")
CLEARANCE_GENE_DIR <- "../../derivatives/yohan_ge_filt/mr10vv0.2/"
CLEARANCE_FILES    <- paste0(CLEARANCE_GENE_DIR, CLEARANCE_GENES, ".csv")
EPIS               <- c("CP", "DG", "CP", "CA1")
PATHOLOGIES        <- c("atrophy", "atrophy", "pSyn", "pSyn")

SIM_INPUT_FILES <- c(
  "../simulations/abm_spread_v.104.87116708452713.spread_rate.0.02122582417807315.dt.0.1.seed.35.injection_amount.97.7903636896988.clearance_gene.Dnajc17.k1.0.03066890459365682.k2.0.42327816662213974..csv",
  "../simulations/retro_abm_spread_v.1.0862247470116264.spread_rate.0.027270129876801334.dt.0.1.seed.40.injection_amount.35.5593507054226.clearance_gene.Twist2.k1.0.7296383313819099.k2.0.059192558911420556..csv",
  "../simulations/abm_spread_v.2.723990344784835.spread_rate.0.05812310477637486.dt.0.1.seed.35.injection_amount.82.23688396384566.clearance_gene.Zmat4...csv",
  "../simulations/abm_spread_v.0.29188965403566514.spread_rate.0.0006938150874085919.dt.0.1.seed.24.injection_amount.54.177718195890705.clearance_gene.Ctnnbip1...csv"
)

MULTIPANEL_NCOL <- 2

# ── Helpers ───────────────────────────────────────────────────────────────────

clearance_gene_col <- function(clearance_file) {
  tools::file_path_sans_ext(basename(clearance_file))
}

load_aba_names <- function(aba_names_file, n_regions) {
  aba              <- read.csv(aba_names_file, header = TRUE)
  aba              <- rbind(aba, aba)
  aba              <- aba[, -1]
  alternating_vals <- rep(c("right", "left"), each = n_regions)
  list(aba_region_names = aba, alternating_vals = alternating_vals)
}

load_panel_inputs <- function(empirical_file, empirical_col,
                              clearance_file, aba, alternating_vals) {
  results_raw       <- read.csv(empirical_file)
  result_rownames   <- results_raw[, 1]
  empirical_col     <- make.names(empirical_col)
  results           <- as.data.frame(results_raw[, empirical_col])
  rownames(results) <- result_rownames

  gene_col                <- clearance_gene_col(clearance_file)
  clearance_df            <- read.csv(clearance_file)
  clearance_df            <- rbind(clearance_df, clearance_df)
  clearance_df$hemisphere <- alternating_vals
  rownames(clearance_df)  <- paste(alternating_vals, aba, sep = " ")
  clearance_df            <- clearance_df[, gene_col, drop = FALSE]

  list(results = results, clearance_df = clearance_df, gene_col = gene_col)
}

prepare_panel_data <- function(sim_file, panel_inputs, aba, alternating_vals) {
  sim_data            <- read.csv(sim_file, header = FALSE)
  sim_data            <- cbind(aba, sim_data)
  sim_data$hemisphere <- alternating_vals
  rownames(sim_data)  <- paste(sim_data[, ncol(sim_data)],
                               sim_data[, 1], sep = " ")
  sim_data <- sim_data[, -1]
  sim_data <- sim_data[, -ncol(sim_data)]

  common_rows <- sort(Reduce(intersect, list(
    rownames(panel_inputs$results),
    rownames(sim_data),
    rownames(panel_inputs$clearance_df)
  )))

  empirical_common <- panel_inputs$results[common_rows, ]
  simulated_common <- sim_data[common_rows, ]

  correlations <- vapply(
    seq_len(ncol(simulated_common)),
    function(i) spearmanr(empirical_common, simulated_common[, i]),
    numeric(1)
  )

  peak_timestep  <- which.max(correlations)
  simulated_peak <- simulated_common[, peak_timestep]
  sim_v_emp      <- data.frame(simulated = simulated_peak,
                               empirical = empirical_common)
  rownames(sim_v_emp) <- common_rows

  gene_col   <- panel_inputs$gene_col
  df         <- sim_v_emp[common_rows, ]
  clr_common <- panel_inputs$clearance_df[common_rows, , drop = FALSE]

  df[[gene_col]]      <- clr_common[[gene_col]]
  df$scaled_simulated <- as.numeric(scale(df$simulated))
  df$scaled_empirical <- as.numeric(scale(df$empirical))
  df$scaled_gene      <- as.numeric(scale(df[[gene_col]]))

  n_labels        <- 3
  df$residual     <- abs(df$scaled_simulated - df$scaled_empirical)
  resid_threshold <- sort(df$residual, decreasing = TRUE)[n_labels]
  df$highlight    <- ifelse(df$residual >= resid_threshold, rownames(df), NA)

  list(df = df, clr_common = clr_common, gene_col = gene_col,
       peak_corr = max(correlations))
}

# ── Single-panel plot builders ────────────────────────────────────────────────

make_sim_vs_emp_panel <- function(df, gene_col, epi, pathology, clearance_gene) {
  z_min      <- min(df$scaled_gene)
  z_max      <- max(df$scaled_gene)
  scaled_mid <- abs(z_min) / (abs(z_min) + z_max)
  axis_label <- if (pathology == "atrophy") "Atrophy" else "pSyn"

  ggplot(data = df,
         aes(y = scaled_simulated, x = scaled_empirical, fill = scaled_gene)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = c("blue", "white", "red"),
      values  = c(0, scaled_mid, 1),
      limits  = c(z_min, z_max),
      breaks  = c(z_min, 0, z_max),
      labels  = round(c(z_min, 0, z_max), 2),
      oob     = scales::squish
    ) +
    theme_minimal() +
    geom_point(size = 3, shape = 21, stroke = 1, color = "black") +
    geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2", fill = NA) +
    geom_point(
      data = subset(df, !is.na(highlight)),
      aes(y = scaled_simulated, x = scaled_empirical, fill = scaled_gene),
      size = 4, shape = 21, stroke = 2, color = "black"
    ) +
    geom_text_repel(
      data               = df,
      aes(x = scaled_empirical, y = scaled_simulated, label = highlight),
      size               = 5,
      color              = "black",
      max.overlaps       = Inf,
      force              = 15,
      force_pull         = 0.5,
      box.padding        = 1.2,
      point.padding      = 0.8,
      min.segment.length = 0,
      seed               = 42
    ) +
    labs(
      x     = paste("Empirical", axis_label),
      y     = paste("Simulated", axis_label),
      title = paste0(epi, " \u2013 ", clearance_gene),
      fill  = clearance_gene
    ) +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 24, face = "bold"),
      axis.title      = element_text(size = 22, face = "bold"),
      legend.title    = element_text(size = 20, face = "bold"),
      legend.text     = element_text(size = 20),
      legend.key.size = unit(1, "cm"),
      axis.text       = element_text(size = 20)
    )
}

make_gene_vs_sim_panel <- function(df, clr_common, gene_col, epi, pathology,
                                   clearance_gene) {
  corr       <- spearmanr(clr_common[[gene_col]], df[, "simulated"])
  axis_label <- if (pathology == "atrophy") "Atrophy" else "pSyn"

  plot_df <- data.frame(
    x = clr_common[[gene_col]],
    y = df[, "simulated"]
  )

  ggplot(data = plot_df, aes(x = x, y = y)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
    geom_text(
      x     = min(plot_df$x),
      y     = max(plot_df$y),
      label = paste0("Spearman r = ", round(corr, 3)),
      hjust = 0, vjust = 1, color = "black", size = 7
    ) +
    labs(
      x     = paste(gene_col, "expression"),
      y     = paste("Simulated", axis_label),
      title = paste0(epi, " \u2013 ", clearance_gene)
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title      = element_text(hjust = 0.5, size = 24, face = "bold"),
      axis.title      = element_text(size = 22, face = "bold"),
      axis.text       = element_text(size = 20)
    )
}

make_gene_vs_emp_panel <- function(df, clr_common, gene_col, epi, pathology,
                                   clearance_gene) {
  corr       <- spearmanr(clr_common[[gene_col]], df[, "empirical"])
  axis_label <- if (pathology == "atrophy") "Atrophy" else "pSyn"

  plot_df <- data.frame(
    x = clr_common[[gene_col]],
    y = df[, "empirical"]
  )

  ggplot(data = plot_df, aes(x = x, y = y)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
    geom_text(
      x     = min(plot_df$x),
      y     = max(plot_df$y),
      label = paste0("Spearman r = ", round(corr, 3)),
      hjust = 0, vjust = 1, color = "black", size = 7
    ) +
    labs(
      x     = paste(gene_col, "expression"),
      y     = paste("Empirical", axis_label),
      title = paste0(epi, " \u2013 ", clearance_gene)
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title      = element_text(hjust = 0.5, size = 24, face = "bold"),
      axis.title      = element_text(size = 22, face = "bold"),
      axis.text       = element_text(size = 20)
    )
}

# ── Multi-panel assembler ─────────────────────────────────────────────────────

build_multipanel_plot <- function(sim_files, empirical_files, empirical_cols,
                                  clearance_files, clearance_genes,
                                  epis, pathologies, aba_shared,
                                  plot_type = "sim_vs_emp", ncol = 2) {
  n_panels <- length(sim_files)
  stopifnot(
    length(empirical_files) == n_panels,
    length(empirical_cols)  == n_panels,
    length(clearance_files) == n_panels,
    length(clearance_genes) == n_panels,
    length(epis)            == n_panels,
    length(pathologies)     == n_panels
  )

  panels <- vector("list", n_panels)

  for (i in seq_len(n_panels)) {
    message(sprintf("  Panel %d/%d: %s – %s", i, n_panels, epis[i], clearance_genes[i]))

    panel_inputs <- load_panel_inputs(
      empirical_file   = empirical_files[i],
      empirical_col    = empirical_cols[i],
      clearance_file   = clearance_files[i],
      aba              = aba_shared$aba_region_names,
      alternating_vals = aba_shared$alternating_vals
    )

    pd <- prepare_panel_data(
      sim_file         = sim_files[i],
      panel_inputs     = panel_inputs,
      aba              = aba_shared$aba_region_names,
      alternating_vals = aba_shared$alternating_vals
    )

    message(sprintf("    Peak Spearman r = %.3f", pd$peak_corr))

    p <- if (plot_type == "sim_vs_emp") {
      make_sim_vs_emp_panel(
        df             = pd$df,
        gene_col       = pd$gene_col,
        epi            = epis[i],
        pathology      = pathologies[i],
        clearance_gene = clearance_genes[i]
      )
    } else if (plot_type == "gene_vs_sim") {
      make_gene_vs_sim_panel(
        df             = pd$df,
        clr_common     = pd$clr_common,
        gene_col       = pd$gene_col,
        epi            = epis[i],
        pathology      = pathologies[i],
        clearance_gene = clearance_genes[i]
      )
    } else {
      make_gene_vs_emp_panel(
        df             = pd$df,
        clr_common     = pd$clr_common,
        gene_col       = pd$gene_col,
        epi            = epis[i],
        pathology      = pathologies[i],
        clearance_gene = clearance_genes[i]
      )
    }

    panels[[i]] <- p
  }

  wrap_plots(panels, ncol = ncol) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 24, face = "bold"))
}

# ── Main ──────────────────────────────────────────────────────────────────────

setwd(WORKING_DIR)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

aba_shared <- load_aba_names(ABA_NAMES_FILE, N_REGIONS)

message("Building sim-vs-empirical multipanel …")
fig1 <- build_multipanel_plot(
  sim_files       = SIM_INPUT_FILES,
  empirical_files = EMPIRICAL_FILES,
  empirical_cols  = EMPIRICAL_COLS,
  clearance_files = CLEARANCE_FILES,
  clearance_genes = CLEARANCE_GENES,
  epis            = EPIS,
  pathologies     = PATHOLOGIES,
  aba_shared      = aba_shared,
  plot_type       = "sim_vs_emp",
  ncol            = MULTIPANEL_NCOL
)
ggsave(file.path(OUTPUT_DIR, "multipanel_sim_vs_emp.png"),
       plot = fig1, width = 16, height = 14, dpi = 500)

message("Building gene-vs-simulated multipanel …")
fig2 <- build_multipanel_plot(
  sim_files       = SIM_INPUT_FILES,
  empirical_files = EMPIRICAL_FILES,
  empirical_cols  = EMPIRICAL_COLS,
  clearance_files = CLEARANCE_FILES,
  clearance_genes = CLEARANCE_GENES,
  epis            = EPIS,
  pathologies     = PATHOLOGIES,
  aba_shared      = aba_shared,
  plot_type       = "gene_vs_sim",
  ncol            = MULTIPANEL_NCOL
)
ggsave(file.path(OUTPUT_DIR, "multipanel_gene_vs_sim.png"),
       plot = fig2, width = 16, height = 14)

message("Building gene-vs-empirical multipanel …")
fig3 <- build_multipanel_plot(
  sim_files       = SIM_INPUT_FILES,
  empirical_files = EMPIRICAL_FILES,
  empirical_cols  = EMPIRICAL_COLS,
  clearance_files = CLEARANCE_FILES,
  clearance_genes = CLEARANCE_GENES,
  epis            = EPIS,
  pathologies     = PATHOLOGIES,
  aba_shared      = aba_shared,
  plot_type       = "gene_vs_emp",
  ncol            = MULTIPANEL_NCOL
)
ggsave(file.path(OUTPUT_DIR, "multipanel_gene_vs_emp.png"),
       plot = fig3, width = 16, height = 14)

message("Done. Figures saved to ", OUTPUT_DIR)