library("pacman")
library("reticulate")
pacman::p_load(dplyr, readr, ggplot2, scales, patchwork, purrr)

pd <- import("pandas")

# ── File paths ─────────────────────────────────────────────────────────────────
ABA_NAMES_FILE <- "../../derivatives/SIR_inputs/csvs/yohan_source_full.csv"
BASELINES_FILE <- "../tables/baselines.csv"
N_REGIONS      <- 209

hipp_retro_sim_file <- "../simulations/retro_abm_spread_v.1.0862247470116264.spread_rate.0.027270129876801334.dt.0.1.seed.40.injection_amount.35.5593507054226.clearance_gene.Twist2.k1.0.7296383313819099.k2.0.059192558911420556..csv"
retro_param_file    <- "../../derivatives/SIR_inputs/params_regionalized_voxel_weights_after_qc_0413_64pct_retro.pkl"

cp_antero_sim_file  <- "../simulations/abm_spread_v.104.87116708452713.spread_rate.0.02122582417807315.dt.0.1.seed.35.injection_amount.97.7903636896988.clearance_gene.Dnajc17.k1.0.03066890459365682.k2.0.42327816662213974..csv"
antero_param_file   <- "../../derivatives/SIR_inputs/params_regionalized_voxel_weights_after_qc_0413_64pct.pkl"

cp_psyn_sim_file  <- "../simulations/abm_spread_v.2.723990344784835.spread_rate.0.05812310477637486.dt.0.1.seed.35.injection_amount.82.23688396384566.clearance_gene.Zmat4...csv"
ca1_psyn_sim_file <- "../simulations/abm_spread_v.0.29188965403566514.spread_rate.0.0006938150874085919.dt.0.1.seed.24.injection_amount.54.177718195890705.clearance_gene.Ctnnbip1...csv"

tstat_hipp_file <- "../../preprocessed/steph_janice_regional_jacobians/hipp_rgn_t_stats_hemiHuPff_hemiPBS.csv"
tstat_cp_file   <- "../../preprocessed/steph_janice_regional_jacobians/rgn_t_stats_full_hemiMsPff_hemiPBS.csv"
TSTAT_COL       <- "Time.3"

psyn_cp_file <- "../../preprocessed/shady_ihc/CP_injection.csv"
psyn_dg_file <- "../../preprocessed/shady_ihc/HIP_injection.csv"
PSYN_COL     <- "Total.Pathology_24.MPI"

ge_dir    <- "../../derivatives/yohan_ge_filt/mr10vv0.2/"
snca_file <- paste0(ge_dir, "Snca.csv")

cl_gene_files <- list(
  cp_atrophy = paste0(ge_dir, "Dnajc17.csv"),
  cp_psyn    = paste0(ge_dir, "Zmat4.csv"),
  dg_atrophy = paste0(ge_dir, "Twist2.csv"),
  ca1_psyn   = paste0(ge_dir, "Ctnnbip1.csv")
)

BASELINE_KEYS <- list(
  cp_atrophy = "cp_ant_redo40_MsPff",
  cp_psyn    = "cp_ant_ihc",
  dg_atrophy = "hipp_ret_redo40",
  ca1_psyn   = "ca1_ret_ihc"
)

# ── Colours ────────────────────────────────────────────────────────────────────
HEMI_COLS  <- c("Right" = "#E8604C", "Left" = "#4CA8E8")
# Epicentre is always right → ipsilateral = right, contralateral = left
SCOPE_COLS <- c("Bilateral"     = "#555555",
                "Ipsilateral"   = "#E8604C",
                "Contralateral" = "#4CA8E8")

# ── Shared themes ──────────────────────────────────────────────────────────────
theme_corr <- function(base_size = 15) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90"),
      plot.title       = element_text(face = "bold", hjust = 0.5, size = base_size + 1),
      plot.subtitle    = element_text(hjust = 0.5, colour = "grey40", size = base_size - 1),
      axis.title       = element_text(face = "bold", size = base_size - 1),
      axis.text        = element_text(colour = "grey30", size = base_size - 2),
      legend.position  = "bottom",
      legend.title     = element_text(size = base_size - 1),
      legend.text      = element_text(size = base_size - 2)
    )
}

theme_bar <- function(base_size = 15) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90"),
      plot.title        = element_text(face = "bold", hjust = 0.5, size = base_size + 1),
      plot.subtitle     = element_text(hjust = 0.5, colour = "grey40", size = base_size - 1),
      axis.title        = element_text(face = "bold", size = base_size - 1),
      axis.text         = element_text(colour = "grey30", size = base_size - 2),
      legend.position   = "bottom",
      legend.title      = element_text(size = base_size - 1),
      legend.text       = element_text(size = base_size - 2)
    )
}

spearman <- function(x, y) cor(as.numeric(x), as.numeric(y), method = "spearman")

# ── Data loaders ───────────────────────────────────────────────────────────────
load_aba_names <- function(aba_names_file, n_regions) {
  aba              <- read.csv(aba_names_file, header = TRUE)
  aba              <- rbind(aba, aba)
  aba              <- aba[, -1]
  alternating_vals <- rep(c("right", "left"), each = n_regions)
  list(aba_region_names = aba, alternating_vals = alternating_vals)
}

load_empirical <- function(csv_file, col) {
  df           <- read.csv(csv_file)
  rownames(df) <- df[, 1]
  col          <- make.names(col)
  as.data.frame(df[, col, drop = FALSE])
}

load_snca <- function(snca_csv, aba_names, alternating_vals) {
  df            <- read.csv(snca_csv)
  col           <- if ("Snca" %in% names(df)) "Snca" else names(df)[ncol(df)]
  df            <- rbind(df, df)
  df$hemisphere <- alternating_vals
  rownames(df)  <- paste(alternating_vals, aba_names, sep = " ")
  df[, col, drop = FALSE]
}

load_params <- function(param_pkl, aba_names, alternating_vals) {
  params <- pd$read_pickle(param_pkl)
  message("Params keys: ", paste(names(params), collapse = ", "))
  rn <- paste(alternating_vals, aba_names, sep = " ")

  conn_strength_df           <- as.data.frame(params$weights)
  colnames(conn_strength_df) <- rn

  distance_df           <- as.data.frame(params$distance)
  colnames(distance_df) <- rn

  region_size_df           <- as.data.frame(unlist(c(params$region_size, params$region_size)))
  rownames(region_size_df) <- rn
  colnames(region_size_df) <- "region_size"

  list(conn_strength = conn_strength_df,
       distance      = distance_df,
       region_size   = region_size_df)
}

load_sim <- function(sim_csv, aba_names, alternating_vals) {
  sim_data            <- read.csv(sim_csv, header = FALSE)
  sim_data            <- cbind(aba_names, sim_data)
  sim_data$hemisphere <- alternating_vals
  rownames(sim_data)  <- paste(sim_data[, ncol(sim_data)], sim_data[, 1], sep = " ")
  sim_data <- sim_data[, -1]
  sim_data <- sim_data[, -ncol(sim_data)]
  sim_data
}

find_peak_timestep <- function(sim_df, empirical_df) {
  common_rows      <- sort(intersect(rownames(empirical_df), rownames(sim_df)))
  print(length(common_rows))
  empirical_common <- empirical_df[common_rows, ]
  simulated_common <- sim_df[common_rows, ]
  cors <- vapply(seq_len(ncol(simulated_common)),
                 function(i) spearman(simulated_common[, i], empirical_common),
                 numeric(1))
  peak_t <- which.max(cors)
  list(peak_t           = peak_t,
       rho              = max(cors, na.rm = TRUE),
       common_rows      = common_rows,
       empirical_common = empirical_common,
       sim_peak         = simulated_common[, peak_t],
       all_cors         = cors)
}

# ── Save helper ────────────────────────────────────────────────────────────────
save_png <- function(p, filename, width = 6, height = 5, dpi = 300) {
  ggsave(file.path("../figures", filename), p, width = width, height = height, dpi = dpi)
  message("Saved: ../figures/", filename)
}

# ── Core scatter ───────────────────────────────────────────────────────────────
make_scatter <- function(x_vec, y_vec,
                         x_lab,
                         y_lab        = "Empirical atrophy (t-statistic)",
                         title_str    = NULL,
                         subtitle_str = NULL,
                         hemi_mask    = NULL,
                         log_x        = FALSE,
                         base_size    = 15) {
  if (log_x) {
    x_vec <- log10(x_vec + .Machine$double.eps)
    x_lab <- paste0(x_lab, " (log10)")
  }

  rho_all <- spearman(x_vec, y_vec)

  if (!is.null(hemi_mask)) {
    rho_R <- spearman(x_vec[hemi_mask],  y_vec[hemi_mask])
    rho_L <- spearman(x_vec[!hemi_mask], y_vec[!hemi_mask])
    annot <- sprintf("\u03c1 (R) = %.3f\n\u03c1 (L) = %.3f\n\u03c1 (all) = %.3f",
                     rho_R, rho_L, rho_all)
    df <- data.frame(x    = x_vec, y = y_vec,
                     hemi = factor(ifelse(hemi_mask, "Right", "Left"),
                                   levels = c("Right", "Left")))
    p <- ggplot(df, aes(x = x, y = y, colour = hemi)) +
      geom_point(alpha = 0.45, size = 1.2) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
      scale_colour_manual(name = "Hemisphere", values = HEMI_COLS,
                          guide = guide_legend(override.aes = list(alpha = 1, size = 2.5)))
  } else {
    annot <- sprintf("\u03c1 = %.3f", rho_all)
    df <- data.frame(x = x_vec, y = y_vec)
    p <- ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.45, size = 1.2, colour = "grey40") +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, colour = "steelblue")
  }

  p +
    annotate("text", x = -Inf, y = Inf,
             hjust = -0.08, vjust = 1.4, label = annot,
             size = 4, fontface = "bold", colour = "grey20", lineheight = 1.3) +
    labs(title = title_str, subtitle = subtitle_str, x = x_lab, y = y_lab) +
    theme_corr(base_size = base_size)
}

# ── Correlation vs. timestep ───────────────────────────────────────────────────
make_cor_timestep_plot <- function(all_cors, peak_t, title_str) {
  df <- data.frame(timestep = seq_along(all_cors), correlation = all_cors)
  ggplot(df, aes(x = timestep, y = correlation)) +
    geom_line(colour = "steelblue", linewidth = 0.7) +
    geom_vline(xintercept = peak_t, linetype = "dashed", colour = "firebrick", linewidth = 0.6) +
    annotate("text", x = peak_t, y = max(all_cors, na.rm = TRUE),
             label = sprintf("peak t = %d\n\u03c1 = %.3f", peak_t, all_cors[peak_t]),
             hjust = -0.1, vjust = 1.1, size = 4.2, fontface = "bold", colour = "firebrick") +
    labs(title = title_str, x = "Timestep", y = "Spearman \u03c1") +
    theme_corr()
}

make_cor_timestep_plot_schematic <- function(all_cors, peak_t, title_str) {
  df <- data.frame(timestep = seq_along(all_cors), correlation = all_cors)
  ggplot(df, aes(x = timestep, y = correlation)) +
    geom_line(colour = "steelblue", linewidth = 0.7) +
    geom_vline(xintercept = peak_t, linetype = "dashed", colour = "firebrick", linewidth = 0.6) +
    labs(title = title_str, x = "Timestep", y = "Spearman \u03c1") +
    theme_corr(base_size = 24)
}

# ── 4-panel input scatter ──────────────────────────────────────────────────────
make_input_4panel <- function(emp_vec, conn_vec, dist_vec, size_vec, snca_vec,
                              hemi_mask, y_lab, panel_title) {
  p1 <- make_scatter(conn_vec, emp_vec, x_lab = "Connection strength", y_lab = y_lab,
                     hemi_mask = hemi_mask, log_x = TRUE)
  p2 <- make_scatter(dist_vec, emp_vec, x_lab = "Distance",            y_lab = y_lab,
                     hemi_mask = hemi_mask)
  p3 <- make_scatter(size_vec, emp_vec, x_lab = "Region size",         y_lab = y_lab,
                     hemi_mask = hemi_mask)
  p4 <- make_scatter(snca_vec, emp_vec, x_lab = "Snca expression",     y_lab = y_lab,
                     hemi_mask = hemi_mask)

  (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = panel_title,
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 17))
    ) &
    theme(legend.position = "none")
}

# ══════════════════════════════════════════════════════════════════════════════
# THREE-LEVEL BARPLOT (bilateral / ipsilateral / contralateral)
# Epicentre is always right hemisphere → ipsilateral = right, contralateral = left
# ══════════════════════════════════════════════════════════════════════════════

# Compute rho for all three scopes for a single predictor vector.
# hemi_mask: logical, TRUE = right (ipsilateral), FALSE = left (contralateral)
# log_transform: apply log10 before correlating
compute_rhos_3scope <- function(pred_vec, emp_vec, hemi_mask, log_transform = FALSE) {
  if (log_transform) pred_vec <- log10(pred_vec + .Machine$double.eps)
  list(
    bilateral     = spearman(pred_vec,               emp_vec),
    ipsilateral   = spearman(pred_vec[hemi_mask],    emp_vec[hemi_mask]),
    contralateral = spearman(pred_vec[!hemi_mask],   emp_vec[!hemi_mask])
  )
}

# Build the long data frame of rhos for all predictors × 3 scopes,
# then draw a horizontal dodged barplot (matching the voxel-level style).
make_rho_barplot_3level <- function(emp_vec, conn_vec, dist_vec, size_vec,
                                    snca_vec, cl_gene_vec, cl_gene_name,
                                    sim_peak_vec, hemi_mask,
                                    panel_title) {

  # Named list of predictors; log = TRUE triggers log10 transform
  predictors <- list(
    list(vec = conn_vec,      label = "Connection\nstrength (log10)", log = TRUE),
    list(vec = dist_vec,      label = "Distance",                     log = FALSE),
    list(vec = size_vec,      label = "Region size",                  log = FALSE),
    list(vec = snca_vec,      label = "Snca",                         log = FALSE),
    list(vec = cl_gene_vec,   label = cl_gene_name,                   log = FALSE),
    list(vec = sim_peak_vec,  label = "Simulated\n(peak)",            log = FALSE)
  )

  rows <- purrr::imap_dfr(predictors, function(p, i) {
    rhos <- compute_rhos_3scope(p$vec, emp_vec, hemi_mask, log_transform = p$log)
    data.frame(
      predictor = p$label,
      scope     = c("Bilateral", "Ipsilateral", "Contralateral"),
      rho       = c(rhos$bilateral, rhos$ipsilateral, rhos$contralateral),
      stringsAsFactors = FALSE
    )
  })


  all_rows <- rows

  # Factor ordering: predictors in display order (top = first predictor), baseline last
  pred_levels <- c(sapply(predictors, `[[`, "label"))
  all_rows$predictor <- factor(all_rows$predictor, levels = rev(pred_levels))
  all_rows$scope     <- factor(all_rows$scope,
                               levels = c("Bilateral", "Ipsilateral", "Contralateral"))

  n_pred <- length(pred_levels)

  ggplot(all_rows, aes(x = rho, y = predictor, fill = scope)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey40") +
    scale_fill_manual(name = "Scope", values = SCOPE_COLS) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2),
                       labels = function(x) sprintf("%.1f", x)) +
    labs(title    = panel_title,
         x        = "Spearman \u03c1",
         y        = NULL,
         caption  = "Ipsilateral = right hemisphere (injection side); Contralateral = left") +
    theme_bar() +
    theme(plot.caption = element_text(size = 10, colour = "grey50"))
}

# ── Full per-condition pipeline ────────────────────────────────────────────────
run_condition <- function(label, display_title, sim_df, emp_df, params,
                          snca_df, epi_index, y_lab, baseline_val,
                          cl_gene_df, cl_gene_name,
                          is_psyn = FALSE) {

  pk           <- find_peak_timestep(sim_df, emp_df)
  common_rows  <- pk$common_rows
  emp_vec      <- as.numeric(pk$empirical_common)
  sim_peak_vec <- as.numeric(pk$sim_peak)

  message(sprintf("%s — peak t = %d, rho = %.3f", label, pk$peak_t, pk$rho))

  snca_vec     <- as.numeric(snca_df[common_rows, , drop = FALSE][, 1])
  cl_gene_vec  <- as.numeric(cl_gene_df[common_rows, , drop = FALSE][, 1])
  conn_vec     <- as.numeric(params$conn_strength[epi_index, common_rows])
  dist_vec     <- as.numeric(params$distance[epi_index, common_rows])
  size_vec     <- as.numeric(params$region_size[common_rows, "region_size"])

  # row names are "right <region>" / "left <region>"; right = ipsilateral (injection side)
  hemi_mask <- grepl("^right", common_rows, ignore.case = TRUE)

  x_sim_lab  <- if (is_psyn) "Simulated I population (peak)" else "Simulated spread (peak)"
  out_suffix <- if (is_psyn) "_psyn" else ""

  # ── Individual scatter plots ────────────────────────────────────────────────
  save_png(make_scatter(snca_vec,     emp_vec, x_lab = "Snca expression",    y_lab = y_lab,
                        hemi_mask = hemi_mask,
                        title_str = paste(display_title, "vs. Snca")),
           paste0(label, out_suffix, "_vs_snca.png"))

  save_png(make_scatter(conn_vec,     emp_vec, x_lab = "Connection strength", y_lab = y_lab,
                        hemi_mask = hemi_mask, log_x = TRUE,
                        title_str = paste(display_title, "vs. connection strength")),
           paste0(label, out_suffix, "_vs_conn.png"))

  save_png(make_scatter(dist_vec,     emp_vec, x_lab = "Distance",            y_lab = y_lab,
                        hemi_mask = hemi_mask,
                        title_str = paste(display_title, "vs. distance")),
           paste0(label, out_suffix, "_vs_dist.png"))

  save_png(make_scatter(size_vec,     emp_vec, x_lab = "Region size",         y_lab = y_lab,
                        hemi_mask = hemi_mask,
                        title_str = paste(display_title, "vs. region size")),
           paste0(label, out_suffix, "_vs_size.png"))

  save_png(make_scatter(sim_peak_vec, emp_vec, x_lab = x_sim_lab,             y_lab = y_lab,
                        hemi_mask = hemi_mask,
                        subtitle_str = sprintf("Peak timestep: %d", pk$peak_t),
                        title_str = paste(display_title, "vs. simulated")),
           paste0(label, out_suffix, "_vs_sim.png"))

  save_png(make_cor_timestep_plot(pk$all_cors, pk$peak_t,
                                  title_str = paste(display_title, "| \u03c1 vs. timestep")),
           paste0(label, out_suffix, "_cor_vs_timestep.png"))

  if (label == "hipp_retro") {
    save_png(make_cor_timestep_plot_schematic(pk$all_cors, pk$peak_t,
                                             title_str = paste(display_title, "| \u03c1 vs. timestep")),
             paste0(label, out_suffix, "_cor_vs_timestep_schematic.png"))
  }

  # ── 4-panel input scatter ──────────────────────────────────────────────────
  p4 <- make_input_4panel(emp_vec, conn_vec, dist_vec, size_vec, snca_vec,
                          hemi_mask, y_lab = y_lab,
                          panel_title = paste(display_title, "| Input predictors"))
  save_png(p4, paste0(label, out_suffix, "_input_4panel.png"), width = 10, height = 8)

  # ── Three-level bilateral/ipsi/contra barplot ──────────────────────────────
  p_bar3 <- make_rho_barplot_3level(
    emp_vec      = emp_vec,
    conn_vec     = conn_vec,
    dist_vec     = dist_vec,
    size_vec     = size_vec,
    snca_vec     = snca_vec,
    cl_gene_vec  = cl_gene_vec,
    cl_gene_name = cl_gene_name,
    sim_peak_vec = sim_peak_vec,
    hemi_mask    = hemi_mask,
    panel_title  = paste(display_title, "| Spearman \u03c1 by hemisphere")
  )
  save_png(p_bar3, paste0(label, out_suffix, "_rho_barplot_3level.png"),
           width = 8, height = 3.5 + 0.45 * 7)   # 7 predictors incl. baseline

  invisible(list(pk = pk, p_bar3 = p_bar3))
}

# ── Cross-modal ────────────────────────────────────────────────────────────────
make_crossmodal_plot <- function(psyn_df, tstat_df, label) {
  common_rows  <- sort(intersect(rownames(psyn_df), rownames(tstat_df)))
  psyn_vec     <- as.numeric(psyn_df[common_rows, ])
  atrophy_vec  <- as.numeric(tstat_df[common_rows, ])
  hemi_mask    <- grepl("^right", common_rows, ignore.case = TRUE)
  title_prefix <- dplyr::case_when(
    grepl("dg|hipp", label, ignore.case = TRUE) ~ "DG",
    grepl("cp",      label, ignore.case = TRUE) ~ "CP",
    TRUE ~ toupper(label)
  )
  save_png(
    make_scatter(psyn_vec, atrophy_vec,
                 x_lab     = "Empirical pSyn",
                 y_lab     = "Empirical atrophy (t-statistic)",
                 title_str = paste(title_prefix, "| Empirical atrophy vs. empirical pSyn"),
                 hemi_mask = hemi_mask),
    paste0(label, "_atrophy_vs_psyn.png")
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════
dir.create("../figures", showWarnings = FALSE)

message("Loading ABA region names …")
aba_shared       <- load_aba_names(ABA_NAMES_FILE, N_REGIONS)
aba_names        <- aba_shared$aba_region_names
alternating_vals <- aba_shared$alternating_vals

cp_epi_index  <- intersect(which(alternating_vals == "right"), which(aba_names == "Caudoputamen"))
dg_epi_index  <- intersect(which(alternating_vals == "right"), which(aba_names == "Dentate gyrus"))
ca1_epi_index <- intersect(which(alternating_vals == "right"), which(aba_names == "Field CA1"))

message("Loading baselines …")
baselines    <- read.csv(BASELINES_FILE, row.names = 1)
get_baseline <- function(key) baselines[baselines$pathology == key, "baseline"]

message("Loading Snca …")
snca_df <- load_snca(snca_file, aba_names, alternating_vals)

message("Loading clearance genes …")
cl_gene_dfs <- lapply(cl_gene_files, function(f) load_snca(f, aba_names, alternating_vals))

message("Loading empirical atrophy t-stat maps …")
tstat_hipp <- load_empirical(tstat_hipp_file, TSTAT_COL)
tstat_cp   <- load_empirical(tstat_cp_file,   TSTAT_COL)

message("Loading empirical pSyn maps …")
psyn_cp <- load_empirical(psyn_cp_file, PSYN_COL)
psyn_dg <- load_empirical(psyn_dg_file, PSYN_COL)

message("Loading params …")
retro_params  <- load_params(retro_param_file,  aba_names, alternating_vals)
antero_params <- load_params(antero_param_file, aba_names, alternating_vals)

message("Loading sim files …")
retro_sim    <- load_sim(hipp_retro_sim_file, aba_names, alternating_vals)
antero_sim   <- load_sim(cp_antero_sim_file,  aba_names, alternating_vals)
cp_psyn_sim  <- load_sim(cp_psyn_sim_file,    aba_names, alternating_vals)
ca1_psyn_sim <- load_sim(ca1_psyn_sim_file,   aba_names, alternating_vals)

# ══════════════════════════════════════════════════════════════════════════════
# (i) CP atrophy
# ══════════════════════════════════════════════════════════════════════════════
message("\n── (i) CP atrophy ──")
res_cp_atrophy <- run_condition(
  label         = "cp_antero",
  display_title = "CP | Empirical atrophy",
  sim_df        = antero_sim,
  emp_df        = tstat_cp,
  params        = antero_params,
  snca_df       = snca_df,
  epi_index     = cp_epi_index,
  y_lab         = "Empirical atrophy (t-statistic)",
  baseline_val  = get_baseline(BASELINE_KEYS$cp_atrophy),
  is_psyn       = FALSE,
  cl_gene_df    = cl_gene_dfs$cp_atrophy,
  cl_gene_name  = "Dnajc17"
)

# ══════════════════════════════════════════════════════════════════════════════
# (ii) CP pSyn
# ══════════════════════════════════════════════════════════════════════════════
message("\n── (ii) CP pSyn ──")
res_cp_psyn <- run_condition(
  label         = "cp_antero",
  display_title = "CP | Empirical pSyn",
  sim_df        = cp_psyn_sim,
  emp_df        = psyn_cp,
  params        = antero_params,
  snca_df       = snca_df,
  epi_index     = cp_epi_index,
  y_lab         = "Empirical pSyn",
  baseline_val  = get_baseline(BASELINE_KEYS$cp_psyn),
  is_psyn       = TRUE,
  cl_gene_df    = cl_gene_dfs$cp_psyn,
  cl_gene_name  = "Zmat4"
)

# ══════════════════════════════════════════════════════════════════════════════
# (iii) DG atrophy
# ══════════════════════════════════════════════════════════════════════════════
message("\n── (iii) DG atrophy ──")
res_dg_atrophy <- run_condition(
  label         = "hipp_retro",
  display_title = "DG | Empirical atrophy",
  sim_df        = retro_sim,
  emp_df        = tstat_hipp,
  params        = retro_params,
  snca_df       = snca_df,
  epi_index     = dg_epi_index,
  y_lab         = "Empirical atrophy (t-statistic)",
  baseline_val  = get_baseline(BASELINE_KEYS$dg_atrophy),
  is_psyn       = FALSE,
  cl_gene_df    = cl_gene_dfs$dg_atrophy,
  cl_gene_name  = "Twist2"
)

# ══════════════════════════════════════════════════════════════════════════════
# (iv) CA1 pSyn
# ══════════════════════════════════════════════════════════════════════════════
message("\n── (iv) CA1 pSyn ──")
res_ca1_psyn <- run_condition(
  label         = "hipp_retro",
  display_title = "CA1 | Empirical pSyn",
  sim_df        = ca1_psyn_sim,
  emp_df        = psyn_dg,
  params        = retro_params,
  snca_df       = snca_df,
  epi_index     = ca1_epi_index,
  y_lab         = "Empirical pSyn",
  baseline_val  = get_baseline(BASELINE_KEYS$ca1_psyn),
  is_psyn       = TRUE,
  cl_gene_df    = cl_gene_dfs$ca1_psyn,
  cl_gene_name  = "Ctnnbip1"
)

# ══════════════════════════════════════════════════════════════════════════════
# Stitched 3-level barplots
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Stitched 3-level barplots ──")

bar_height <- 3.5 + 0.45 * 7   # consistent with individual saves

save_png(
  res_cp_atrophy$p_bar3 + res_dg_atrophy$p_bar3 +
    plot_annotation(
      title = "Atrophy | Spearman \u03c1 by hemisphere",
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18))
    ) +
    plot_layout(guides = "collect") & theme(legend.position = "bottom"),
  "barplot_atrophy_3level_stitched.png", width = 16, height = bar_height
)

save_png(
  res_cp_psyn$p_bar3 + res_ca1_psyn$p_bar3 +
    plot_annotation(
      title = "pSyn | Spearman \u03c1 by hemisphere",
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18))
    ) +
    plot_layout(guides = "collect") & theme(legend.position = "bottom"),
  "barplot_psyn_3level_stitched.png", width = 16, height = bar_height
)

# ══════════════════════════════════════════════════════════════════════════════
# Cross-modal
# ══════════════════════════════════════════════════════════════════════════════
message("\n── Cross-modal plots ──")
make_crossmodal_plot(psyn_cp, tstat_cp,   "cp")
make_crossmodal_plot(psyn_dg, tstat_hipp, "dg")

message("\nAll done. Figures saved to ../figures/")