# ============================================================
# VIDEO: Simulated spread over time (regional brain maps)
# ============================================================
library(av)
source("plot_all_sir_input_output_vols.R")

# ── CONFIG ────────────────────────────────────────────────────────────────────
ABA_NAMES_FILE <- "../../derivatives/SIR_inputs/csvs/yohan_source_full.csv"
N_REGIONS      <- 209
N_TIMESTEPS    <- 1000L
VIDEO_FPS      <- 20L
SLICES         <- seq(-7.5, 5, 1)
BASE_FRAME_DIR <- "../figures/video_frames/"
BASE_VIDEO_DIR <- "../figures/"

EMPIRICAL_FILES <- c(
  "../../preprocessed/steph_janice_regional_jacobians/rgn_t_stats_full_hemiMsPff_hemiPBS.csv",
  "../../preprocessed/steph_janice_regional_jacobians/hipp_rgn_t_stats_hemiHuPff_hemiPBS.csv",
  "../../preprocessed/shady_ihc/CP_injection.csv",
  "../../preprocessed/shady_ihc/HIP_injection.csv"
)
EMPIRICAL_COLS <- c("Time.3", "Time.3", "Total.Pathology_24.MPI", "Total.Pathology_24.MPI")
EPIS           <- c("CP", "DG", "CP", "CA1")
PATHOLOGIES    <- c("atrophy", "atrophy", "pSyn", "pSyn")
SIM_INPUT_FILES <- c(
  "../simulations/abm_spread_v.104.87116708452713.spread_rate.0.02122582417807315.dt.0.1.seed.35.injection_amount.97.7903636896988.clearance_gene.Dnajc17.k1.0.03066890459365682.k2.0.42327816662213974..csv",
  "../simulations/retro_abm_spread_v.1.0862247470116264.spread_rate.0.027270129876801334.dt.0.1.seed.40.injection_amount.35.5593507054226.clearance_gene.Twist2.k1.0.7296383313819099.k2.0.059192558911420556..csv",
  "../simulations/abm_spread_v.2.723990344784835.spread_rate.0.05812310477637486.dt.0.1.seed.35.injection_amount.82.23688396384566.clearance_gene.Zmat4...csv",
  "../simulations/abm_spread_v.0.29188965403566514.spread_rate.0.0006938150874085919.dt.0.1.seed.24.injection_amount.54.177718195890705.clearance_gene.Ctnnbip1...csv"
)

# ── SHARED: ABA names + anatomy raster (loaded once) ─────────────────────────
load_aba_names <- function(aba_names_file, n_regions) {
  aba              <- read.csv(aba_names_file, header = TRUE)
  aba              <- rbind(aba, aba)
  aba              <- aba[, -1]
  alternating_vals <- rep(c("right", "left"), each = n_regions)
  list(aba_region_names = aba, alternating_vals = alternating_vals)
}

# ── Shared white-background brain plot theme ──────────────────────────────────
theme_brain_white <- function(base_size = 28) {
  theme_void(base_size = base_size) +
    theme(
      panel.spacing    = unit(0.01, "npc"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.title       = element_text(size = 26, face = "bold",
                                      hjust = 0.5, colour = "black"),
      plot.subtitle    = element_text(size = 22, hjust = 0.5, colour = "black"),
      strip.text       = element_text(size = 20, face = "bold", colour = "black"),
      legend.text      = element_text(colour = "black", size = 18),
      legend.title     = element_text(colour = "black", size = 24, face = "bold"),
   )
}

if (!exists("aba_shared")) {
  aba_shared <- load_aba_names(ABA_NAMES_FILE, N_REGIONS)
}
aba_names_vec    <- aba_shared$aba_region_names
alternating_vals <- aba_shared$alternating_vals

message("Pre-loading anatomy raster …")
anat_df <- prepare_masked_anatomy(allen_template_path, allen_mask_path,
                                  "y", SLICES)[[2]] %>%
  dplyr::filter(mask_value == 1)

# ── CORE FUNCTION ─────────────────────────────────────────────────────────────
render_video <- function(sim_file, empirical_file, empirical_col,
                         epi, pathology, condition_tag) {

  frame_dir <- file.path(BASE_FRAME_DIR, condition_tag)
  video_out <- file.path(BASE_VIDEO_DIR, paste0(condition_tag, "_spread.mp4"))
  dir.create(frame_dir, showWarnings = FALSE, recursive = TRUE)

  # ── Load sim (418 rows × N_TIMESTEPS cols) ─────────────────────────────────
  message(sprintf("[%s] Loading sim …", condition_tag))
  sim_raw            <- read.csv(sim_file, header = FALSE)
  sim_raw            <- cbind(aba_names_vec, sim_raw)
  sim_raw$hemisphere <- alternating_vals
  rownames(sim_raw)  <- paste(sim_raw[, ncol(sim_raw)],
                              sim_raw[, 1], sep = " ")
  sim_raw <- sim_raw[, -1]
  sim_raw <- sim_raw[, -ncol(sim_raw)]

  # ── Load empirical + compute Spearman r at every timestep ──────────────────
  message(sprintf("[%s] Computing Spearman r …", condition_tag))
  emp_raw       <- read.csv(empirical_file)
  emp_col_clean <- make.names(empirical_col)
  emp_named     <- setNames(emp_raw[[emp_col_clean]], emp_raw[, 1])
  common_rows   <- intersect(names(emp_named), rownames(sim_raw))
  emp_common    <- emp_named[common_rows]

  spearman_r <- vapply(
    seq_len(N_TIMESTEPS),
    function(t) cor(as.numeric(sim_raw[common_rows, t]),
                    as.numeric(emp_common),
                    method = "spearman"),
    numeric(1)
  )

  peak_t   <- which.max(spearman_r)
  peak_rho <- spearman_r[peak_t]
  message(sprintf("[%s] Peak r = %.3f at t = %d — rendering frames 1:%d",
                  condition_tag, peak_rho, peak_t, peak_t))

  # ── Helper: one timestep → overlay data frame ───────────────────────────────
  timestep_to_overlay_df <- function(t) {
    sim_vals        <- as.numeric(sim_raw[, t])
    names(sim_vals) <- rownames(sim_raw)

    rh_vals <- sim_vals[1:N_REGIONS]
    lh_vals <- sim_vals[(N_REGIONS + 1):(2 * N_REGIONS)]

    make_df <- function(vals) {
      sids <- region_names_to_sids(names(vals))
      data.frame(`structure-id` = as.character(sids),
                  value          = as.numeric(vals),
                  check.names    = FALSE)
    }

    vol <- make_bilateral_label_vol(
      lh_df         = make_df(lh_vals),
      rh_df         = make_df(rh_vals),
      exclude_sid   = NULL,
      lh_label_file = aba_label_file,
      rh_label_file = rh_aba_label_file
    )

    tmp_mnc <- tempfile(fileext = ".mnc")
    mincWriteVolume(vol, tmp_mnc, like = aba_label_filepath)
    ov <- prepare_masked_anatomy(tmp_mnc, allen_mask_path, "y", SLICES)[[2]] %>%
      dplyr::filter(mask_value == 1, intensity != 0)
    unlink(tmp_mnc)
    ov
  }

  # ── Colour scale parameters (pathology-specific) ───────────────────────────
  if (pathology == "atrophy") {
    scale_colours <- c("orange", "white", "cyan")
    legend_lab    <- "Atrophy (z)"
  } else {
    scale_colours <- c("#ffffcc", "#7a0177")   # yellow → purple, positive only
    legend_lab    <- "pSyn (a.u.)"
  }

  # ── Render frames 1 … peak_t only ──────────────────────────────────────────
  axis_label <- if (pathology == "atrophy") "Atrophy (z)" else "pSyn (a.u.)"
  for (t in seq_len(1)) {
  #for (t in seq_len(peak_t)) {

    if (t %% 50 == 0) message(sprintf("  [%s] Frame %d / %d", condition_tag, t, peak_t))

    overlay_df <- timestep_to_overlay_df(t)
    rho_t      <- spearman_r[t]
    is_peak    <- (t == peak_t)

    # z-score per frame; fall back to raw values if scale is degenerate
    z_vals_candidate <- tryCatch(
      scale(overlay_df$intensity)[, 1],
      warning = function(w) NULL,
      error   = function(e) NULL
    )
    degenerate <- is.null(z_vals_candidate) ||
                  all(is.na(z_vals_candidate)) ||
                  is.infinite(max(abs(z_vals_candidate), na.rm = TRUE))

    z_vals <- if (t == 1 || degenerate) overlay_df$intensity else z_vals_candidate
    overlay_df$intensity <- z_vals
    frame_lim <- max(abs(z_vals), na.rm = TRUE)
    if (frame_lim == 0 || is.na(frame_lim) || is.infinite(frame_lim)) frame_lim <- 1

    rho_label <- sprintf("%s | %s     Timestep %d / %d     Spearman \u03c1 = %.3f%s",
                         epi, pathology, t, peak_t, rho_t,
                         if (is_peak) "   \u2605 PEAK" else "")

    # Build fill scale: diverging (atrophy) or sequential (pSyn)
    if (pathology == "atrophy") {
      frame_fill_scale <- scale_fill_gradientn(
        name    = legend_lab,
        colours = scale_colours,
        values  = scales::rescale(c(-frame_lim, 0, frame_lim)),
        limits  = c(-frame_lim, frame_lim),
        oob     = scales::squish,
        guide   = guide_colourbar(
          barwidth       = 0.5, barheight = 4.5,
          title.position = "top", title.hjust = 0.5,
          title.theme    = element_text(size = 17, colour = "black",
                                        face = "bold", lineheight = 1.1),
          label.theme    = element_text(size = 16, colour = "grey30")
        )
      )
    } else {
      frame_fill_scale <- scale_fill_gradientn(
        name    = legend_lab,
        colours = scale_colours,
        limits  = c(0, frame_lim),
        oob     = scales::squish,
        guide   = guide_colourbar(
          barwidth       = 0.5, barheight = 4.5,
          title.position = "top", title.hjust = 0.5,
          title.theme    = element_text(size = 17, colour = "black",
                                        face = "bold", lineheight = 1.1),
          label.theme    = element_text(size = 16, colour = "grey30")
        )
      )
    }

    p <- ggplot(anat_df, aes(x = x, y = z)) +
      geom_raster(aes(fill = intensity), interpolate = TRUE) +
      scale_fill_gradient(low = "black", high = "white",
                          oob = scales::squish, guide = "none") +
      ggnewscale::new_scale_fill() +
      geom_raster(data = overlay_df, aes(fill = intensity), alpha = 1) +
      frame_fill_scale +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      facet_wrap(~ slice_world, ncol = 13,
                 labeller = labeller(slice_world = function(x) paste0("Slice: ", x))) +
      coord_fixed(ratio = 1) +
      theme_brain_white() +
      theme(
        plot.title = element_text(
          size   = 26,
          face   = if (is_peak) "bold" else "plain",
          hjust  = 0.5,
          colour = if (is_peak) "#B8860B" else "black"
        )
      ) +
      labs(title = rho_label)

    frame_path <- sprintf("%s/frame_%04d.png", frame_dir, t)
    ggsave(frame_path, p, width = 24, height = 2, dpi = 120)
  }

  # ── Stitch ─────────────────────────────────────────────────────────────────
  message(sprintf("[%s] Stitching %d frames → %s", condition_tag, peak_t, video_out))
  frame_paths <- sprintf("%s/frame_%04d.png", frame_dir, seq_len(peak_t))
  av::av_encode_video(frame_paths, video_out, framerate = VIDEO_FPS)
  message(sprintf("[%s] Done.", condition_tag))

  # ── Peak-timestep comparison plot (simulated vs empirical, 2 rows) ──────────
  message(sprintf("[%s] Rendering peak comparison plot …", condition_tag))

  peak_overlay_df <- timestep_to_overlay_df(peak_t)

  # Simulated: z-score (with same degenerate guard)
  z_cand <- tryCatch(scale(peak_overlay_df$intensity)[, 1],
                     warning = function(w) NULL, error = function(e) NULL)
  degen  <- is.null(z_cand) || all(is.na(z_cand)) ||
            is.infinite(max(abs(z_cand), na.rm = TRUE))
  peak_overlay_df$intensity <- if (degen) peak_overlay_df$intensity else z_cand
  sim_lim <- max(abs(peak_overlay_df$intensity), na.rm = TRUE)
  if (sim_lim == 0 || is.na(sim_lim) || is.infinite(sim_lim)) sim_lim <- 1

  # Empirical: load and map to voxels using same bilateral pipeline
  emp_vals_full <- setNames(emp_raw[[emp_col_clean]], emp_raw[, 1])

  # For atrophy use load_atrophy_named_vec ordering (LH first);
  # for pSyn rows are already in the right order
  if (pathology == "atrophy") {
    n_half   <- length(emp_vals_full) / 2
    emp_lh   <- emp_vals_full[1:n_half]       # left hemisphere is first for empirical atrophy
    emp_rh   <- emp_vals_full[(n_half + 1):length(emp_vals_full)]
  } else {
    n_half   <- length(emp_vals_full) / 2
    emp_rh   <- emp_vals_full[1:n_half]
    emp_lh   <- emp_vals_full[(n_half + 1):length(emp_vals_full)]
  }

  make_sid_df <- function(vals) {
    sids <- region_names_to_sids(names(vals))
    data.frame(`structure-id` = as.character(sids),
                value          = as.numeric(vals), check.names = FALSE)
  }
  emp_vol <- make_bilateral_label_vol(
    lh_df = make_sid_df(emp_lh), rh_df = make_sid_df(emp_rh),
    exclude_sid = NULL, lh_label_file = aba_label_file,
    rh_label_file = rh_aba_label_file
  )
  tmp_emp <- tempfile(fileext = ".mnc")
  mincWriteVolume(emp_vol, tmp_emp, like = aba_label_filepath)
  emp_overlay_df <- prepare_masked_anatomy(tmp_emp, allen_mask_path, "y", SLICES)[[2]] %>%
    dplyr::filter(mask_value == 1, intensity != 0)
  unlink(tmp_emp)
  emp_lim <- max(abs(emp_overlay_df$intensity), na.rm = TRUE)
  if (emp_lim == 0 || is.na(emp_lim) || is.infinite(emp_lim)) emp_lim <- 1

  # Tag rows for faceting into two strips
  peak_overlay_df$source <- "Simulated (peak)"
  emp_overlay_df$source  <- "Empirical"

  # Filter and scale empirical atrophy values
  emp_overlay_df$intensity <- as.numeric(emp_overlay_df$intensity)
  if (pathology == "atrophy") {
    emp_overlay_df <- emp_overlay_df %>% dplyr::filter(abs(intensity) > 0.05)
    emp_overlay_df$intensity <- scale(as.numeric(emp_overlay_df$intensity))
  }

  combined_df <- rbind(peak_overlay_df, emp_overlay_df)
  combined_df$source <- factor(combined_df$source,
                               levels = c("Simulated (peak)", "Empirical"))

  if (pathology == "atrophy") {
    comp_colours  <- c("orange", "white", "cyan")
    comp_values   <- NULL
    comp_name_sim <- "Simulated Atrophy (z)"
    comp_name_emp <- "Empirical Atrophy (z)"
  } else {
    comp_colours  <- c("#ffffcc", "#7a0177")
    comp_values   <- NULL
    comp_name_sim <- "Simulated I (a.u.)"
    comp_name_emp <- "Empirical pSyn"
  }

  # ── Split combined_df into the two sources ──────────────────────────────────
  sim_reg_df <- combined_df %>% filter(source == "Simulated (peak)")
  emp_reg_df <- combined_df %>% filter(source == "Empirical")

  sim_reg_limits <- range(sim_reg_df$intensity, na.rm = TRUE)
  emp_reg_limits <- range(emp_reg_df$intensity, na.rm = TRUE)

  make_sym_values <- function(lims) {
    scales::rescale(c(lims[1], 0, lims[2]), from = c(lims[1], lims[2]))
  }

  make_reg_row <- function(overlay_df, lims, vals = NULL, row_title, fmt_fn) {
    ggplot(anat_df, aes(x = x, y = z)) +
      geom_raster(aes(fill = intensity), interpolate = TRUE) +
      scale_fill_gradient(low = "black", high = "white",
                          oob = scales::squish, guide = "none") +
      ggnewscale::new_scale_fill() +
      geom_raster(data = overlay_df, aes(fill = intensity), alpha = 0.75) +
      scale_fill_gradientn(
        name    = row_title,
        colours = comp_colours,
        values  = vals,
        limits  = lims,
        oob     = scales::squish,
        labels  = fmt_fn,
        guide   = guide_colourbar(
          barwidth       = 0.5, barheight = 4.5,
          title.position = "top", title.hjust = 0.5,
          title.theme    = element_text(size = 17, colour = "black",
                                        face = "bold", lineheight = 1.1),
          label.theme    = element_text(size = 16, colour = "grey30")
        )
      ) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      facet_wrap(~ slice_world, nrow = 1,
                 labeller = labeller(
                   slice_world = function(x) paste0("Slice: ", x))) +
      coord_fixed(ratio = 1) +
      theme_brain_white() +
      theme(
        legend.position = "right",
        strip.text.x    = element_text(size = 14, face = "bold", colour = "grey30")
      )
  }

  p_sim_reg <- make_reg_row(
    sim_reg_df, sim_reg_limits,
    vals      = if (pathology == "atrophy") make_sym_values(sim_reg_limits) else NULL,
    row_title = comp_name_sim,
    fmt_fn    = function(x) sprintf("%.2f", x)
  )

  p_emp_reg <- make_reg_row(
    emp_reg_df, emp_reg_limits,
    vals      = if (pathology == "atrophy") make_sym_values(emp_reg_limits) else NULL,
    row_title = comp_name_emp,
    fmt_fn    = function(x) sprintf("%.2f", x)
  )

  # ── Stitch and annotate ──────────────────────────────────────────────────────
  p_comp <- (p_sim_reg / p_emp_reg) +
    plot_annotation(
      title    = sprintf("%s | %s — Simulated vs. Empirical at peak (t = %d)",
                         epi, pathology, peak_t),
      subtitle = sprintf("Spearman \u03c1 = %.3f", peak_rho),
      theme    = theme(
        plot.title      = element_text(size = 24, face = "bold",
                                       colour = "black", hjust = 0.5, vjust = 0),
        plot.subtitle   = element_text(size = 19, colour = "grey40", hjust = 0.5),
        plot.background = element_rect(fill = "white", colour = NA),
        strip.text      = element_text(size = 14, face = "bold", colour = "grey30")
      )
    )

  comp_path <- file.path(BASE_VIDEO_DIR,
                         sprintf("%s_peak_comparison.png", condition_tag))
  ggsave(comp_path, p_comp, width = 24, height = 4, dpi = 300)
  message(sprintf("[%s] Comparison plot saved: %s", condition_tag, comp_path))
}

# ── RUN ALL CONDITIONS ────────────────────────────────────────────────────────
condition_tags <- paste0(tolower(EPIS), "_", PATHOLOGIES)
# → "cp_atrophy", "dg_atrophy", "cp_psyn", "ca1_psyn"
n_conditions <- 4
for (i in c(1:n_conditions)) {
  render_video(
    sim_file       = SIM_INPUT_FILES[i],
    empirical_file = EMPIRICAL_FILES[i],
    empirical_col  = EMPIRICAL_COLS[i],
    epi            = EPIS[i],
    pathology      = PATHOLOGIES[i],
    condition_tag  = condition_tags[i]
  )
}

message("All videos complete.")