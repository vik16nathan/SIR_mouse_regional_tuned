library("pacman")
library("RMINC")
pacman::p_load(dplyr, readr, readxl, stringr, ggrepel, tidyr, ggplot2, ggnewscale, scales, robustbase)
setwd(".")
source("../ggslicer/plotting_functions/plotting_functions.R")
source("../ggslicer/plotting_functions/plotting_functions_labels.R")

# ============================================================
# PATHS
# ============================================================

allen_input_dir <- "../../preprocessed/allen_template_inputs/"
connectome_dir  <- "../../derivatives/SIR_inputs/csvs/"
vol_output_dir  <- "../../derivatives/SIR_input_vols/"
figure_dir      <- "../figures/"

aba_label_filepath    <- paste0(allen_input_dir, "AMBA_25um_resampled_50um_int.mnc")
aba_label_filepath_RH <- paste0(allen_input_dir, "AMBA_25um_resampled_50um_int_RH.mnc")
aba_region_filepath   <- paste0(allen_input_dir, "allen_ccfv3_tree_wang_2020_s2.xlsx")
allen_template_path   <- paste0(allen_input_dir, "average_template_50.mnc")
allen_mask_path       <- paste0(allen_input_dir, "mask_50um.mnc")

ABA_NAMES_FILE <- "../../derivatives/SIR_inputs/csvs/yohan_source_full.csv"
N_REGIONS      <- 209

tstat_hipp_file <- "../../preprocessed/steph_janice_regional_jacobians/hipp_rgn_t_stats_hemiHuPff_hemiPBS.csv"
tstat_cp_file   <- "../../preprocessed/steph_janice_regional_jacobians/rgn_t_stats_full_hemiMsPff_hemiPBS.csv"
TSTAT_COL       <- "Time.3"

psyn_cp_file <- "../../preprocessed/shady_ihc/CP_injection.csv"
psyn_ca1_file <- "../../preprocessed/shady_ihc/HIP_injection.csv"
PSYN_COL     <- "Total.Pathology_24.MPI"

snca_file <- "../../derivatives/yohan_ge_filt/mr10vv0.2/Snca.csv"

# ============================================================
# LOAD ATLAS INPUTS
# ============================================================

aba_label_file    <- round(mincGetVolume(aba_label_filepath))
rh_aba_label_file <- round(mincGetVolume(aba_label_filepath_RH))

aba_region_labels <- as.data.frame(read_excel(aba_region_filepath))
colnames(aba_region_labels) <- aba_region_labels[1, ]
aba_region_labels <- aba_region_labels[2:nrow(aba_region_labels), ]

mesoscale_rgn_to_label_file_rgn_dict <- as.data.frame(
  read.csv("../../preprocessed/allen_template_inputs/mesoscale_oh_rgn_to_label_file_rgn_dict.csv")
)
mesoscale_to_label_file <- setNames(
  lapply(mesoscale_rgn_to_label_file_rgn_dict[, "label.file.region"],
         function(x) as.integer(strsplit(x, ",")[[1]])),
  as.character(mesoscale_rgn_to_label_file_rgn_dict[, "mesoscale.region"])
)

label_to_idx <- split(seq_along(aba_label_file), aba_label_file)

# ============================================================
# CORE ATLAS HELPERS
# ============================================================

load_connectome <- function(filepath) {
  df <- as.data.frame(read.csv(filepath))
  rownames(df) <- df[, 1]
  df <- df[, 2:ncol(df)]
  rownames(df) <- aba_region_labels$`structure ID`[
    match(rownames(df), aba_region_labels$abbreviation)
  ]
  colnames(df) <- c(paste0(rownames(df), ".R"), paste0(rownames(df), ".L"))
  return(df)
}

abbrev_to_sid <- function(abbrev) {
  aba_region_labels$`structure ID`[match(abbrev, aba_region_labels$abbreviation)]
}

fullname_to_sid <- function(name) {
  aba_region_labels[which(aba_region_labels[["full structure name"]] == name), "structure ID"]
}

# ============================================================
# VOLUME MAPPING
# ============================================================

map_values_to_label_vol <- function(structure_IDs_with_values,
                                    aba_label_file,
                                    mesoscale_to_label_file,
                                    label_to_idx = split(seq_along(aba_label_file), aba_label_file)) {
  count_vol <- rep(0, length(aba_label_file))
  sids      <- as.character(structure_IDs_with_values$`structure-id`)
  values    <- as.numeric(structure_IDs_with_values$`value`)
  is_parent <- sids %in% names(mesoscale_to_label_file)

  direct_sids   <- sids[!is_parent]
  direct_values <- values[!is_parent]
  for (i in seq_along(direct_sids)) {
    idx <- label_to_idx[[direct_sids[i]]]
    if (!is.null(idx)) count_vol[idx] <- direct_values[i]
  }

  parent_sids   <- sids[is_parent]
  parent_values <- values[is_parent]
  for (i in seq_along(parent_sids)) {
    child_ids <- mesoscale_to_label_file[[parent_sids[i]]]
    idx <- unlist(label_to_idx[as.character(child_ids)], use.names = FALSE)
    if (length(idx)) count_vol[idx] <- parent_values[i]
  }
  return(count_vol)
}

make_bilateral_label_vol <- function(lh_df,
                                     rh_df,
                                     exclude_sid   = NULL,
                                     lh_label_file = aba_label_file,
                                     rh_label_file = rh_aba_label_file) {
  if (!is.null(exclude_sid))
    rh_df$value[rh_df$`structure-id` == exclude_sid] <- 0

  vol       <- map_values_to_label_vol(lh_df, lh_label_file, mesoscale_to_label_file, label_to_idx)
  rh_voxels <- which(rh_label_file > 0)
  rh_vol    <- map_values_to_label_vol(rh_df, lh_label_file, mesoscale_to_label_file, label_to_idx)
  vol[rh_voxels] <- rh_vol[rh_voxels]
  return(vol)
}

#' Write any 418-row ("structure-id", "value") data frame to a .mnc volume.
#'
#' This is the general entry point for all non-connectome inputs.
#' The data frame must have exactly two columns named "structure-id" and "value",
#' with one row per region (up to 418). No outlier filtering is applied here.
#'
#' @param sid_value_df  Data frame with columns "structure-id" and "value".
#' @param out_mnc_name  Filename (no path) for the output .mnc volume.
#' @return              The voxel volume vector (invisibly).
sid_value_df_to_mnc <- function(sid_value_df, out_mnc_name) {
  vol <- map_values_to_label_vol(sid_value_df, aba_label_file,
                                 mesoscale_to_label_file, label_to_idx)
  out_path <- paste0(vol_output_dir, out_mnc_name)
  mincWriteVolume(vol, out_path, like = aba_label_filepath)
  message(sprintf("Written: %s", out_path))
  invisible(vol)
}

# ============================================================
# CONNECTOME HELPERS
# ============================================================

connectome_row_to_sid_value_df <- function(connectome, source_abbrev,
                                           hemisphere = c("LH", "RH")) {
  hemisphere <- match.arg(hemisphere)
  src_sid    <- abbrev_to_sid(source_abbrev)
  src_row    <- which(rownames(connectome) == src_sid)
  n          <- nrow(connectome)
  conn_row   <- connectome[src_row, ]

  vals <- if (hemisphere == "LH") conn_row[, (n + 1):ncol(conn_row)] else conn_row[, 1:n]
  colnames(vals) <- rownames(connectome)

  data.frame(`structure-id` = colnames(vals),
              value          = as.numeric(vals),
              check.names    = FALSE)
}

apply_outlier_threshold <- function(df, n_total = 209, label = "") {
  vals    <- as.numeric(df$value)
  stats   <- robustbase::adjboxStats(vals)
  upper   <- stats$stats[5]
  n_above <- sum(vals > upper, na.rm = TRUE)
  message(sprintf("[%s] adjbox upper = %.4f | %d / %d regions clipped (%.1f%%)",
                  label, upper, n_above, n_total, 100 * n_above / n_total))
  df$value <- pmin(vals, upper)
  list(df = df, upper = upper)
}

#' Threshold both hemispheres jointly, build bilateral volume, write .mnc.
#' Returns the shared upper threshold for use in plot_sid_value_vol().
make_outlier_thresh_conn_vol <- function(conn_filename, rgn_abbr) {
  connectome <- load_connectome(paste0(connectome_dir, conn_filename))

  lh_df <- connectome_row_to_sid_value_df(connectome, rgn_abbr, "LH")
  rh_df <- connectome_row_to_sid_value_df(connectome, rgn_abbr, "RH")

  combined_vals <- c(lh_df$value, rh_df$value)
  stats   <- robustbase::adjboxStats(combined_vals)
  upper   <- stats$stats[5]
  n_above <- sum(combined_vals > upper, na.rm = TRUE)
  message(sprintf("[%s LH+RH] adjbox upper = %.4f | %d / %d regions clipped (%.1f%%)",
                  rgn_abbr, upper, n_above, length(combined_vals),
                  100 * n_above / length(combined_vals)))

  lh_df$value <- pmin(lh_df$value, upper)
  rh_df$value <- pmin(rh_df$value, upper)

  vol <- make_bilateral_label_vol(lh_df, rh_df,
                                  exclude_sid = abbrev_to_sid(rgn_abbr))
  mincWriteVolume(vol,
                  paste0(vol_output_dir, rgn_abbr, "_conn_vol_rgns.mnc"),
                  like = aba_label_filepath)
  return(upper)
}

# ============================================================
# INPUT DATA → SID/VALUE DATA FRAME HELPERS
# ============================================================

#' Convert a named vector of region names → structure IDs.
#' Region names must match "full structure name" in aba_region_labels exactly.
#' Strips leading "right " / "left " prefixes before lookup.
region_names_to_sids <- function(names) {
  clean <- sub("^(right|left) ", "", names, ignore.case = TRUE)
  sids  <- sapply(clean, fullname_to_sid)
  return(sids)
}

region_abbrev_to_sid <- function(names) {
  clean <- sub("^(right|left) ", "", names, ignore.case = TRUE)
  sids  <- sapply(clean, abbrev_to_sid)
  return(sids)
}

# ============================================================
# BILATERAL INPUT HELPERS  (atrophy, pSyn)
# ============================================================

#' Split a 418-row named vector into RH (rows 1–209) and LH (rows 210–418)
#' sid/value data frames, mirroring the connectome convention.
#'
#' @param named_values  Named numeric vector, length 418.
#'                      Names must be "right <Region>" / "left <Region>".
#' @return              List with $rh and $lh, each a ("structure-id","value") df.
n_sources <- 209
bilateral_named_vec_to_rhlh_dfs <- function(named_values) {
  rh_vals <- named_values[1:n_sources]
  lh_vals <- named_values[(n_sources+1):(2*n_sources)]

  make_df <- function(vals) {
    sids <- region_names_to_sids(names(vals))
    data.frame(`structure-id` = as.character(sids),
                value          = as.numeric(vals),
                check.names    = FALSE)
  }
  list(rh = make_df(rh_vals), lh = make_df(lh_vals))
}

#' Write a bilateral 418-row input to a .mnc volume using RH masking.
#'
#' @param named_values  Named numeric vector, length 418 ("right/left <Name>").
#' @param out_mnc_name  Filename (no path) for the output .mnc volume.
#' @return              The upper limit of the data range (for plot scaling).
bilateral_named_vec_to_mnc <- function(named_values, out_mnc_name) {
  dfs <- bilateral_named_vec_to_rhlh_dfs(named_values)

  vol <- make_bilateral_label_vol(
    lh_df         = dfs$lh,
    rh_df         = dfs$rh,
    exclude_sid   = NULL,
    lh_label_file = aba_label_file,
    rh_label_file = rh_aba_label_file
  )
  out_path <- paste0(vol_output_dir, out_mnc_name)
  mincWriteVolume(vol, out_path, like = aba_label_filepath)
  message(sprintf("Written: %s", out_path))
  invisible(vol)
}

# ── Loaders ───────────────────────────────────────────────────────────────────

load_atrophy_named_vec <- function(filepath, tstat_col) {
  df <- as.data.frame(read.csv(filepath))
  df_colnames <- colnames(df)
  ####switch right and left order so RH regions are BEFORE LH regions!!
  df <- rbind(df[((nrow(df)/2):nrow(df)),], df[(1:(nrow(df)/2)),])
  colnames(df) <- df_colnames
  setNames(df[[tstat_col]], df[, 1])   # named numeric vector, length 418
}

load_psyn_named_vec <- function(filepath, psyn_col) {
  df <- as.data.frame(read.csv(filepath))
  setNames(df[[psyn_col]], df[, 1])
}


#' Build a 418-row sid/value data frame from a unilateral input (209 regions),
#' copying the same values to both hemispheres.
#'
#' @param named_values  Named numeric vector (length 209), names are bare region
#'                      names (no "right"/"left" prefix) OR "right <Name>" format.
#' @param z_score       If TRUE, z-score the values before mapping.
#' @return              Data frame with columns "structure-id" and "value".
unilateral_named_vec_to_sid_value_df <- function(named_values, z_score = FALSE) {
  if (z_score) named_values <- scale(named_values)[, 1]

  sids <- region_names_to_sids(names(named_values))

  # Duplicate: each region appears once for RH and once for LH
  data.frame(`structure-id` = as.character(c(sids, sids)),
              value          = as.numeric(c(named_values, named_values)),
              check.names    = FALSE)
}

# ── Load and convert Snca (209 regions, z-scored, copied bilaterally) ─────────
#' Snca: 209 regions, z-scored, mirrored to both hemispheres.
load_snca_sid_value_df <- function() {
  snca_raw <- as.data.frame(read.csv(snca_file))
  rownames(snca_raw) <- snca_raw[, 1]
  vals <- scale(snca_raw$Snca)[, 1]
  names(vals) <- rownames(snca_raw)
  sids <- region_abbrev_to_sid(names(vals))
  # Mirror identically to both hemispheres — no RH masking needed
  data.frame(`structure-id` = as.character(c(sids, sids)),
              value          = as.numeric(c(vals, vals)),
              check.names    = FALSE)
}

# ── Load and convert atrophy t-stats (418 rows, bilateral, no filtering) ──────
load_atrophy_sid_value_df <- function(filepath, tstat_col) {
  df   <- as.data.frame(read.csv(filepath))
  t_values <- df[[tstat_col]]
  vals <- setNames(df[[tstat_col]], df[, 1])             # first col = region names

  bilateral_named_vec_to_sid_value_df(vals)
}

# ── Load and convert pSyn (418 rows, bilateral, no filtering) ─────────────────
load_psyn_sid_value_df <- function(filepath, psyn_col) {
  df   <- as.data.frame(read.csv(filepath))
  vals <- setNames(df[[psyn_col]], df[, 1])
  bilateral_named_vec_to_sid_value_df(vals)
}

# ============================================================
# PLOTTING
# ============================================================

plot_sid_value_vol <- function(sid_value_df,
                                out_mnc_name,
                                fig_name,
                                title        = "",
                                legend_label = "Value",
                                upper        = NULL,
                                lower        = NULL,
                                slices       = seq(-7.5, 5, 1),
                                low_col      = "blue",
                                mid_col      = NULL,
                                high_col     = "red",
                                threshold_outliers = FALSE,
                                thresh = 0) {

  # Optionally build volume from df
  if (!is.null(sid_value_df)) {
    if (threshold_outliers) {
      thresh <- apply_outlier_threshold(sid_value_df, label = title)
      upper  <- thresh$upper
      df_use <- thresh$df
    } else {
      df_use <- sid_value_df
      if (is.null(upper)) upper <- max(as.numeric(sid_value_df$value), na.rm = TRUE)
      if (is.null(lower)) lower <- min(as.numeric(sid_value_df$value), na.rm = TRUE)
    }
    sid_value_df_to_mnc(df_use, out_mnc_name)
  }

  anat_df <- prepare_masked_anatomy(allen_template_path, allen_mask_path,
                                    "y", slices)[[2]] %>% filter(mask_value == 1)

  overlay_df <- prepare_masked_anatomy(paste0(vol_output_dir, out_mnc_name),
                                       allen_mask_path, "y", slices)[[2]] %>%
    filter(mask_value == 1, abs(as.numeric(intensity)) > thresh)

  

  # Build colour scale — three-point if mid_col supplied, two-point otherwise
  if (!is.null(mid_col)) {
    fill_scale <- scale_fill_gradient2(
      name   = legend_label,
      low    = low_col,
      mid    = mid_col,
      high   = high_col,
      midpoint = 0,
      limits = c(lower, upper),
      oob    = squish,
      guide  = "colourbar"
    )
  } else {
    fill_scale <- scale_fill_gradient(
      name   = legend_label,
      low    = low_col,
      high   = high_col,
      limits = c(lower, upper),
      oob    = squish,
      guide  = "colourbar"
    )
  }

  p <- ggplot(anat_df, aes(x = x, y = z)) +
    geom_raster(aes(fill = intensity), interpolate = TRUE) +
    scale_fill_gradient(low = "black", high = "white", oob = squish, guide = "none") +
    ggnewscale::new_scale_fill() +
    geom_raster(data = overlay_df, aes(fill = intensity), alpha=0.7) +
    fill_scale +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~ slice_world, ncol = 4,
               labeller = labeller(slice_world = function(x) paste0("Slice: ", x))) +
    coord_fixed(ratio = 1) +
    theme_void(base_size = 30) +
    theme(
      panel.spacing = unit(0, "npc"),
      strip.text    = element_text(size = 16, face = "bold"),
      plot.title    = element_text(size = 24, face = "bold",   hjust = 0.5),
      plot.subtitle = element_text(size = 20, face = "italic", hjust = 0.5)
    ) +
    labs(title = title)

  out_fig <- paste0(figure_dir, fig_name)
  ggsave(out_fig, p, width = 24, height = 16, dpi = 300)
  message(sprintf("Saved: %s", out_fig))

  invisible(list(plot = p, upper = upper, lower = lower))
}

plot_connectome_heatmaps_ipsi_contra <- function(connectome,
                                                 aba_region_labels,
                                                 major_division_dict,
                                                 fig_prefix,
                                                 title = "Connectome",
                                                 log_scale = TRUE) {
  library(pheatmap)
  n <- nrow(connectome)
  conn_ipsi   <- connectome[, 1:n]
  conn_contra <- connectome[, (n + 1):(2 * n)]

  find_major_division <- function(label) {
    sample_string <- aba_region_labels[
      which(aba_region_labels[, "structure ID"] == label), "structure_id_path"
    ]
    region_hierarchy_list <- as.numeric(unlist(strsplit(sample_string, "/")))
    colnames(major_division_dict)[which(major_division_dict %in% region_hierarchy_list)]
  }

  row_ids        <- as.numeric(rownames(connectome))
  col_ids_ipsi   <- as.numeric(gsub("\\.R|\\.L", "", colnames(conn_ipsi)))
  col_ids_contra <- as.numeric(gsub("\\.R|\\.L", "", colnames(conn_contra)))

  row_md        <- sapply(row_ids, find_major_division)
  col_md_ipsi   <- sapply(col_ids_ipsi, find_major_division)
  col_md_contra <- sapply(col_ids_contra, find_major_division)

  desired_order  <- colnames(major_division_dict)
  row_idx        <- order(match(row_md, desired_order))
  col_idx_ipsi   <- order(match(col_md_ipsi, desired_order))
  col_idx_contra <- order(match(col_md_contra, desired_order))

  conn_ipsi   <- as.matrix(conn_ipsi[row_idx, col_idx_ipsi])
  conn_contra <- as.matrix(conn_contra[row_idx, col_idx_contra])
  row_md        <- row_md[row_idx]
  col_md_ipsi   <- col_md_ipsi[col_idx_ipsi]
  col_md_contra <- col_md_contra[col_idx_contra]

  fmd_row        <- factor(row_md,        levels = desired_order)
  fmd_col_ipsi   <- factor(col_md_ipsi,   levels = desired_order)
  fmd_col_contra <- factor(col_md_contra, levels = desired_order)

  annotation_row        <- data.frame(major_division = fmd_row);        rownames(annotation_row)        <- rownames(conn_ipsi)
  annotation_col_ipsi   <- data.frame(major_division = fmd_col_ipsi);   rownames(annotation_col_ipsi)   <- colnames(conn_ipsi)
  annotation_col_contra <- data.frame(major_division = fmd_col_contra); rownames(annotation_col_contra) <- colnames(conn_contra)

  iwant_hex <- c("#a83537","#4fc79c","#63348a","#73c161","#6280d6",
                 "#d3a046","#ca78cd","#4c792a","#bb467a","#aab248","#a96126","#ff846b")
  div_lvls  <- unique(c(as.character(fmd_row), as.character(fmd_col_ipsi), as.character(fmd_col_contra)))
  div_cols  <- setNames(iwant_hex[seq_along(div_lvls)], div_lvls)
  annotation_colors <- list(major_division = div_cols)

  gaps_row       <- cumsum(table(fmd_row))[-length(levels(fmd_row))]
  gaps_col_ipsi  <- cumsum(table(fmd_col_ipsi))[-length(levels(fmd_col_ipsi))]
  gaps_col_contra<- cumsum(table(fmd_col_contra))[-length(levels(fmd_col_contra))]

  transform_mat <- function(mat) if (log_scale) log10(mat + 1e-6) else mat
  cols <- colorRampPalette(c("white", "lightblue", "red"))(100)

  p_ipsi <- pheatmap(transform_mat(conn_ipsi),
                     cluster_rows = FALSE, cluster_cols = FALSE, color = cols,
                     border_color = NA, show_rownames = FALSE, show_colnames = FALSE,
                     annotation_row = annotation_row, annotation_col = annotation_col_ipsi,
                     annotation_colors = annotation_colors,
                     annotation_names_row = FALSE, annotation_names_col = FALSE,
                     gaps_row = gaps_row, gaps_col = gaps_col_ipsi,
                     main = paste0(title, " (Ipsilateral)"), fontsize = 30)
  ggsave(paste0(fig_prefix, "_ipsi.png"), p_ipsi, width = 20, height = 16, dpi = 300)

  p_contra <- pheatmap(transform_mat(conn_contra),
                       cluster_rows = FALSE, cluster_cols = FALSE, color = cols,
                       border_color = NA, show_rownames = FALSE, show_colnames = FALSE,
                       annotation_row = annotation_row, annotation_col = annotation_col_contra,
                       annotation_colors = annotation_colors,
                       annotation_names_row = FALSE, annotation_names_col = FALSE,
                       gaps_row = gaps_row, gaps_col = gaps_col_contra,
                       main = paste0(title, " (Contralateral)"), fontsize = 30)
  ggsave(paste0(fig_prefix, "_contra.png"), p_contra, width = 20, height = 16, dpi = 300)

  message("Saved ipsi + contra heatmaps")
  invisible(list(ipsi = p_ipsi, contra = p_contra))
}

# ============================================================
# CONNECTOME PLOTS
# ============================================================

major_division_dict <- data.frame(
  Isocortex=315, OLF=698, HPF=1089, CTXsp=703,
  STR=477, PAL=803, Thal=549, Hypothal=1097,
  Midbrain=313, Pons=771, Medulla=354, CB=512
)

conn_files <- list(
  rebuilt_rvm = "params_homogeneous_weights_after_qc_0413_64pct.csv",
  rebuilt_hm  = "params_regionalized_voxel_weights_after_qc_0413_64pct.csv",
  original_oh = "params_nature_yohan_ccfv3.csv"
)
conn_files_retro <- list(
  rebuilt_rvm = "params_homogeneous_weights_after_qc_0413_64pct_retro.csv",
  rebuilt_hm  = "params_regionalized_voxel_weights_after_qc_0413_64pct_retro.csv",
  original_oh = "params_nature_retro_yohan_ccfv3.csv"
)

rvm_connectome      <- load_connectome(paste0(connectome_dir, conn_files$rebuilt_rvm))
rvm_connectome_retro <- load_connectome(paste0(connectome_dir, conn_files_retro$rebuilt_rvm))

plot_connectome_heatmaps_ipsi_contra(rvm_connectome, aba_region_labels, major_division_dict,
                                     paste0(figure_dir, "rvm_connectome"),
                                     title = "Rebuilt RVM Connectome")
plot_connectome_heatmaps_ipsi_contra(rvm_connectome_retro, aba_region_labels, major_division_dict,
                                     paste0(figure_dir, "rvm_connectome_retro"),
                                     title = "Rebuilt RVM Connectome (Retro)")

# ============================================================
# CONNECTIVITY VOLUMES (outlier-thresholded)
# ============================================================

upper_cp  <- make_outlier_thresh_conn_vol(conn_files$rebuilt_rvm,       "CP")
upper_dg  <- make_outlier_thresh_conn_vol(conn_files_retro$rebuilt_rvm, "DG")
upper_ca1 <- make_outlier_thresh_conn_vol(conn_files_retro$rebuilt_rvm, "CA1")

plot_sid_value_vol(sid_value_df = NULL, out_mnc_name = "CP_conn_vol_rgns.mnc",
                   fig_name = "conn_to_cp.png",
                   title = "Rebuilt RVM Conn. Strength to CP",
                   legend_label = "Conn. strength to CP", upper = upper_cp,
                   lower = 0, low_col = "yellow", high_col = "red")

plot_sid_value_vol(sid_value_df = NULL, out_mnc_name = "DG_conn_vol_rgns.mnc",
                   fig_name = "conn_to_dg_retro.png",
                   title = "Rebuilt RVM Retro. Conn. to DG",
                   legend_label = "Conn. strength to DG", upper = upper_dg,
                   lower = 0, low_col = "yellow", high_col = "red")

plot_sid_value_vol(sid_value_df = NULL, out_mnc_name = "CA1_conn_vol_rgns.mnc",
                   fig_name = "conn_to_ca1_retro.png",
                   title = "Rebuilt RVM Retro. Conn. to CA1",
                   legend_label = "Conn. strength to CA1", upper = upper_ca1,
                   lower = 0, low_col = "yellow", high_col = "red")

# ============================================================
# SNCA / ATROPHY / PSYN VOLUMES  (no outlier filtering)
# ============================================================

# ── Snca (209 regions, z-scored, mirrored bilaterally) ────────────────────────
snca_df <- load_snca_sid_value_df()
snca_lim <- max(abs(range(as.numeric(snca_df$value), na.rm = TRUE)))

plot_sid_value_vol(
  sid_value_df = snca_df,
  out_mnc_name = "snca_zscore.mnc",
  fig_name     = "snca_zscore.png",
  title        = "Snca expression (z-score)",
  legend_label = "z-score",
  lower        = -snca_lim,
  upper        =  snca_lim,
  low_col      = "blue",
  mid_col      = "white",
  high_col     = "red"
)

# ── Atrophy — CP injection (bilateral t-stats, diverging scale) ───────────────
# ── Atrophy (bilateral, diverging, RH/LH correct) ─────────────────────────────
atrophy_cp_vals   <- load_atrophy_named_vec(tstat_cp_file,   TSTAT_COL)
atrophy_hipp_vals <- load_atrophy_named_vec(tstat_hipp_file, TSTAT_COL)

bilateral_named_vec_to_mnc(atrophy_cp_vals,   "atrophy_cp.mnc")
bilateral_named_vec_to_mnc(atrophy_hipp_vals, "atrophy_hipp.mnc")

atrophy_lim_cp   <- max(abs(range(atrophy_cp_vals,   na.rm = TRUE)))
atrophy_lim_hipp <- max(abs(range(atrophy_hipp_vals, na.rm = TRUE)))

plot_sid_value_vol(sid_value_df = NULL, out_mnc_name = "atrophy_cp.mnc",
                   fig_name = "atrophy_cp.png",
                   title = "Atrophy — CP injection (t-stat)", legend_label = "t-stat",
                   lower = -atrophy_lim_cp,   upper = atrophy_lim_cp,
                   low_col = "orange", mid_col = "white", high_col = "cyan", thresh=0.05)

plot_sid_value_vol(sid_value_df = NULL, out_mnc_name = "atrophy_hipp.mnc",
                   fig_name = "atrophy_hipp.png",
                   title = "Atrophy — Hippocampus injection (t-stat)", legend_label = "t-stat",
                   lower = -atrophy_lim_hipp, upper = atrophy_lim_hipp,
                   low_col = "orange", mid_col = "white", high_col = "cyan", thresh=0.05)

# ── pSyn (bilateral, positive only, RH/LH correct) ────────────────────────────
psyn_cp_vals  <- load_psyn_named_vec(psyn_cp_file,  PSYN_COL)
psyn_ca1_vals  <- load_psyn_named_vec(psyn_ca1_file,  PSYN_COL)

bilateral_named_vec_to_mnc(psyn_cp_vals,  "psyn_cp.mnc")
bilateral_named_vec_to_mnc(psyn_ca1_vals,  "psyn_ca1.mnc")

plot_sid_value_vol(sid_value_df = NULL, out_mnc_name = "psyn_cp.mnc",
                   fig_name = "psyn_cp.png",
                   title = "pSyn — CP injection", legend_label = "pSyn (avg. rating)",
                   lower = 0, upper = 4, low_col = "#ffffcc", high_col = "#7a0177")

plot_sid_value_vol(sid_value_df = NULL, out_mnc_name = "psyn_ca1.mnc",
                   fig_name = "psyn_ca1.png",
                   title = "pSyn — DG injection", legend_label = "pSyn (avg. rating)",
                   lower = 0, upper = 4, low_col = "#ffffcc", high_col = "#7a0177")