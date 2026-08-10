source("/data/chamal/projects/natvik/sir_voxel/analysis/ggslicer/plotting_functions/plotting_functions.R")
source("/data/chamal/projects/natvik/sir_voxel/analysis/ggslicer/plotting_functions/plotting_functions_labels.R")
library("tidyverse")
library("RMINC")
library("readxl")
library("dplyr")
library("readr")
library("scales")
library("ggnewscale")
library("patchwork")
aba_label_file_dir <- "/data/chamal/projects/yohan/common/allenbrain/data/mouse_atlas/ccfv3/"
aba_label_path_original <- "/data/chamal/projects/yohan/common/allenbrain/data/mouse_atlas/ccfv3/AMBA_relabeled_25um.mnc"
aba_label_file_original <- mincGetVolume(aba_label_path_original)
aba_regions <- unique(round(aba_label_file_original))
aba_regions <- aba_regions[which(aba_regions != 0)]


###load 200um ABA label file###

aba_label_path <- paste0(aba_label_file_dir, "AMBA_relabeled_25um_resampled_200um.mnc")
aba_label_file <- mincGetVolume(aba_label_path)
aba_label_file_rounded <- round(aba_label_file)
aba_label_file_integer_voxels <- which(aba_label_file_rounded %in% aba_regions)

###load CP and DG empirical atrophy data########
results_dir <- "/data/chamal/projects/natvik/sir_extended/analysis/"
cp_rgn_atrophy <- as.data.frame(read.csv(paste0(results_dir, "rgn_t_stats_full_hemiMsPff_hemiPBS.csv")))
hipp_rgn_atrophy <- as.data.frame(read.csv(paste0(results_dir, "hipp_rgn_t_stats_hemiMsPff_hemiPBS.csv")))

cp_rgn_atrophy_90dpi <- as.data.frame(cp_rgn_atrophy[, "Time.3"])
rownames(cp_rgn_atrophy_90dpi) <- cp_rgn_atrophy[,"X"]

dg_rgn_atrophy_90dpi <- as.data.frame(hipp_rgn_atrophy[, "Time.3"])
rownames(dg_rgn_atrophy_90dpi) <- hipp_rgn_atrophy[,"X"]

###reuse plot_shady_ihc.R####

###load from previous file
output_dir <- "/data/chamal/projects/natvik/sir_voxel/analysis/"

aba_world_coord_df <- readRDS(paste0(output_dir,"allen_voxel_index_to_world_coords_200um.Rds"))

allen_label_names <- as.data.frame(read_csv("/data/chamal/projects/yohan/common/allenbrain/data/mouse_atlas/ccfv3/AMBA_relabeled_defs_byname.csv"))

#load region hierarchy (filter for subregions too)
aba_region_filepath <- paste0(output_dir, "allen_ccfv3_tree_wang_2020_s2.xlsx")
aba_region_df <- as.data.frame(read_excel(aba_region_filepath))
colnames(aba_region_df) <- aba_region_df[1,]
aba_region_df <- aba_region_df[2:nrow(aba_region_df),]

##get rid of commas for easier merging by names
aba_region_df[,"full structure name"] <- gsub(",","",aba_region_df[,"full structure name"])
aba_region_df[which(aba_region_df[,"full structure name"] == "Nucleus of reuniens"), "full structure name"] <- "Nucleus of reunions"

###note: can't resolve the two parts of the subiculum using this 
#extract all the label names from the label file template at 25 um
##look at all 672 structures, iterating to find parents if needed
nreg <- 211 ###NOT 213 mesoscale regions bc subiculum 
mesoscale_region_names <- substring(rownames(cp_rgn_atrophy_90dpi)[1:nreg],6,)
allen_labels_total <- aba_region_df %>% filter(aba_region_df[,"structure ID"] %in% aba_regions) 
mesoscale_region_numbers <- aba_region_df[which(aba_region_df[,"full structure name"] %in% mesoscale_region_names),"structure ID"]
allen_label_names_total <- allen_labels_total[,"full structure name"]
aba_parents_total <- allen_labels_total[,"parent_id"]
aba_parent_names_total <- aba_region_df %>% filter(aba_region_df[,"structure ID"] %in% aba_parents_total)
aba_parent_names_total <- aba_parent_names_total[, "full structure name"]
#get rid of L/R hemisphere labels, and then convert to region numbers
length(intersect(mesoscale_region_names, allen_label_names_total)) ###153
length(intersect(mesoscale_region_names, aba_parent_names_total)) ###58

##########TWO MISSING REGIONS OUT OF 213: Subiculum dorsal part and Subiculum ventral part#################
region_cp_label_output <- rep(0, length(aba_label_file))
region_dg_label_output <- rep(0, length(aba_label_file))

for(structure in aba_regions) {
  region_name <- aba_region_df[which(aba_region_df[,"structure ID"] == structure),"full structure name"]
  parent_id <-aba_region_df[which(aba_region_df[,"structure ID"] == structure),"parent_id"]
  parent_name <- aba_region_df[which(aba_region_df[,"structure ID"] == parent_id), "full structure name"]
  
  ###option 1: structure is in mesoscale regions
  ###option 2: structure's parent in mesoscale regions
  if(structure %in% mesoscale_region_numbers) {
    region_voxels <- which(aba_label_file_rounded == structure)
    region_world_coords <- aba_world_coord_df[which(aba_world_coord_df[,"voxel_index"] %in% region_voxels),] 
    lh_voxels <- region_world_coords[which(region_world_coords[,"x"] < 0), "voxel_index"]
    rh_voxels <- region_world_coords[which(region_world_coords[,"x"] > 0), "voxel_index"]
    
    ##write outputs
    region_cp_label_output[lh_voxels] <- cp_rgn_atrophy_90dpi[paste0("left ",region_name),]
    region_cp_label_output[rh_voxels] <- cp_rgn_atrophy_90dpi[paste0("right ",region_name),]
    
    region_dg_label_output[lh_voxels] <- dg_rgn_atrophy_90dpi[paste0("left ",region_name),]
    region_dg_label_output[rh_voxels] <- dg_rgn_atrophy_90dpi[paste0("right ",region_name),]
    
  }
  if(length(parent_id) == 0) {
    next
  }
  if(parent_id %in% mesoscale_region_numbers) {
    ###repeat code from above, except use the parent's value 
    region_voxels <- which(aba_label_file_rounded == structure)
    region_world_coords <- aba_world_coord_df[which(aba_world_coord_df[,"voxel_index"] %in% region_voxels),] 
    lh_voxels <- region_world_coords[which(region_world_coords[,"x"] < 0), "voxel_index"]
    rh_voxels <- region_world_coords[which(region_world_coords[,"x"] > 0), "voxel_index"]
    
    ##write outputs
    region_cp_label_output[lh_voxels] <- cp_rgn_atrophy_90dpi[paste0("left ",parent_name),]
    region_cp_label_output[rh_voxels] <- cp_rgn_atrophy_90dpi[paste0("right ",parent_name),]
    
    region_dg_label_output[lh_voxels] <- dg_rgn_atrophy_90dpi[paste0("left ",parent_name),]
    region_dg_label_output[rh_voxels] <- dg_rgn_atrophy_90dpi[paste0("right ",parent_name),]
    
  }
}

##write output volumes
output_dir <- "/data/chamal/projects/natvik/sir_extended/analysis/"
setwd(output_dir)
mincWriteVolume(region_cp_label_output, "cp_atrophy_90dpi_parcellated.mnc", like=aba_label_path)
mincWriteVolume(region_dg_label_output, "dg_atrophy_90dpi_parcellated.mnc", like=aba_label_path)
####################################################################################################
###read and visualize; start with CP
####note: these maps are PBS - PFF (positive stat = volume loss)
allen_template_path_200um <- paste0(aba_label_file_dir, "average_template_200.mnc")
allen_mask_path_200um <- paste0(aba_label_file_dir, "mask_200um.mnc")
allen_200um_template_df <- prepare_masked_anatomy(allen_template_path_200um, allen_mask_path_200um, "y", seq(-7.5, 5.5, 1.5))[[2]]
setwd(output_dir)
region_atrophy_df <- prepare_masked_anatomy("cp_atrophy_90dpi_parcellated.mnc", allen_mask_path_200um, "y", seq(-7.5, 5.5, 1.5))[[2]]

allen_200um_template_df <- allen_200um_template_df %>% filter(mask_value == 1)
region_atrophy_df <- region_atrophy_df %>% filter(mask_value == 1) %>% filter(intensity > 0)

# Combine "slice_world" with "y" for unique facet labels
region_atrophy_df <- region_atrophy_df %>% 
  mutate(facet_label = paste0("y = ", y))

z_min <- min(region_atrophy_df$intensity) ###this is actually a t-stat but this code was snipped from a 
##similar analysis with z-scores
z_max <- max(region_atrophy_df$intensity)
z_mean <- mean(region_atrophy_df$intensity)
scaled_min <- 0
scaled_mid <- z_mean / (z_max - z_min)
scaled_max <- 1



ggplot(data = allen_200um_template_df, mapping = aes(x = x, y = z)) +
  geom_raster(mapping = aes(fill = intensity), interpolate = TRUE) +
  scale_fill_gradient(
    low = 'black',
    high = 'white',
    oob = scales::squish,
    guide = 'none'
  ) +
  ggnewscale::new_scale_fill() + 
  geom_raster(data = region_atrophy_df, mapping = aes(fill = intensity)) +
  scale_fill_gradientn(
    name = "Atrophy (t)",
    colours = c("darkorange", "white", "blue"),
    values = c(scaled_min, scaled_mid, scaled_max),
    limits = c(z_min, z_max),
    breaks = c(z_min, z_mean, z_max),
    labels = round(c(z_min, z_mean, z_max), 2),
    oob = scales::squish
  )+ 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~ slice_world, ncol = 9, labeller = labeller(slice_world = function(x) paste0("Slice: ", x))) +  
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(
    legend.title=element_text(size=16),
    legend.text=element_text(size=14),
    panel.spacing = unit(0, "npc"),                     # Keep spacing minimal
    strip.text = element_text(size = 16, face = "bold"), # Increase facet label text size
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5), # Make title bigger and center it
    plot.subtitle = element_text(size = 20, hjust = 0.5, face = "italic") # Make subtitle bigger and center it
  ) +
  labs(
    title = "Hemi-PBS vs. Hemi-PFF, CP injection", # Add title
  ) 
