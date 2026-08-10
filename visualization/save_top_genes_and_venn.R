# Install ggVennDiagram if you haven't already
# if (!requireNamespace("ggVennDiagram", quietly = TRUE)) install.packages("ggVennDiagram")
# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("yanlinlin82/ggvenn")

library(ggplot2)
library(dplyr)
library(ggVennDiagram)
library(ggvenn)

setwd(".")

# ── Directories ───────────────────────────────────────────────────────────────
dir.create("../figures", showWarnings = FALSE)
dir.create("../tables",  showWarnings = FALSE)

baselines <- as.data.frame(read.csv("../tables/baselines.csv"))

# ── Helper: load & filter one gene list ──────────────────────────────────────
load_gene_list <- function(path, baseline_label) {
  dat <- read.csv(path, header = FALSE)[, c(1, 2)]
  colnames(dat) <- c("gene", "spearman_r")
  dat$gene      <- as.character(dat$gene)
  dat$baseline  <- baselines[baselines$pathology == baseline_label, "baseline"]
  dat %>% filter(spearman_r > baseline + 0.01)
}

# ── Helper: two-way Venn diagram ─────────────────────────────────────────────
make_venn2 <- function(list1, list2,
                       label1, label2,
                       title, out_file) {
  genes <- list(
    setNames(list(list1$gene, list2$gene), c(label1, label2))
  )[[1]]

  p <- ggVennDiagram(genes,
                     label       = "count",
                     label_alpha = 0,
                     label_size  = 13) +
    scale_fill_gradientn(colours = c("white", "black", "lightgreen")) +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 24, face = "bold",
                                  margin = margin(t = 15)),
      plot.margin  = margin(30, 30, 20, 20),
      legend.position = "none",
      plot.background = element_blank()
    ) +
    scale_x_continuous(expand = expansion(mult = 0.3))

  # Increase set-name label size (layer 3) and nudge outwards
  p$layers[[3]]$aes_params$size <- 7

  ggsave(file.path("../figures", out_file), plot = p,
         width = 10, height = 6, dpi = 300)
  message("Saved: ../figures/", out_file)
  invisible(p)
}

# ── Helper: intersection & unique calculations ────────────────────────────────
calc_intersection <- function(df1, df2) {
  shared <- intersect(df1$gene, df2$gene)
  merge(df1[df1$gene %in% shared, ],
        df2[df2$gene %in% shared, ],
        by = "gene", suffixes = c("_df1", "_df2"))
}

calc_unique <- function(df, others) {
  excl <- unlist(lapply(others, `[[`, "gene"))
  df[!(df$gene %in% excl), ]
}

##############################################################################
# ── TWO-WAY VENN DIAGRAMS ────────────────────────────────────────────────────
##############################################################################

# ── hipp_ret: DG Atrophy (redo40, 90 dpi) vs. CA1 IHC (24 mpi) ──────────────
hipp_ret_atrophy <- load_gene_list(
  "../sir_result_csvs/clearance/atrophy/Time_3hipp_rgn_t_stats_hemiHuPff_hemiPBS/retro/final_output.csv",
  "hipp_ret_redo40"
)

ca1_ret_ihc <- load_gene_list(
  "../sir_result_csvs/clearance/pSyn/Total_Pathology_24_MPIHIP_injection/retro/final_output.csv",
  "ca1_ret_ihc"
)

make_venn2(
  list1    = hipp_ret_atrophy,
  list2    = ca1_ret_ihc,
  label1   = "Atrophy",
  label2   = "pSyn",
  title    = "Hipp. Atrophy vs. pSyn",
  out_file = "hipp_ret_venn_diagram.png"
)

# ── cp_ant: CP Atrophy (redo40, 90 dpi) vs. CP IHC (24 mpi) ─────────────────
cp_ant_atrophy <- load_gene_list(
  "../sir_result_csvs/clearance/atrophy/Time_3cp_rgn_t_stats_full_hemiMsPff_hemiPBS/antero/final_output.csv",
  "cp_ant_redo40"
)

cp_ant_ihc <- load_gene_list(
  "../sir_result_csvs/clearance/pSyn/Total_Pathology_24_MPICP_injection/antero/final_output.csv",
  "cp_ant_ihc"
)

make_venn2(
  list1    = cp_ant_atrophy,
  list2    = cp_ant_ihc,
  label1   = "Atrophy",
  label2   = "pSyn",
  title    = "CP Atrophy vs. pSyn",
  out_file = "cp_ant_venn_diagram.png"
)

##############################################################################
# ── FOUR-WAY VENN DIAGRAM ────────────────────────────────────────────────────
##############################################################################

hipp_atrophy <- load_gene_list(
  "../sir_result_csvs/clearance/atrophy/Time_3hipp_rgn_t_stats_hemiHuPff_hemiPBS/retro/final_output.csv",
  "hipp_ret_redo40"
)

ca1_ihc <- load_gene_list(
  "../sir_result_csvs/clearance/pSyn/Total_Pathology_24_MPIHIP_injection/retro/final_output.csv",
  "ca1_ret_ihc"
)

cp_atrophy <- load_gene_list(
  "../sir_result_csvs/clearance/atrophy/Time_3cp_rgn_t_stats_full_hemiMsPff_hemiPBS/antero/final_output.csv",
  "cp_ant_redo40"
)

cp_ihc <- load_gene_list(
  "../sir_result_csvs/clearance/pSyn/Total_Pathology_24_MPICP_injection/antero/final_output.csv",
  "cp_ant_ihc"
)

# ── Pairwise & four-way intersections ────────────────────────────────────────
cp_atrophy_and_hipp_atrophy <- calc_intersection(cp_atrophy, hipp_atrophy)
cp_ihc_and_ca1_ihc          <- calc_intersection(cp_ihc,     ca1_ihc)
cp_atrophy_and_cp_ihc        <- calc_intersection(cp_atrophy, cp_ihc)
hipp_atrophy_and_ca1_ihc     <- calc_intersection(hipp_atrophy, ca1_ihc)

four_way_genes <- Reduce(intersect,
                         list(cp_atrophy$gene, cp_ihc$gene,
                              hipp_atrophy$gene, ca1_ihc$gene))
four_way <- cp_atrophy[cp_atrophy$gene %in% four_way_genes, ]

# ── Unique sets ───────────────────────────────────────────────────────────────
unique_cp_atrophy   <- calc_unique(cp_atrophy,   list(cp_ihc))
unique_cp_ihc       <- calc_unique(cp_ihc,       list(cp_atrophy))
unique_hipp_atrophy <- calc_unique(hipp_atrophy, list(ca1_ihc))
unique_ca1_ihc      <- calc_unique(ca1_ihc,      list(hipp_atrophy))

# ── Save all CSVs to tables/ ──────────────────────────────────────────────────
write.csv(unique_cp_atrophy,          "../tables/unique_cp_atrophy_genes_with_correlation.csv",   row.names = FALSE)
write.csv(unique_cp_ihc,              "../tables/unique_cp_ihc_genes_with_correlation.csv",        row.names = FALSE)
write.csv(unique_hipp_atrophy,        "../tables/unique_hipp_atrophy_genes_with_correlation.csv",  row.names = FALSE)
write.csv(unique_ca1_ihc,             "../tables/unique_ca1_ihc_genes_with_correlation.csv",       row.names = FALSE)
write.csv(cp_ihc_and_ca1_ihc,         "../tables/cp_ihc_and_ca1_ihc_genes_with_correlation.csv",   row.names = FALSE)
write.csv(cp_atrophy_and_hipp_atrophy,"../tables/cp_atrophy_and_hipp_atrophy_genes_with_correlation.csv", row.names = FALSE)
write.csv(cp_atrophy_and_cp_ihc,      "../tables/cp_atrophy_and_cp_ihc_genes_with_correlation.csv", row.names = FALSE)
write.csv(hipp_atrophy_and_ca1_ihc,   "../tables/hipp_atrophy_and_ca1_ihc_genes_with_correlation.csv", row.names = FALSE)
write.csv(four_way,                   "../tables/four_way_ca1_intersection_genes_with_correlation.csv", row.names = FALSE)
message("All intersection CSVs saved to ../tables/")

# ── Four-way Venn plot ────────────────────────────────────────────────────────
genes_4way <- list(
  "CP Ant. Atrophy" = cp_atrophy$gene,
  "CP Ant. IHC"     = cp_ihc$gene,
  "DG Ret. Atrophy" = hipp_atrophy$gene,
  "CA1 Ret. IHC"    = ca1_ihc$gene
)

p4 <- ggvenn(
  genes_4way,
  fill_color      = c("#0073C2FF", "green", "lightyellow", "#CD534CFF"),
  stroke_size     = 0.5,
  set_name_size   = 7,
  text_size       = 7,
  show_percentage = FALSE
) +
  labs(title = "CP Ant. & DG Ret. — Atrophy vs. IHC") +
  theme(
    plot.title  = element_text(hjust = 0.5, size = 24, face = "bold",
                               margin = margin(t = 15)),
    plot.margin = margin(30, 30, 20, 20)
  )

ggsave("../figures/four_way_venn_diagram.png", plot = p4,
       width = 12, height = 10, dpi = 300)
message("Saved: ../figures/four_way_venn_diagram.png")
