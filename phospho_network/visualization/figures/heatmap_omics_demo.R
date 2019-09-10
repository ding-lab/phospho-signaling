# Yige Wu @ WashU 2019 Jan
## plot a heatmap with genomics data and proteomics data and kinase-substrate score status for given pairs



# source ------------------------------------------------------------------
code_top_dir <- "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/"
path2phospho_plotting_shared <- paste0(code_top_dir, "phospho_network/phospho_network_plotting.R")
source(path2phospho_plotting_shared)

plotpheatmap <- function(mat_value, color.palette, col_anno, ann_colors, width, height, fn, if_cluster_row, if_cluster_col) {
  result <- tryCatch({
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno,
                           scale = "row",
                           na_col = "white", 
                           cluster_rows=if_cluster_row, cluster_cols=if_cluster_col, show_colnames = F, annotation_colors = ann_colors)
  }, error = function(err) {
    print(print(paste("MY_ERROR:  ",err)))
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno,
                           scale = "row",
                           na_col = "white", 
                           cluster_rows=if_cluster_row, cluster_cols=F, show_colnames = F, annotation_colors = ann_colors)
    return(my_heatmap)
  }, warning = function(war) {
    print(print(paste("MY_WARNING:  ", war)))
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno,
                           scale = "row",
                           na_col = "white", 
                           cluster_rows=if_cluster_row, cluster_cols=F, show_colnames = F, annotation_colors = ann_colors)
    
    return(my_heatmap)
  }, finally = {
    save_pheatmap_pdf(x = my_heatmap, 
                      filename = fn, 
                      width = width, height = height)
    my_heatmap <- NULL
  })
  return(result)
}

# save_pheatmap_pdf(x = my_heatmap,
#                   filename = fn,
#                   width = 20, height = 4.5)


# set variables -----------------------------------------------------------
## variables for inputting outlier score table
enzyme_type <- "kinase"
reg_nonNA <- 20
outlier_sd <- 1.5

## plotting paramters
cap <- 3
breaks = seq(-(cap),cap, by=0.2)
## add color palette
color.palette <- colorRampPalette(rev(brewer.pal(10, "RdBu")))(length(breaks))


# AKT1 signaling ----------------------------------------------------------
cancer_tmp <- "BRCA"
mut_genes <- c(SMGs[[cancer_tmp]])
cna_genes <- c("")
rna_genes <- c("")
pro_genes <- c("")
phog_genes <- c("")
pho_genes_uniq <- c( "AKT1",  "GSK3B")
pho_genes <- pho_genes_uniq
rsds <- c("S129", "S9")
row_order <- c(paste0(pho_genes, "_", rsds))
fig_width <- 20
fig_height <- 5
nonNA_cutoff <- 0
version_tmp <- 1
if_cluster_row_tmp <- F
if_cluster_col_tmp <- F

# MTOR signaling ----------------------------------------------------------
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MTOR",  "AKT1S1")
# pho_genes <- pho_genes_uniq
# rsds <- c("S2481", "S183")
# row_order <- c(paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 3
# nonNA_cutoff <- 0
# version_tmp <- 1
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T
# 
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MTOR",  "MAF1")
# pho_genes <- pho_genes_uniq
# rsds <- c("S2481", "S75")
# row_order <- c(paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 3
# nonNA_cutoff <- 0
# version_tmp <- 1
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T
# 
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MTOR",  "MAF1")
# pho_genes <- pho_genes_uniq
# rsds <- c("S2481", "T64")
# row_order <- c(paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 3
# nonNA_cutoff <- 0
# version_tmp <- 1
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T
# 
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MTOR",  "RPTOR")
# pho_genes <- pho_genes_uniq
# rsds <- c("S2481", "S863")
# row_order <- c(paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 3
# nonNA_cutoff <- 0
# version_tmp <- 1
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T
# 
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MTOR",  "UVRAG")
# pho_genes <- pho_genes_uniq
# rsds <- c("S2481", "S571")
# row_order <- c(paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 3
# nonNA_cutoff <- 0
# version_tmp <- 1
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T
# 
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MTOR",  "ZNRF2")
# pho_genes <- pho_genes_uniq
# rsds <- c("S2481", "S145")
# row_order <- c(paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 3
# nonNA_cutoff <- 0
# version_tmp <- 1
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T
# 
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MTOR",  "UVRAG",  "ZNRF2")
# pho_genes <- pho_genes_uniq
# rsds <- c("S2481", "S571", "S145")
# row_order <- c(paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 3
# nonNA_cutoff <- 0
# version_tmp <- 1
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T

# plot VHL phosphosites for CCRCC-----------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("VHL")
# rna_genes <- c("VHL")
# pro_genes <- c("VHL", "HIF1A", "EPAS1")
# phog_genes <- c("")
# pho_genes <- c(as.vector(pho_smg_tab$Gene[pho_smg_tab$Gene %in% pro_genes]))
# rsds <- c(as.vector(pho_smg_tab$Phosphosite[pho_smg_tab$Gene %in% pro_genes]))
# row_order <- c("VHL_RNA", "VHL_PRO", "HIF1A_PRO", "EPAS1_PRO",
#                paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 6
# nonNA_cutoff <- 0

# plot PBRM1 phosphosites for CCRCC-----------------------------------------------
# mut_genes <- c(SMGs[["CCRCC"]], "ATM")
# cna_genes <- c("ATM", "PBRM1")
# rna_genes <- c("ATM", "PBRM1")
# pro_genes <- c("ATM", "PBRM1")
# phog_genes <- c("")
# pho_genes_uniq <- c( "ATM", "PBRM1")
# pho_genes <- pho_genes_uniq
# rsds <- c("S1983", "S33S34S35")
# row_order <- c("ATM_RNA", "ATM_PRO", "ATM_S1983",
#                "PBRM1_RNA", "PBRM1_PRO", "PBRM1_S33S34S35")
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 20
# version_tmp <- 1
# 
# mut_genes <- c(SMGs[["CCRCC"]], "ATM")
# cna_genes <- c("ATM", "PBRM1")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("")
# pho_genes_uniq <- c( "ATM", "PBRM1")
# pho_genes <- c(as.vector(pho_tab$Gene[pho_tab$Gene %in% pho_genes_uniq]))
# rsds <- c(as.vector(pho_tab$Phosphosite[pho_tab$Gene %in% pho_genes_uniq]))
# row_order <- c(paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 6
# nonNA_cutoff <- 30
# version_tmp <- 2

# plot MET - SHC for CCRCC-----------------------------------------------
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET", "SHC1")
# rna_genes <- c("MET", "SHC1")
# pro_genes <- c("MET", "SHC1")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MET", "SHC1")
# pho_genes <- pho_genes_uniq
# rsds <- c("T977", "Y317")
# row_order <- c("MET_PRO", "MET_T977",
#               "SHC1_PRO", "SHC1_Y317")
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 20
# version_tmp <- 1
# if_cluster_row_tmp <- F
# if_cluster_col_tmp <- T

# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET", "SHC1")
# rna_genes <- c("MET", "SHC1")
# pro_genes <- c("MET", "SHC1")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MET", "SHC1")
# pho_genes <- pho_genes_uniq
# rsds <- c("T992", "Y239Y240")
# row_order <- c("MET_PRO", "MET_T992",
#               "SHC1_PRO", "SHC1_Y239Y240")
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 20
# version_tmp <- 3
# if_cluster_row_tmp <- F
# if_cluster_col_tmp <- T

# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET", "SHC1")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MET", "SHC1")
# pho_genes <- c(as.vector(pho_tab$Gene[pho_tab$Gene %in% pho_genes_uniq]))
# rsds <- c(as.vector(pho_tab$Phosphosite[pho_tab$Gene %in% pho_genes_uniq]))
# row_order <- c(paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 6
# nonNA_cutoff <- 30
# version_tmp <- 2
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- F
# 


# plot MET - PTEN - SHC ---------------------------------------------------
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET", "PTEN", "SHC1")
# rna_genes <- c("SHC1")
# pro_genes <- c("SHC1")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MET", "PTEN", "SHC1")
# pho_genes <- pho_genes_uniq
# rsds <- c("T977", "S385", "Y317")
# row_order <- c("MET_T977", "PTEN_S385", "SHC1_Y317", "SHC1_RNA", "SHC1_PRO")
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 20
# version_tmp <- 1
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T

# plot MET - CTTN  ---------------------------------------------------
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET", "PTEN", "FLT4")
# rna_genes <- c("CTTN")
# pro_genes <- c("CTTN")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MET", "CTTN")
# pho_genes <- pho_genes_uniq
# rsds <- c("T977", "T401Y421")
# row_order <- c(paste0(pho_genes, "_", rsds),
#                "CTTN_RNA", "CTTN_PRO")
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 20
# version_tmp <- 1
# if_cluster_row_tmp <- F
# if_cluster_col_tmp <- T

# plot PBRM1 - complex partners  ---------------------------------------------------
# mut_genes <- c(SMGs[["CCRCC"]])
# cna_genes <- c("")
# rna_genes <- c("PBRM1", "SMARCC1", "SMARCA2", "ARID2", "SMARCD3", "WDR77")
# pro_genes <- c("")
# phog_genes <- c("")
# pho_genes_uniq <- c( "")
# pho_genes <- pho_genes_uniq
# rsds <- c("")
# row_order <- c(paste0(rna_genes, "_", "RNA"))
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 20
# version_tmp <- 1
# if_cluster_row_tmp <- F
# if_cluster_col_tmp <- T

# plot PBRM1 - complex partners  ---------------------------------------------------
# mut_genes <- c(SMGs[["CCRCC"]])
# cna_genes <- c("PBRM1")
# rna_genes <- c("PBRM1", "SMARCC1")
# pro_genes <- rna_genes
# phog_genes <- c("")
# pho_genes_uniq <- c( "")
# pho_genes <- pho_genes_uniq
# rsds <- c("")
# row_order <- c(paste0(rna_genes, "_", "RNA"),
#                paste0(rna_genes, "_", "PRO"))
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 20
# version_tmp <- 1
# if_cluster_row_tmp <- F
# if_cluster_col_tmp <- T

# plot PBRM1 - PTM pathway  ---------------------------------------------------
# ptm_genes <- fread(input = "./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/analysis_results/phospho_network/validation/screen_PBRM1_mutation_impact_ccRCC/RNA_not_silent_PBRM1_DE_down_capped_genes/R-HSA-597592.csv", skip = 2)
# mut_genes <- c(SMGs[["CCRCC"]])
# cna_genes <- c("PBRM1")
# rna_genes <- c("PBRM1", as.vector(ptm_genes$`Gene Symbol`))
# pro_genes <- c("")
# phog_genes <- c("")
# pho_genes_uniq <- c( "")
# pho_genes <- pho_genes_uniq
# rsds <- c("")
# row_order <- c(paste0(rna_genes, "_", "RNA"))
# fig_width <- 20
# fig_height <- 8
# nonNA_cutoff <- 20
# version_tmp <- 1
# if_cluster_row_tmp <- F
# if_cluster_col_tmp <- T

# plot MET - CTNNB1  ---------------------------------------------------
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET", "FLT4")
# rna_genes <- c("CTNNB1")
# pro_genes <- c("CTNNB1")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MET",  "CTNNB1")
# pho_genes <- pho_genes_uniq
# rsds <- c("T977", "Y716")
# row_order <- c(paste0(pho_genes, "_", rsds),
#                "CTNNB1_RNA", "CTNNB1_PRO")
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 20
# version_tmp <- 1
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T
# 
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET", "FLT4")
# rna_genes <- c("CTNNB1", "MYC")
# pro_genes <- c("CTNNB1")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MET",  "CTNNB1")
# pho_genes <- pho_genes_uniq
# rsds <- c("T977", "Y716")
# row_order <- c(paste0(pho_genes, "_", rsds),
#                "CTNNB1_RNA", "CTNNB1_PRO", "MYC_RNA")
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 20
# version_tmp <- 2
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T
# 
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET", "FLT4")
# rna_genes <- c("CTNNB1")
# pro_genes <- c("CTNNB1", "MYC")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MET",  "CTNNB1")
# pho_genes <- pho_genes_uniq
# rsds <- c("T977", "Y716")
# row_order <- c(paste0(pho_genes, "_", rsds),
#                "CTNNB1_RNA", "CTNNB1_PRO", "MYC_PRO")
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 20
# version_tmp <- 3
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T

# plot MET - CTNND1  ---------------------------------------------------
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET", "FLT4")
# rna_genes <- c("CTNND1")
# pro_genes <- c("CTNND1")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MET",  "CTNND1")
# pho_genes <- pho_genes_uniq
# rsds <- c("T977", "Y221")
# row_order <- c(paste0(pho_genes, "_", rsds),
#                "CTNND1_RNA", "CTNND1_PRO")
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 20
# version_tmp <- 1
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T

# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET", "FLT4")
# rna_genes <- c("CTNND1")
# pro_genes <- c("CTNND1")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MET",  "CTNND1")
# pho_genes <- pho_genes_uniq
# rsds <- c("T977", "Y228")
# row_order <- c(paste0(pho_genes, "_", rsds),
#                "CTNND1_RNA", "CTNND1_PRO")
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 20
# version_tmp <- 2
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T

# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET", "FLT4")
# rna_genes <- c("CTNND1")
# pro_genes <- c("MET", "CTNND1", "CTNNB1")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MET",  "CTNND1", "CTNNB1")
# pho_genes <- pho_genes_uniq
# rsds <- c("T977", "Y228", "Y716")
# row_order <- c("MET_PRO", "CTNND1_PRO", "CTNNB1_PRO",
#                paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 5
# nonNA_cutoff <- 20
# version_tmp <- 3
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T

# plot MET - EGFR  ---------------------------------------------------
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET", "FLT4")
# rna_genes <- c("EGFR")
# pro_genes <- c("EGFR")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MET",  "EGFR")
# pho_genes <- pho_genes_uniq
# rsds <- c("T992", "Y1119")
# row_order <- c(paste0(pho_genes, "_", rsds),
#                "EGFR_RNA", "EGFR_PRO")
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 20
# version_tmp <- 1
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T

# plot MET - SHC1  ---------------------------------------------------
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET", "PTEN", "FLT1", "KDR", "FLT4")
# rna_genes <- c("SHC1")
# pro_genes <- c("SHC1")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MET", "FLT1", "KDR", "FLT4", "PTEN", "SHC1")
# pho_genes <- pho_genes_uniq
# rsds <- c("T977", "S1207", "S1235", "S1362", "S385", "Y317")
# row_order <- c(paste0(pho_genes, "_", rsds), 
#                "SHC1_RNA", "SHC1_PRO")
# fig_width <- 20
# fig_height <- 6
# nonNA_cutoff <- 20
# version_tmp <- 1
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T
# 
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET", "PTEN", "FLT1", "KDR", "FLT4")
# rna_genes <- c("SHC1")
# pro_genes <- c("SHC1")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MET", "FLT1", "KDR", "FLT4", "PTEN", "SHC1")
# pho_genes <- pho_genes_uniq
# rsds <- c("T977", "S1207", "S1235", "S1362", "S385", "Y317")
# row_order <- c(paste0(pho_genes, "_", rsds), 
#                "SHC1_RNA", "SHC1_PRO")
# fig_width <- 20
# fig_height <- 5
# nonNA_cutoff <- 50
# version_tmp <- 2
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T

# plot FLT4 - SHC1  ---------------------------------------------------
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET", "FLT4")
# rna_genes <- c("")
# pro_genes <- c("FLT4", "SHC1")
# phog_genes <- c("")
# pho_genes_uniq <- c("FLT4", "SHC1")
# pho_genes <- pho_genes_uniq
# rsds <- c("S1362", "Y317")
# row_order <- c("FLT4_PRO", "FLT4_S1362",
#                "SHC1_PRO", "SHC1_Y317")
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 20
# version_tmp <- 1
# if_cluster_row_tmp <- F
# if_cluster_col_tmp <- T

# plot SHC1 - RAS signaling -----------------------------------------------
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET", "FLT4")
# rna_genes <- c("SHC1")
# pro_genes <- c("SHC1", "MAP2K1")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MET", "FLT4", "SHC1", "MAP2K1", "MAPK1")
# pho_genes <- pho_genes_uniq
# rsds <- c("T977", "S1362", "Y317", "S218", "T185Y187")
# row_order <- c(paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 30
# version_tmp <- 1
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T

# plot MET - SHC1 - RAS signaling -----------------------------------------------
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET")
# rna_genes <- c("SHC1")
# pro_genes <- c("SHC1")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MET",  "SHC1", "MAP2K1", "MAPK1")
# pho_genes <- pho_genes_uniq
# rsds <- c("T977", "Y317", "S218", "T185Y187")
# row_order <- c(paste0(pho_genes, "_", rsds), 
#                "SHC1_RNA", "SHC1_PRO")
# fig_width <- 20
# fig_height <- 5
# nonNA_cutoff <- 30
# version_tmp <- 1
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T


# plot MET - GAB1 for CCRCC-----------------------------------------------
# mut_genes <- c(SMGs[["CCRCC"]], "MET")
# cna_genes <- c("MET", "GAB1")
# rna_genes <- c("MET", "GAB1")
# pro_genes <- c("MET", "GAB1")
# phog_genes <- c("")
# pho_genes_uniq <- c( "MET", "GAB1")
# pho_genes <- pho_genes_uniq
# rsds <- c("T977", "Y265")
# row_order <- c("MET_RNA", "MET_PRO", "MET_T977",
#                "GAB1_RNA", "GAB1_PRO", "GAB1_Y265")
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 20
# version_tmp <- 1
# if_cluster_row_tmp <- F
# if_cluster_col_tmp <- T

# plot TCEB1 phosphosites for CCRCC-----------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("TCEB1")
# rna_genes <- c("TCEB1")
# pro_genes <- c("TCEB1")
# phog_genes <- c("")
# pho_genes <- c(as.vector(pho_tab$Gene[pho_tab$Gene %in% pro_genes]))
# rsds <- c(as.vector(pho_tab$Phosphosite[pho_tab$Gene %in% pro_genes]))
# row_order <- c("TCEB1_RNA", "TCEB1_PRO",
#                paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 0


# plot PTEN phosphosites for CCRCC-----------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("PTEN")
# rna_genes <- c("PTEN")
# pro_genes <- c("PTEN")
# phog_genes <- c("")
# pho_genes <- c(as.vector(pho_smg_tab$Gene[pho_smg_tab$Gene == pro_genes]), "SHC1", "PTK2", "CREB1")
# rsds <- c(as.vector(pho_smg_tab$Phosphosite[pho_smg_tab$Gene == pro_genes]), "Y317", "Y576", "S271")
# row_order <- c("PTEN_RNA", "PTEN_PRO",
#                paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 5

# plot SETD2 phosphosites for CCRCC-----------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("SETD2")
# rna_genes <- c("SETD2")
# pro_genes <- c("SETD2")
# phog_genes <- c("")
# pho_genes <- c(as.vector(pho_smg_tab$Gene[pho_smg_tab$Gene == pro_genes]))
# rsds <- c(as.vector(pho_smg_tab$Phosphosite[pho_smg_tab$Gene == pro_genes]))
# row_order <- c("SETD2_RNA", "SETD2_PRO",
#                paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 5
# nonNA_cutoff <- 50

# plot PBRM1 phosphosites for CCRCC-----------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("PBRM1")
# rna_genes <- c("PBRM1")
# pro_genes <- c("PBRM1")
# phog_genes <- c("")
# pho_genes <- c(as.vector(pho_smg_tab$Gene[pho_smg_tab$Gene == pro_genes]))
# rsds <- c(as.vector(pho_smg_tab$Phosphosite[pho_smg_tab$Gene == pro_genes]))
# row_order <- c("PBRM1_RNA", "PBRM1_PRO",
#                paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 5
# nonNA_cutoff <- 50

# plot KDM5C phosphosites for CCRCC-----------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("KDM5C")
# rna_genes <- c("KDM5C")
# pro_genes <- c("KDM5C")
# phog_genes <- c("")
# pho_genes <- c(as.vector(pho_smg_tab$Gene[pho_smg_tab$Gene == pro_genes]))
# rsds <- c(as.vector(pho_smg_tab$Phosphosite[pho_smg_tab$Gene == pro_genes]))
# row_order <- c("KDM5C_RNA", "KDM5C_PRO",
#                paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 5
# nonNA_cutoff <- 0

# plot BAP1 phosphosites for CCRCC-----------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("BAP1")
# rna_genes <- c("BAP1")
# pro_genes <- c("BAP1", "BRCA1", "BARD1")
# phog_genes <- c("")
# pho_genes <- c(as.vector(pho_smg_tab$Gene[pho_smg_tab$Gene == rna_genes]))
# rsds <- c(as.vector(pho_smg_tab$Phosphosite[pho_smg_tab$Gene == rna_genes]))
# row_order <- c("BAP1_RNA", "BAP1_PRO",
#                paste0(pho_genes, "_", rsds), 
#                "BRCA1_PRO", "BARD1_PRO")
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 15

# plot TP53 phosphosites for CCRCC-----------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("TP53", "MDM4", "CDKN2A")
# rna_genes <- c("TP53")
# pro_genes <- c("TP53", "MDM4", "CDKN2A")
# phog_genes <- c("")
# pho_genes <- c(as.vector(pho_smg_tab$Gene[pho_smg_tab$Gene == rna_genes]))
# rsds <- c(as.vector(pho_smg_tab$Phosphosite[pho_smg_tab$Gene == rna_genes]))
# row_order <- c("TP53_RNA", "TP53_PRO",
#                paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 5
# nonNA_cutoff <- 0

# plot MTOR phosphosites for CCRCC-----------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("MTOR")
# rna_genes <- c("MTOR")
# pro_genes <- c("MTOR")
# phog_genes <- c("")
# pho_genes <- c(as.vector(pho_smg_tab$Gene[pho_smg_tab$Gene == pro_genes]))
# rsds <- c(as.vector(pho_smg_tab$Phosphosite[pho_smg_tab$Gene == pro_genes]))
# row_order <- c("MTOR_RNA", "MTOR_PRO",
#                paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 50

# plot PIK3CA phosphosites for CCRCC-----------------------------------------------
# mut_genes <- c(SMGs[["CCRCC"]], "PIK3CA")
# cna_genes <- c("PIK3CA")
# rna_genes <- c("PIK3CA")
# pro_genes <- c("PIK3CA")
# phog_genes <- c("")
# pho_genes <- c(as.vector(pho_smg_tab$Gene[pho_smg_tab$Gene == pro_genes]))
# rsds <- c(as.vector(pho_smg_tab$Phosphosite[pho_smg_tab$Gene == pro_genes]))
# row_order <- c("PIK3CA_RNA", "PIK3CA_PRO",
#                paste0(pho_genes, "_", rsds))
# fig_width <- 20
# fig_height <- 4
# nonNA_cutoff <- 0

# bussiness ------------------------------------------------------------------
geneA <- paste(head(unique(c(mut_genes, cna_genes)), 5), collapse = "_")
geneB <- paste(head(unique(c(rna_genes, pro_genes, pho_genes)), 5), collapse = "_")
phosphosite <- paste0(rsds, collapse = "_")

# for (cancer in c("UCEC", "BRCA", "CCRCC", "CO", "OV")) {
for (cancer in cancer_tmp) {
  subdir1 <- paste0(makeOutDir(), cancer, "/")
  dir.create(subdir1)
  
  fn <- paste0(subdir1, paste(geneA, geneB, phosphosite, sep = "_"), "_", cancer, ".pdf")
  
  if (!file.exists(fn)) {
    ann_colors <- list()
    
    # input data first because different for each cancer type --------------------------------------------------------------
    ## input mutation matrix
    maf <- loadMaf(cancer = cancer, maf_files = maf_files)
    mut_mat <- generate_somatic_mutation_matrix(pair_tab = mut_genes, maf = maf)
    # mut_mat <- fread(paste0("./cptac2p/analysis_results/phospho_network/genoalt/tables/test_mut_impact_proteome/", cancer, "_somatic_mutation.txt"), data.table = F)
    ## mutation needs to show both geneA and geneB
    mut_mat <- mut_mat[mut_mat$Hugo_Symbol %in% mut_genes,]
    
    ## input CNA matrix
    cna_tab <- loadCNAstatus(cancer = cancer)
    cna_tab <- cna_tab[cna_tab$gene %in% cna_genes, ]
    
    
    ## load RNA
    rna_tab <- loadRNA(cancer = cancer)
    rna_tab <- rna_tab[rna_tab$gene %in% rna_genes,]
    
    ## input protein data
    if (cancer %in% c("BRCA", "OV", "CO")) {
      pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
      pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
      phog_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "collapsed_PHO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
      
    } else if (cancer == "UCEC") {
      pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
      pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
      phog_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "collapsed_PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
      
    } else if (cancer == "CCRCC") {
      pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
      pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
      phog_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "collapsed_PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
    } else if (cancer == "LIHC") {
      pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
      pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
      phog_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "collapsed_PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
    }
    pro_tab <- pro_tab[pro_tab$Gene %in% pro_genes,]
    pho_tab <- pho_tab[pho_tab$Gene %in% pho_genes & pho_tab$Phosphosite %in% rsds,]
    phog_tab <- phog_tab[phog_tab$Gene %in% phog_genes,]
    
    # make the annotation columns for each sample -----------------------------
    partIDs <- colnames(pho_tab)[!(colnames(pho_tab) %in% c("Gene", "Phosphosite", "Peptide_ID"))]
    col_anno <- data.frame(partID = partIDs)
    
    if (nrow(mut_mat) > 0){
      mut_mat.m <- melt(mut_mat, id.vars = "Hugo_Symbol")
      mut_mat.m %>% head()
      colnames(mut_mat.m) <- c("Gene", "partID", "variant_class")
      
      ## distinguish by missense and truncation
      mut_mat.m$variant_class[is.na(mut_mat.m$variant_class)] <- ""
      mut_mat.m$variant_class_sim <- "other_mutation"
      mut_mat.m$variant_class_sim[mut_mat.m$variant_class == ""] <- "wild_type"
      mut_mat.m$variant_class_sim[mut_mat.m$variant_class  == "Silent"] <- "silent"
      mut_mat.m$variant_class_sim[grepl(x = mut_mat.m$variant_class, pattern = "Missense_Mutation")] <- "missense"
      mut_mat.m$variant_class_sim[grepl(x = mut_mat.m$variant_class, pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del")] <- "truncation"
      mut_mat.m$variant_class_sim[sapply(X = mut_mat.m$variant_class, FUN = function(v) (grepl(pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del", x = v) & grepl(pattern = "Missense_Mutation", x = v)))] <- "missense&truncation"
      
      for (gene in unique(mut_mat.m$Gene[mut_mat.m$variant_class_sim != "wild_type"])) {
        mut_mat2merge <- mut_mat.m[mut_mat.m$Gene == gene, c("partID", "variant_class_sim")]
        colnames(mut_mat2merge) <- c("partID", paste0("mutation.", gene))
        col_anno <- merge(col_anno, mut_mat2merge, by = c("partID"), all.x = T)
      }
    } else {
      print("no mutation!")
    }
    
    ## CNA needs to show both geneA and geneB; levels: amplification, deletion, neutral
    if (nrow(cna_tab) > 0) {
      cna_tab.m <- melt(cna_tab, id.vars = "gene")
      colnames(cna_tab.m) <- c("gene", "partID", "CNA")
      
      for (gene in intersect(cna_genes, unique(cna_tab.m$gene[cna_tab.m$CNA != "neutral"]))) {
        cna_mat2merge <- cna_tab.m[cna_tab.m$gene == gene, c("partID", "CNA")]
        colnames(cna_mat2merge) <- c("partID", paste0("CNA.", gene))
        col_anno <- merge(col_anno, cna_mat2merge, by = c("partID"), all.x = T)
      }
    } else {
      print("no CNA!")
    }
    
    ## subtype info, need to adapt when subtype info is absent
    if (cancer %in% c("BRCA", "CO", "UCEC", "CCRCC")) {
      if (cancer == "BRCA") {
        subtypes <- partID2pam50(patientID_vector = partIDs, pam50_map = loadPAM50Map())
      } else if (cancer == "CO") {
        subtypes <- partID2MSI(patientID_vector = partIDs, subtype_map = loadMSIMap())
      } else if (cancer == "UCEC") {
        subtypes <- partID2UCECsubtype(patientID_vector = partIDs)
      } else if (cancer == "CCRCC") {
        subtypes <- partID2CCRCCImmune(patientID_vector = partIDs)
      }
      subtypes2merge <- data.frame(partID = partIDs, subtype = subtypes)
      col_anno <- merge(col_anno, subtypes2merge, by = c("partID"), all.x = T)
    }
    
    ## order samples
    col_anno %>% head()
    
    for (gene in cna_genes) {
      if (paste0("CNA.", gene) %in% colnames(col_anno)) {
        col_anno <- col_anno[order(col_anno[, paste0("CNA.", gene)], decreasing = T),]
        ann_colors[[paste0("CNA.", gene)]] <-  c(amplification = "#E41A1C", deletion = "#377EB8", "neutral" = "grey")
      }
    }
    for (gene in mut_genes) {
      if (paste0("mutation.", gene) %in% colnames(col_anno)) {
        col_anno <- col_anno[order(col_anno[, paste0("mutation.", gene)], decreasing = T),]
        ann_colors[[paste0("mutation.", gene)]] <- c(missense = "#E41A1C", truncation = "#377EB8", wild_type = "white", "missense&truncation" = "#6A3D9A", other_mutation = "#FF7F00", silent = "#33A02C")
      }
    }
    if ("subtype" %in% colnames(col_anno)) {
      col_anno <- col_anno[order(col_anno$subtype, decreasing = T),]
      subtype_colors <- set2[1:length(unique(subtypes))]
      names(subtype_colors) <- unique(subtypes)
      ann_colors[["subtype"]] <- subtype_colors
    }
    if ("esscore_outlier" %in% colnames(col_anno)) {
      ann_colors[["esscore_outlier"]] <- c("TRUE" = "#E41A1C", "FALSE" = "grey")
      ann_colors[["escore_outlier"]] <- c("TRUE" = "#E41A1C", "FALSE" = "grey")
      ann_colors[["sscore_outlier"]] <- c("TRUE" = "#E41A1C", "FALSE" = "grey")
      tmp <- color.palette; names(tmp) <- seq(from = -2, to = 2, length.out = length(color.palette))
      ann_colors[["esscore"]] <- tmp
      # col_anno <- col_anno[order(col_anno$esscore, decreasing = T),]
    }
    
    
    # col_anno <- col_anno[order(col_anno$variant_class_sim),]
    col_anno %>% head()
    rownames(col_anno) <- col_anno$partID
    
    # make the matrix of values showing in heatmap ----------------------------
    sup_tab_can <- NULL
    
    if (nrow(rna_tab) > 0) {
      rna_tab.m <- melt(rna_tab, id.vars = "gene")
      colnames(rna_tab.m) <- c("Gene", "partID", "exp_value")
      rna_tab.m$Phosphosite <- "RNA"
      sup_tab_can <- rbind(sup_tab_can, rna_tab.m[,c("Gene", "Phosphosite", "partID", "exp_value")])
    }
    
    if (nrow(pro_tab) > 0) {
      pro_tab.m <- melt(pro_tab, id.vars = "Gene")
      pro_tab.m %>% head()
      colnames(pro_tab.m) <- c("Gene", "partID", "exp_value")
      pro_tab.m$Phosphosite <- "PRO"
      sup_tab_can <- rbind(sup_tab_can, pro_tab.m[,c("Gene", "Phosphosite", "partID", "exp_value")])
    }
    
    if (nrow(pho_tab) > 0) {
      pho_tab <- pho_tab[!duplicated(pho_tab[, c("Gene", "Phosphosite")]),!(colnames(pho_tab) %in% c("Peptide_ID"))]
      pho_tab.m <- melt(pho_tab, id.vars = c("Gene", "Phosphosite"))
      pho_tab.m %>% head()
      colnames(pho_tab.m) <- c("Gene", "Phosphosite", "partID", "exp_value")
      sup_tab_can <- rbind(sup_tab_can, pho_tab.m[,c("Gene", "Phosphosite", "partID", "exp_value")])
    }
    
    if (!is.null(phog_tab)) {
      if (nrow(phog_tab) > 0) {
        phog_tab.m <- melt(phog_tab, id.vars = "Gene")
        phog_tab.m %>% head()
        colnames(phog_tab.m) <- c("Gene", "partID", "exp_value")
        phog_tab.m$Phosphosite <- "collapsed_PHO"
        sup_tab_can <- rbind(sup_tab_can, phog_tab.m[,c("Gene", "Phosphosite", "partID", "exp_value")])
      }
      
    }
    
    sup_tab_can$id_row <- paste0(sup_tab_can$Gene, "_", sup_tab_can$Phosphosite)
    sup_tab_can$exp_value <- as.numeric(as.vector(sup_tab_can$exp_value))
    sup_tab_can <- unique(sup_tab_can)
    
    ## make the matrix for the heatmap body
    df_value <- dcast(data = sup_tab_can, id_row ~ partID, value.var = "exp_value")
    
    df_value %>% head()
    mat_value <- as.matrix(df_value[,-1])
    rownames(mat_value) <- df_value$id_row
    head(mat_value)
    # mat_value <- mat_value[,colSums(!is.na(mat_value)) > 0]
    mat_value %>% head()
    
    ## order the matrix column
    partIDs_ordered <- intersect(as.vector(rownames(col_anno)), colnames(mat_value))
    col_anno$partID <- NULL
    mat_value <- mat_value[, partIDs_ordered]
    
    ## order the matrix rows
    # row_order <- c(paste0(rep(unique(c(rna_genes, pro_genes, phog_genes)), 3), "_", rep(c("RNA", "PRO", "collapsed_PHO"), length(unique(c(rna_genes, pro_genes, pho_genes))) + c(0,0,0))), paste0(pho_genes, "_", rsds))
    if (length(row_order) > 1) {
      mat_value <- mat_value[intersect(row_order, rownames(mat_value)),]
      mat_value <- mat_value[rowSums(!is.na(mat_value)) >= nonNA_cutoff, ]
      mat_value <- mat_value[,colSums(!is.na(mat_value)) >= 1]
    } else {
      mat_value <- matrix(data = mat_value, nrow = 1, dimnames = list(row_order, names(mat_value)))
    }
    # fn <- paste0(subdir1, paste(unique(c(mut_genes, cna_genes, rna_genes, pro_genes, pho_genes_uniq)), collapse = "_"), "_", cancer, "_V", version_tmp, ".pdf")
    fn <- paste0(subdir1, paste(geneA, geneB, phosphosite, sep = "_"), "_", cancer,  "_V", version_tmp, ".pdf")
    
    # plotting ----------------------------------------------------------------
    
    ## reformat the outlier status column
    if ("esscore_outlier" %in% colnames(col_anno)) {
      tmp <- vector(mode = "character", length = nrow(col_anno))
      tmp[is.na(col_anno$esscore_outlier)] <- NA
      tmp[!is.na(col_anno$esscore_outlier) & col_anno$esscore_outlier == TRUE] <- "TRUE"
      tmp[!is.na(col_anno$esscore_outlier) & col_anno$esscore_outlier == FALSE] <- "FALSE"
      col_anno$esscore_outlier <- tmp
    }
    
    # plotpheatmap(mat_value, color.palette, col_anno, ann_colors, 20, 9, fn = fn)
    plotpheatmap(mat_value = mat_value, color.palette, col_anno, ann_colors, width = fig_width , height = fig_height, fn = fn, if_cluster_row = if_cluster_row_tmp, if_cluster_col = if_cluster_col_tmp)
    
  }
}


