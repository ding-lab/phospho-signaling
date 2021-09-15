# Yige Wu @ WashU 2020 Nov
## reference: https://www.nature.com/articles/bjc2014215

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "readxl"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
packages = c(
  "ggplot2",
  "ggrepel",
  "ComplexHeatmap",
  "circlize",
  "RColorBrewer"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory (optional)
source("./phospho-signaling_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------
## input the mutation data
maf_df <- fread(data.table = F, input = "./analysis_results/preprocess_files/tables/patch_pda_maf/20201125.v1/HTAN_PDAC_somatic_dnp.patched.maf")
vaf_df <- fread(data.table = F, input = "./Resources/PDA/Data/data_for_heatmap.tsv")
## input phosphosite data
phosite_df <- fread(data.table = F, input = "./Resources/PDA/Data/Formatted/HTAN_PDA_phosphosite_matrix-log2_ratios-MD_norm_MAD_scaling_Formatted_NA20.tsv")
## input treatment data and VAF
metadata_df <- readxl::read_excel(path = "./Resources/PDA/metadata.xlsx")

# set parameters ----------------------------------------------------------
## set plotting parameters
num_nonna <- 10
## set genes to plot
gene2pathway_plot_df <- data.frame(Gene = c("BRAF",
                                            "MAP2K1", "MAP2K2", 
                                            "MAPK3", "MAPK1",
                                            "PDPK1", "AKT1", "GSK3B", "AKT1S1", "RPS6KB1", "EIF4EBP1"),
                                   Pathway = c(rep("Raf/Mek/Erk", 5),
                                               rep("PI3K/Pdk1/Akt", 6)))

# make data matrix-------------------------------------
phosite_df <- phosite_df %>%
  mutate(Phosphosite.Name = paste0(Gene, "_", Phosphosite.Index.short))
## get ids
colnames(phosite_df)
colnames_dataid <- c("Gene", "Phosphosite.Index", "Flanking.Sequence.Phosphosites.6mer", "Modifications", "Phosphosite.Index.short", "Phosphosite.Name")
colnames_sampleid <- colnames(phosite_df)[!(colnames(phosite_df) %in% colnames_dataid)]
colnames_sampleid
## get the data frame
plot_data_df <- phosite_df %>%
  filter(Gene %in% c("KRAS", 
                     "BRAF",
                     "MAP2K1", "MAP2K2", 
                     "MAPK3", "MAPK1",
                     "PDPK1", "AKT1", "AKT1S1", "RPS6KB1", "EIF4EBP1")) %>%
  filter(Phosphosite.Name %in% c("PDPK1_S241", "AKT1_S473", "AKT1S1_T246", "RPS6KB1_S240S244", "EIF4EBP1_T37T46",
                                 "MAPK3_T202", "MAPK3_Y204", "MAPK3_T202Y204", "MAPK1_T185Y187", "MAPK1_T185", "MAPK1_Y187")) %>%
  filter(Phosphosite.Index.short != "") %>%
  arrange(Gene)
## filter by non-na number
number_nonna <- rowSums(!is.na(plot_data_df[,colnames_sampleid]))
plot_data_df <- plot_data_df[number_nonna >= num_nonna,]
## remove duplicates
plot_data_df <- plot_data_df[!duplicated(plot_data_df$Phosphosite.Name),]; nrow(plot_data_df)
## transform data frame to matrix
plot_data_mat <- as.matrix(plot_data_df[, colnames_sampleid])
rownames(plot_data_mat) <- plot_data_df$Phosphosite.Name
## get ids
genesymbols_plot_vec <- plot_data_df$Gene

# make colors -------------------------------------------------------------
## get color corresponding to values
# summary(unlist(plot_data_mat))
colors_heatmapbody = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))
## make colors for the mutation
colors_ifmut <- c("TRUE" = "purple", "FALSE" = "white")
## make colors for the treatment
colors_treatment <- RColorBrewer::brewer.pal(n = 4, name = "Set1")
names(colors_treatment) <- c("Chemo-RT", "FOLFIRINOX", "Naive", "Gem_Abrax")
swatch(colors_treatment)
## make colors for KRAS VAF
colors_kras_vaf <- colorRamp2(c(0, 0.3), c("white", "brown"))
## make colors for purity
colors_purity <- colorRamp2(seq(0, 70, 10), brewer.pal(n = 8, name = "YlGnBu"))
## make colors for the moffit subtype
colors_moffit_subtype <- RColorBrewer::brewer.pal(n = 6, name = "Set1")[c(5,6)]
colors_moffit_subtype <- c(colors_moffit_subtype, "grey70")
names(colors_moffit_subtype) <- c("Basal-like", "Classical", "NA")

# make column annotation --------------------------------------------------
aachange_df <- generate_variant_aachange_matrix(pair_tab = c("KRAS", "TP53", "SMAD4", "CDKN2A"), maf = maf_df)
sampleids_miss <- colnames_sampleid[!(colnames_sampleid %in% colnames(aachange_df))]
sampleids_overlap <- intersect(colnames_sampleid, colnames(aachange_df))
## get KRAS mutation
kras_aachange_vec <- aachange_df[aachange_df$Hugo_Symbol == "KRAS", sampleids_overlap]
kras_aachange_vec[sampleids_miss] <- ""; kras_aachange_vec <- kras_aachange_vec[colnames_sampleid]
## get treatment data
treatment_vec <- mapvalues(x = colnames_sampleid, from = metadata_df$Sample, to = as.vector(metadata_df$Treatment))
table(treatment_vec)
## get KRAS vaf data
kras_vaf_vec <- mapvalues(x = colnames_sampleid, from = vaf_df[vaf_df$Hugo_Symbol == "KRAS", "Sample"], to = as.vector(vaf_df[vaf_df$Hugo_Symbol == "KRAS", "VAF"]))
kras_vaf_vec[kras_vaf_vec == colnames_sampleid] <- "0"
kras_vaf_vec <- as.numeric(kras_vaf_vec)
kras_vaf_vec[is.na(kras_vaf_vec)] <- 0
## get TP53 vaf data
tp53_vaf_vec <- mapvalues(x = colnames_sampleid, from = vaf_df[vaf_df$Hugo_Symbol == "TP53", "Sample"], to = as.vector(vaf_df[vaf_df$Hugo_Symbol == "TP53", "VAF"]))
tp53_vaf_vec[tp53_vaf_vec == colnames_sampleid] <- "0"
tp53_vaf_vec <- as.numeric(tp53_vaf_vec)
tp53_vaf_vec[is.na(tp53_vaf_vec)] <- 0
## get moffit subtype
moffit_subtype_vec <- mapvalues(x = colnames_sampleid, from = metadata_df$Sample, to = as.vector(metadata_df$Moffitt_Subtype))
## get KRAS vaf data
purity_vec <- mapvalues(x = colnames_sampleid, from = metadata_df$Sample, to = as.vector(metadata_df$scRNA_tumor_percentage))
purity_vec <- as.numeric(purity_vec)
## make 
col_anno_obj <- HeatmapAnnotation(Purity = anno_simple(x = purity_vec, col = colors_purity),
                                  KRAS_Mut = anno_text(x = kras_aachange_vec, gp = gpar(border = "black")),
                                  KRAS_VAF = anno_simple(x = kras_vaf_vec, col = colors_kras_vaf),
                                  TP53_VAF = anno_simple(x = tp53_vaf_vec, col = colors_kras_vaf),
                                  Treatment = anno_simple(x = treatment_vec, col = colors_treatment[treatment_vec]),
                                  Moffitt_Subtype = anno_simple(x = moffit_subtype_vec, col = colors_moffit_subtype[moffit_subtype_vec]))

# make row split ----------------------------------------------------------
row_split_vec <- mapvalues(x = genesymbols_plot_vec, from = gene2pathway_plot_df$Gene, to = as.vector(gene2pathway_plot_df$Pathway))
row_split_factor <- factor(x = row_split_vec, levels = c("Raf/Mek/Erk", "PI3K/Pdk1/Akt"))

# make heatmap body -------------------------------------------------------
p <- Heatmap(matrix = plot_data_mat,
             col = colors_heatmapbody,
             ## column
             cluster_columns = T,
             show_column_names = T, column_names_side = "top",
             top_annotation = col_anno_obj,
             ## row
             row_split = row_split_factor,
             cluster_rows = F,
             show_heatmap_legend = F)
p

# make legend list --------------------------------------------------------
list_lgd = list(
  Legend(col_fun = colors_heatmapbody, 
         title = "log2Intensity\n(Sample/Reference)", 
         title_gp = gpar(fontsize = 10),
         labels_gp = gpar(fontsize = 10),
         legend_width = unit(4, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_purity, 
         title = "scRNA tumor percentage", 
         title_gp = gpar(fontsize = 10),
         labels_gp = gpar(fontsize = 10),
         legend_width = unit(3, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_kras_vaf, 
         title = "KRAS VAF", 
         title_gp = gpar(fontsize = 10),
         labels_gp = gpar(fontsize = 10),
         legend_width = unit(3, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"),
  Legend(labels = names(colors_ifmut),
         title = "TP53/SMAD4\nMutation Status", 
         legend_gp = gpar(fill = colors_ifmut), 
         labels_gp = gpar(fontsize=10), title_gp = gpar(12),
         grid_height = unit(0.5, "cm"),
         grid_width = unit(0.5, "cm")),
  Legend(labels = names(colors_treatment),
         title = "Treatment Status", 
         legend_gp = gpar(fill = colors_treatment), 
         labels_gp = gpar(fontsize=10), title_gp = gpar(fontsize = 12),
         grid_height = unit(0.5, "cm"),
         grid_width = unit(0.5, "cm")),
  Legend(labels = names(colors_moffit_subtype),
         title = "Moffitt Subtype", 
         legend_gp = gpar(fill = colors_moffit_subtype), 
         labels_gp = gpar(fontsize=10), title_gp = gpar(fontsize = 12),
         grid_height = unit(0.5, "cm"),
         grid_width = unit(0.5, "cm")))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "KRAS_signaling", ".png")
png(file2write, width = 1000, height = 1200, res = 150)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()
file2write <- paste0(dir_out, "KRAS_signaling", ".pdf")
pdf(file2write, width = 10, height = 7, useDingbats = F)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()



