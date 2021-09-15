# Yige Wu @ WashU 2020 Nov
## reference: https://www.nature.com/articles/bjc2014215

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/"
setwd(dir_base)
source("./phospho-signaling_analysis/load_pkgs.R")
source("./phospho-signaling_analysis/functions.R")
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
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------
## input the mutation data
maf_df <- fread(data.table = F, input = "./analysis_results/preprocess_files/tables/patch_pda_maf/20201125.v1/HTAN_PDAC_somatic_dnp.patched.maf")
## input phosphosite data
protein_df <- fread(data.table = F, input = "./Resources/PDA/Data/HTAN_PDA_proteome_matrix-log2_ratios-MD_norm_MAD_scaling_Formatted_NA20.txt", fill=TRUE)
## input treatment data and VAF
metadata_df <- readxl::read_excel(path = "./Resources/PDA/metadata.xlsx")

# set parameters ----------------------------------------------------------
## set plotting parameters
num_nonna <- 10
## set genes to plot
gene2pathway_plot_df <- data.frame(Gene = c("TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4",
                                            "KRAS", "RELA", "RELB", "GNAQ"),
                                   Pathway = c(rep("TGF-beta", 6),
                                               rep("RAS-RAL", 4)))

# make data matrix-------------------------------------
## get ids
colnames(protein_df)
colnames_dataid <- c("Gene")
colnames_sampleid <- colnames(protein_df)[!(colnames(protein_df) %in% colnames_dataid)]
colnames_sampleid
## get the data frame
plot_data_df <- protein_df %>%
  filter(Gene %in% gene2pathway_plot_df$Gene) %>%
  arrange(Gene)
## filter by non-na number
number_nonna <- rowSums(!is.na(plot_data_df[,colnames_sampleid]))
plot_data_df <- plot_data_df[number_nonna >= num_nonna,]
## remove duplicates
plot_data_df <- plot_data_df[!duplicated(plot_data_df$Gene),]; nrow(plot_data_df)
## transform data frame to matrix
plot_data_mat <- as.matrix(plot_data_df[, colnames_sampleid])
rownames(plot_data_mat) <- plot_data_df$Gene
## get ids
genesymbols_plot_vec <- plot_data_df$Gene

# make colors -------------------------------------------------------------
## get color corresponding to values
summary(as.vector(plot_data_mat))
colors_heatmapbody = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
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
aachange_df <- generate_variant_aachange_matrix(pair_tab = c("KRAS", "TP53", "SMAD4", "CDKN2A", "GNAQ"), maf = maf_df)
sampleids_miss <- colnames_sampleid[!(colnames_sampleid %in% colnames(aachange_df))]
sampleids_overlap <- intersect(colnames_sampleid, colnames(aachange_df))
## get KRAS mutation
kras_aachange_vec <- aachange_df[aachange_df$Hugo_Symbol == "KRAS", sampleids_overlap]
kras_aachange_vec[sampleids_miss] <- ""; kras_aachange_vec <- kras_aachange_vec[colnames_sampleid]
## get TP53 mutation
tp53_mut_vec <- as.character(aachange_df[aachange_df$Hugo_Symbol == "TP53", sampleids_overlap] != ""); names(tp53_mut_vec) <- sampleids_overlap
tp53_mut_vec[sampleids_miss] <- "FALSE"; tp53_mut_vec <- tp53_mut_vec[colnames_sampleid]
## get SMAD4 mutation
smad4_mut_vec <- as.vector(aachange_df[aachange_df$Hugo_Symbol == "SMAD4", sampleids_overlap] != ""); names(smad4_mut_vec) <- sampleids_overlap
smad4_mut_vec[sampleids_miss] <- "FALSE"; smad4_mut_vec <- smad4_mut_vec[colnames_sampleid]
## get GNAQ mutation
gnaq_mut_vec <- as.vector(aachange_df[aachange_df$Hugo_Symbol == "GNAQ", sampleids_overlap] != ""); names(gnaq_mut_vec) <- sampleids_overlap
gnaq_mut_vec[sampleids_miss] <- "FALSE"; gnaq_mut_vec <- gnaq_mut_vec[colnames_sampleid]
## get CDKN2A mutation
cdkn2a_mut_vec <- as.vector(aachange_df[aachange_df$Hugo_Symbol == "CDKN2A", sampleids_overlap] != ""); names(cdkn2a_mut_vec) <- sampleids_overlap
cdkn2a_mut_vec[sampleids_miss] <- "FALSE"; cdkn2a_mut_vec <- cdkn2a_mut_vec[colnames_sampleid]
## get treatment data
treatment_vec <- mapvalues(x = colnames_sampleid, from = metadata_df$Sample, to = as.vector(metadata_df$Treatment))
table(treatment_vec)
## get KRAS vaf data
kras_vaf_vec <- mapvalues(x = colnames_sampleid, from = metadata_df$Sample, to = as.vector(metadata_df$KRAS_VAF))
kras_vaf_vec <- as.numeric(kras_vaf_vec)
kras_vaf_vec[is.na(kras_vaf_vec)] <- 0
## get moffit subtype
moffit_subtype_vec <- mapvalues(x = colnames_sampleid, from = metadata_df$Sample, to = as.vector(metadata_df$Moffitt_Subtype))
## get KRAS vaf data
purity_vec <- mapvalues(x = colnames_sampleid, from = metadata_df$Sample, to = as.vector(metadata_df$scRNA_tumor_percentage))
purity_vec <- as.numeric(purity_vec)
## make 
col_anno_obj <- HeatmapAnnotation(Purity = anno_simple(x = purity_vec, col = colors_purity),
                                  KRAS_Mut = anno_text(x = kras_aachange_vec, gp = gpar(border = "black")),
                                  KRAS_VAF = anno_simple(x = kras_vaf_vec, col = colors_kras_vaf),
                                  TP53_Mut = anno_simple(x = tp53_mut_vec, col = colors_ifmut[tp53_mut_vec]),
                                  SMAD4_Mut = anno_simple(x = smad4_mut_vec, col = colors_ifmut[smad4_mut_vec]),
                                  GNAQ_Mut = anno_simple(x = gnaq_mut_vec, col = colors_ifmut[gnaq_mut_vec]),
                                  Treatment = anno_simple(x = treatment_vec, col = colors_treatment[treatment_vec]),
                                  Moffitt_Subtype = anno_simple(x = moffit_subtype_vec, col = colors_moffit_subtype[moffit_subtype_vec]))

# make row split ----------------------------------------------------------
row_split_vec <- mapvalues(x = genesymbols_plot_vec, from = gene2pathway_plot_df$Gene, to = as.vector(gene2pathway_plot_df$Pathway))
row_split_factor <- factor(x = row_split_vec, levels = unique(gene2pathway_plot_df$Pathway))

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
file2write <- paste0(dir_out, "TGF_RAL", ".png")
png(file2write, width = 1000, height = 1200, res = 150)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()
file2write <- paste0(dir_out, "TGF_RAL", ".pdf")
pdf(file2write, width = 10, height = 7, useDingbats = F)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()



