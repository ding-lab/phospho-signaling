# Yige Wu @ WashU 2019 Jan
## plot a heatmap with genomics data and proteomics data and kinase-substrate score status for given pairs



# source ------------------------------------------------------------------
setwd(dir = "~/Box Sync/")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
library(dplyr)

plotpheatmap <- function(mat_value, color.palette, col_anno, ann_colors, width, height, fn) {
  result <- tryCatch({
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno,
                           scale = "row",
                           na_col = "white", 
                           cluster_rows=F, cluster_cols=T, show_colnames = F, annotation_colors = ann_colors)
  }, error = function(err) {
    print(print(paste("MY_ERROR:  ",err)))
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno,
                           scale = "row",
                           na_col = "white", 
                           cluster_rows=F, cluster_cols=F, show_colnames = F, annotation_colors = ann_colors)
    return(my_heatmap)
  }, warning = function(war) {
    print(print(paste("MY_WARNING:  ", war)))
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno,
                           scale = "row",
                           na_col = "white", 
                           cluster_rows=F, cluster_cols=F, show_colnames = F, annotation_colors = ann_colors)
    
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

# plot VEGF pathway for CCRCC-----------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("HIF1A", "CDKN2A", "PTEN", "NEGR1", "QKI", "CADM2", "PTPRD", "NRXN3",
#                "PRKCI", "MECOM", "MDM4", "MYC", "JAK2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("KDR", "FLT1", "FLT4")
# pho_genes <- c("SHC1")
# rsds <- c("Y317")
# row_order <- c("KDR_collapsed_PHO", "FLT1_collapsed_PHO", "FLT4_collapsed_PHO", "SHC1_Y317")

# plot VEGF pathway for CCRCC-----------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("HIF1A", "CDKN2A", "PTEN", "NEGR1", "QKI", "CADM2", "PTPRD", "NRXN3",
#                "PRKCI", "MECOM", "MDM4", "MYC", "JAK2")
# rna_genes <- c("")
# pro_genes <- c("KDR", "FLT1", "FLT4")
# phog_genes <- c("KDR", "FLT1", "FLT4")
# pho_genes <- c("SHC1")
# rsds <- c("Y317")
# row_order <- c("KDR_PRO", "KDR_collapsed_PHO", 
#                "FLT1_PRO", "FLT1_collapsed_PHO", 
#                "FLT4_PRO", "FLT4_collapsed_PHO", "SHC1_Y317")

# plot VEGF pathway for CCRCC-----------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("HIF1A", "CDKN2A", "PTEN", "NEGR1", "QKI", "CADM2", "PTPRD", "NRXN3",
#                "PRKCI", "MECOM", "MDM4", "MYC", "JAK2")
# rna_genes <- c("")
# pro_genes <- c("KDR", "FLT1", "FLT4")
# phog_genes <- c("KDR", "FLT1", "FLT4")
# pho_genes <- c("SHC1")
# rsds <- c("Y317")
# row_order <- c("KDR_PRO", "KDR_collapsed_PHO", 
#                "FLT1_PRO", "FLT1_collapsed_PHO", 
#                "FLT4_PRO", "FLT4_collapsed_PHO", "SHC1_Y317")

# plot hypoxia pathway for CCRCC-----------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("HIF1A", "CDKN2A", "PTEN", "NEGR1", "QKI", "CADM2", "PTPRD", "NRXN3",
#                "PRKCI", "MECOM", "MDM4", "MYC", "JAK2")
# rna_genes <- c("")
# pro_genes <- c("HIF2A")
# phog_genes <- c("MAPK1")
# pho_genes <- c("EP300")
# rsds <- c("S12")
# row_order <- c("HIF2A_PRO", "MAPK1_collapsed_PHO", "EP300_S12")

# plot MTOR pathway for CCRCC-----------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("HIF1A", "CDKN2A", "PTEN", "NEGR1", "QKI", "CADM2", "PTPRD", "NRXN3",
#                "PRKCI", "MECOM", "MDM4", "MYC", "JAK2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("MTOR")
# pho_genes <- c("EIF4EBP1")
# rsds <- c("S65")
# row_order <- c("MTOR_collapsed_PHO", "EIF4EBP1_S65")

# plot 71406        Pyruvate metabolism and Citric Acid (TCA) cycle pathway for CCRCC-----------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("HIF1A", "CDKN2A", "PTEN", "NEGR1", "QKI", "CADM2", "PTPRD", "NRXN3",
#                "PRKCI", "MECOM", "MDM4", "MYC", "JAK2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("PDK3")
# pho_genes <- c("PDHA1")
# rsds <- c("T269S270")
# row_order <- c("PDK3_collapsed_PHO", "PDHA1_T269S270")

# plot 71406        Pyruvate metabolism and Citric Acid (TCA) cycle pathway for CCRCC-----------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("HIF1A", "CDKN2A", "PTEN", "NEGR1", "QKI", "CADM2", "PTPRD", "NRXN3",
#                "PRKCI", "MECOM", "MDM4", "MYC", "JAK2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("MAPK1", "PKN1")
# pho_genes <- c("PKM", "PGAM1")
# rsds <- c("S151", "S14")
# row_order <- c("MAPK1_collapsed_PHO", "PKM_S151",
#                "PKN1_collapsed_PHO", "PGAM1_S14")

# # Plot cell cycle for CCRCC---------------------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("HIF1A", "CDKN2A", "PTEN", "NEGR1", "QKI", "CADM2", "PTPRD", "NRXN3",
#                "PRKCI", "MECOM", "MDM4", "MYC", "JAK2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("CDK1")
# pho_genes <- c("MCM4", "RB1", "RRM2", "TOP2A")
# rsds <- c("S120", "S249", "S20", "S1106")
# row_order <- c("CDK1_collapsed_PHO", "MCM4_S120", "RB1_S249", "RRM2_S20", "TOP2A_S1106")

# # Plot NF-kb for ccRCC ---------------------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("HIF1A", "CDKN2A", "PTEN", "NEGR1", "QKI", "CADM2", "PTPRD", "NRXN3",
#                "PRKCI", "MECOM", "MDM4", "MYC", "JAK2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("IKBKG", "MAP2K4", "MAPK14", "MAPKAPK2")
# pho_genes <- c("RELA", "MAPK14", "ATF2", "MAPKAPK2", "TAB3")
# rsds <- c("S202", "T241", "S90", "T222", "S101")
# row_order <- c("IKBKG_collapsed_PHO", "RELA_S202",
#                "MAP2K4_collapsed_PHO", "MAPK14_T241",
#                "MAPK14_collapsed_PHO", "ATF2_S90", "MAPKAPK2_T222",
#                "MAPKAPK2_collapsed_PHO", "TAB3_S101")

# Plot Cytokine Signaling in Immune system ---------------------------------------------------------
# mut_genes <- SMGs[["CCRCC"]]
# cna_genes <- c("HIF1A", "CDKN2A", "PTEN", "NEGR1", "QKI", "CADM2", "PTPRD", "NRXN3",
#                "PRKCI", "MECOM", "MDM4", "MYC", "JAK2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("TYK2", "LYN", "PTPN6", "JAK2")
# pho_genes <- c("STAT1", "STAT5A", "JAK2", "PTPN11", "GAB2", "STAT3")
# rsds <- c("S727", "S128", "T522", "T222", "Y546", "Y705")
# row_order <- c("TYK2_collapsed_PHO", "STAT1_S727",
#                "LYN_collapsed_PHO", "STAT5A_S128",
#                "PTPN6_collapsed_PHO", "JAK2_T522", "STAT3_Y705",
#                "JAK2_collapsed_PHO", "PTPN11_Y546", "GAB2_S330")

# Plot Cytokine Signaling in Immune system ---------------------------------------------------------
mut_genes <- SMGs[["CCRCC"]]
cna_genes <- c("HIF1A", "CDKN2A", "PTEN", "NEGR1", "QKI", "CADM2", "PTPRD", "NRXN3",
               "PRKCI", "MECOM", "MDM4", "MYC", "JAK2")
rna_genes <- c("")
pro_genes <- c("")
phog_genes <- c("TYK2")
pho_genes <- c("STAT1")
rsds <- c("S727")
row_order <- c("TYK2_collapsed_PHO", "STAT1_S727")

# Plot 68882        Mitotic Anaphase for UCEC ---------------------------------------------------------
# mut_genes <- SMGs[["UCEC"]]
# cna_genes <- c("LRP1B",
#                "FGFR3", "FGFR1", "ERBB2", "IGF1R", "SOX17", "FGFR3", "CCNE1", "MYC")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("CDK1")
# pho_genes <- c("TOP2B","TP53BP1")
# rsds <- c("T1292", "S1431", "S41")
# row_order <- c("CDK1_collapsed_PHO", "TOP2B_T1292", "TP53BP1_S1431")

# # Plot 68882        Mitotic Anaphase for UCEC ---------------------------------------------------------
# mut_genes <- SMGs[["UCEC"]]
# cna_genes <- c("LRP1B",
#                "FGFR3", "FGFR1", "ERBB2", "IGF1R", "SOX17", "FGFR3", "CCNE1", "MYC")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("PLK1", "BUB1")
# pho_genes <- c("CLIP1","BUB1B", "CDC20")
# rsds <- c("T287", "S670", "S41")
# row_order <- c("PLK1_collapsed_PHO", "CLIP1_T287", "BUB1B_S670",
#                "BUB1_collapsed_PHO", "CDC20_S41")

# # Plot 187687        Signalling to ERKs for UCEC ---------------------------------------------------------
# mut_genes <- SMGs[["UCEC"]]
# cna_genes <- c("LRP1B",
#                "FGFR3", "FGFR1", "ERBB2", "IGF1R", "SOX17", "FGFR3", "CCNE1", "MYC")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("SRC", "MAPK1", "MAPK3")
# pho_genes <- c("RAF1", "BRAF", "MAPKAPK2", "SOS1")
# rsds <- c("S12", "S151", "Y225", "S1137")
# row_order <- c("MAPK1_collapsed_PHO", "BRAF_S151",
#                "MAPK3_collapsed_PHO", "SOS1", "S1137")

# # Plot 187687        Signalling to ERKs for UCEC ---------------------------------------------------------
# mut_genes <- SMGs[["UCEC"]]
# cna_genes <- c("LRP1B",
#                "FGFR3", "FGFR1", "ERBB2", "IGF1R", "SOX17", "FGFR3", "CCNE1", "MYC")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("EGFR", "PTEN", "MAPK1")
# pho_genes <- c("SHC1", "SHC1", "EGFR")
# rsds <- c("Y428", "S139", "T693")
# row_order <- c("MAPK1_collapsed_PHO", "EGFR_T693",
#                "EGFR_collapsed_PHO", "SHC1_Y428", 
#                "PTEN_collapsed_PHO", "SHC1_S139")


# plot MTOR pathway for UCEC -----------------------------------------------
# mut_genes <- SMGs[["UCEC"]]
# cna_genes <- c("LRP1B",
#                "FGFR3", "FGFR1", "ERBB2", "IGF1R", "SOX17", "CCNE1", "MYC")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("MTOR", "AKT3","PRKAA2")
# pho_genes <- c("AKT1", "RPTOR", "TSC2")
# rsds <- c("S126", "T706", "S1452")
# row_order <- c("MTOR_collapsed_PHO", "AKT1_S126",
#                "AKT3_collapsed_PHO", "RPTOR_T706",
#                "PRKAA2_collapsed_PHO", "TSC2_S1452")

# plot WNT pathway for UCEC -----------------------------------------------
# mut_genes <- SMGs[["UCEC"]]
# cna_genes <- c("LRP1B",
#                "FGFR3", "FGFR1", "ERBB2", "IGF1R", "SOX17", "CCNE1", "MYC")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("GSK3B", "MAP4K5")
# pho_genes <- c("CTNNB1", "APC")
# rsds <- c("S552", "S2392")
# row_order <- c("GSK3B_collapsed_PHO", "CTNNB1_S552",
#                "MAP4K5_collapsed_PHO", "APC_S2392")

# plot MTOR pathway for UCEC -----------------------------------------------
# mut_genes <- SMGs[["UCEC"]]
# cna_genes <- c("LRP1B",
#                "FGFR3", "FGFR1", "ERBB2", "IGF1R", "SOX17", "FGFR3", "CCNE1", "MYC")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("EGFR", "STK11")
# pho_genes <- c("SHC1", "PTEN")
# rsds <- c("Y428", "S467")
# row_order <- c("EGFR_collapsed_PHO", "SHC1_Y428",
#                "STK11_collapsed_PHO", "PTEN_S467")

# plot 198753\tERK/MAPK targets pathway for COAD -----------------------------------------------
# mut_genes <- SMGs[["CO"]]
# cna_genes <- c("SMAD4", "TP53", "FHIT", "RBFOX1", "WWOX", "APC", "PTEN", "SMAD3", "TCF7L2", "USP12",
#                "KLF5", "WHSC1L1", "MYC", "ERBB2", "INS", "IGF2", "TH")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("PRKACA",  "PRKD1")
# pho_genes <- c("CTNNB1")
# rsds <- c("S191")
# row_order <- c("PRKACA_collapsed_PHO", "CTNNB1_S191",
#                "PRKD1_collapsed_PHO",
#                "PAK4_collapsed_PHO")

# # plot 198753\tERK/MAPK targets pathway for COAD -----------------------------------------------
# mut_genes <- SMGs[["CO"]]
# cna_genes <- c("SMAD4", "TP53", "FHIT", "RBFOX1", "WWOX", "APC", "PTEN", "SMAD3", "TCF7L2", "USP12",
#                "KLF5", "WHSC1L1", "MYC", "ERBB2", "INS", "IGF2", "TH")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("MAPK1")
# pho_genes <- c("RPS6KA1")
# rsds <- c("T582")
# row_order <- c("MAPK1_collapsed_PHO", "RPS6KA1_T582")

# # plot 177929\tSignaling by EGFR pathway for COAD -----------------------------------------------
# mut_genes <- SMGs[["CO"]]
# cna_genes <- c("SMAD4", "TP53", "FHIT", "RBFOX1", "WWOX", "APC", "PTEN", "SMAD3", "TCF7L2", "USP12",
#                "KLF5", "WHSC1L1", "MYC", "ERBB2", "INS", "IGF2", "TH")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("SRC", "RAF1", "MAPK1")
# pho_genes <- c("RAF1", "MAP2K2", "EGFR")
# rsds <- c("S43", "T396", "S946")
# row_order <- c("SRC_collapsed_PHO", "RAF1_S43",
#                "RAF1_collapsed_PHO", "MAP2K2_T396",
#                "MAPK1_collapsed_PHO", "EGFR_S946")

# plot 177929\tSignaling by EGFR pathway for COAD -----------------------------------------------
# mut_genes <- SMGs[["CO"]]
# cna_genes <- c("SMAD4", "TP53", "FHIT", "RBFOX1", "WWOX", "APC", "PTEN", "SMAD3", "TCF7L2", "USP12",
#                "KLF5", "WHSC1L1", "MYC", "ERBB2", "INS", "IGF2", "TH")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("PTPN2", "MAPK1")
# pho_genes <- c("EGFR", "EGFR")
# rsds <- c("S994", "S946")
# row_order <- c("PTPN2_collapsed_PHO", "EGFR_S994",
#                "MAPK1_collapsed_PHO", "EGFR_S946")

# plot 167044        Signalling to RAS pathway for COAD -----------------------------------------------
# mut_genes <- SMGs[["CO"]]
# cna_genes <- c("SMAD4", "TP53", "FHIT", "RBFOX1", "WWOX", "APC", "PTEN", "SMAD3", "TCF7L2", "USP12",
#                "KLF5", "WHSC1L1", "MYC", "ERBB2", "INS", "IGF2", "TH")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("SRC", "RAF1", "MAP2K6")
# pho_genes <- c("ARAF", "MAP2K2", "MAPK14")
# rsds <- c("S157", "T396")
# row_order <- c("SRC_collapsed_PHO", "ARAF_S157",
#                "RAF1_collapsed_PHO", "MAP2K2_T396",
#                "MAP2K6_collapsed_PHO", "MAPK14_T180Y182")

# plot 167044        Signalling to RAS pathway for COAD -----------------------------------------------
# mut_genes <- SMGs[["CO"]]
# cna_genes <- c("SMAD4", "TP53", "FHIT", "RBFOX1", "WWOX", "APC", "PTEN", "SMAD3", "TCF7L2", "USP12",
#                "KLF5", "WHSC1L1", "MYC", "ERBB2", "INS", "IGF2", "TH")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("MAPK3", "RPS6KA5", "LYN", "PTPN2")
# pho_genes <- c("STAT3", "STAT3", "STAT5A", "STAT3")
# rsds <- c("S726", "S726", "Y694", "S726")
# row_order <- c("MAPK3_collapsed_PHO", "STAT3_S726",
#                "RPS6KA5_collapsed_PHO", "STAT3_S726",
#                "PTPN2_collapsed_PHO", "STAT3_S726")

# plot 109704        PI3K Cascade pathway for COAD -----------------------------------------------
# mut_genes <- SMGs[["CO"]]
# cna_genes <- c("SMAD4", "TP53", "FHIT", "RBFOX1", "WWOX", "APC", "PTEN", "SMAD3", "TCF7L2", "USP12",
#                "KLF5", "WHSC1L1", "MYC", "ERBB2", "INS", "IGF2", "TH")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("AKT1")
# pho_genes <- c("AKT1S1")
# rsds <- c("S203")
# row_order <- c("AKT1_collapsed_PHO", "AKT1S1_S203")

# # plot 177929\tSignaling by EGFR pathway for COAD -----------------------------------------------
# mut_genes <- SMGs[["CO"]]
# cna_genes <- c("SMAD4", "TP53", "FHIT", "RBFOX1", "WWOX", "APC", "PTEN", "SMAD3", "TCF7L2", "USP12",
#                "KLF5", "WHSC1L1", "MYC", "ERBB2", "INS", "IGF2", "TH")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("PTPN2", "MAPK3", "LYN")
# pho_genes <- c("STAT3", "STAT3", "STAT5A")
# rsds <- c("S726", "S726", "Y694")
# row_order <- c("PTPN2_collapsed_PHO", "STAT3_S726",
#                "MAPK3_collapsed_PHO", "STAT3_T396",
#                "LYN_collapsed_PHO", "STAT5A_Y694")

# # plot mTOR pathway for COAD -----------------------------------------------
# mut_genes <- SMGs[["CO"]]
# cna_genes <- c("SMAD4", "TP53", "FHIT", "RBFOX1", "WWOX", "APC", "PTEN", "SMAD3", "TCF7L2", "USP12",
#                "KLF5", "WHSC1L1", "MYC", "ERBB2", "INS", "IGF2", "TH")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("PRKAB1")
# pho_genes <- c("MTOR")
# rsds <- c("S1261")
# row_order <- c("PRKAB1_collapsed_PHO", "MTOR_S1261")

# # plot metabolism pathway for COAD -----------------------------------------------
# mut_genes <- SMGs[["CO"]]
# cna_genes <- c("SMAD4", "TP53", "FHIT", "RBFOX1", "WWOX", "APC", "PTEN", "SMAD3", "TCF7L2", "USP12",
#                "KLF5", "WHSC1L1", "MYC", "ERBB2", "INS", "IGF2", "TH")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("PRKACA", "AKT1", "PRKCA", "PPP1CA")
# pho_genes <- c("ACACA", "NCOR1", "STXBP1", "PRKD1")
# rsds <- c("S66", "S2151", "S594", "S205")
# row_order <- c("PRKACA_collapsed_PHO", "ACACA_S66",
#                "AKT1_collapsed_PHO", "NCOR1_S2151",
#                "PRKCA_collapsed_PHO", "STXBP1_S594",
#                "PPP1CA_collapsed_PHO", "PRKD1_S205")

# plot 157118        Signaling by NOTCH pathway for OV -----------------------------------------------
# mut_genes <- SMGs[["OV"]]
# cna_genes <- c("LRP1B", "PRIM2", "CDKN2A", "PTEN", "RB1", "CREBBP", "WWOX", "ANKRD11", "MAP2K4", "NF1",
#                "AURKAIP1", "MYCL1", "MCL1", "XPR1", "PPP1CB", "PAX8", "CD47", "MECOM", "TACC3", "ANKRD17", "TERT", "SKP2", "ID4", "IKBKB", "SOX17", "DEPTOR", "MYC", "ALG8", "SC5DL", "KRAS", "ERBB3", "FGS2", "METTL17", "AKT1", "ERBB2", "HOXB", "TAF4B", "CCNE1", "BCL2L1", "ZMYND8", "MTMR3")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("GSK3B")
# pho_genes <- c("NOTCH2")
# rsds <- c("S2115")
# row_order <- c("GSK3B_collapsed_PHO", "NOTCH2_S2115")

# plot 389357        CD28 dependent PI3K/Akt signaling pathway for OV -----------------------------------------------
# mut_genes <- SMGs[["OV"]]
# cna_genes <- c("LRP1B", "PRIM2", "CDKN2A", "PTEN", "RB1", "CREBBP", "WWOX", "ANKRD11", "MAP2K4", "NF1",
#                "AURKAIP1", "MYCL1", "MCL1", "XPR1", "PPP1CB", "PAX8", "CD47", "MECOM", "TACC3", "ANKRD17", "TERT", "SKP2", "ID4", "IKBKB", "SOX17", "DEPTOR", "MYC", "ALG8", "SC5DL", "KRAS", "ERBB3", "FGS2", "METTL17", "AKT1", "ERBB2", "HOXB", "TAF4B", "CCNE1", "BCL2L1", "ZMYND8", "MTMR3")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("AKT1")
# pho_genes <- c("EIF4EBP1")
# rsds <- c("T70")
# row_order <- c("AKT1_collapsed_PHO", "EIF4EBP1_T70")

# plot 389357        CD28 dependent PI3K/Akt signaling pathway for OV -----------------------------------------------
# mut_genes <- SMGs[["OV"]]
# cna_genes <- c("LRP1B", "PRIM2", "CDKN2A", "PTEN", "RB1", "CREBBP", "WWOX", "ANKRD11", "MAP2K4", "NF1",
#                "AURKAIP1", "MYCL1", "MCL1", "XPR1", "PPP1CB", "PAX8", "CD47", "MECOM", "TACC3", "ANKRD17", "TERT", "SKP2", "ID4", "IKBKB", "SOX17", "DEPTOR", "MYC", "ALG8", "SC5DL", "KRAS", "ERBB3", "FGS2", "METTL17", "AKT1", "ERBB2", "HOXB", "TAF4B", "CCNE1", "BCL2L1", "ZMYND8", "MTMR3")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("AKT1")
# pho_genes <- c("EIF4EBP1")
# rsds <- c("T70")
# row_order <- c("AKT1_collapsed_PHO", "EIF4EBP1_T70")

# plot 69563        p53-Dependent G1 DNA Damage Response pathway for OV -----------------------------------------------
# mut_genes <- SMGs[["OV"]]
# cna_genes <- c("LRP1B", "PRIM2", "CDKN2A", "PTEN", "RB1", "CREBBP", "WWOX", "ANKRD11", "MAP2K4", "NF1",
#                "AURKAIP1", "MYCL1", "MCL1", "XPR1", "PPP1CB", "PAX8", "CD47", "MECOM", "TACC3", "ANKRD17", "TERT", "SKP2", "ID4", "IKBKB", "SOX17", "DEPTOR", "MYC", "ALG8", "SC5DL", "KRAS", "ERBB3", "FGS2", "METTL17", "AKT1", "ERBB2", "HOXB", "TAF4B", "CCNE1", "BCL2L1", "ZMYND8", "MTMR3")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("CDK2")
# pho_genes <- c("RB1", "TP53BP1")
# rsds <- c("T826", "S1367")
# row_order <- c("CDK2_collapsed_PHO", "RB1_T826", "TP53BP1_S1367")

# plot 453279\tMitotic G1-G1/S phases pathway for OV -----------------------------------------------
# mut_genes <- SMGs[["OV"]]
# cna_genes <- c("LRP1B", "PRIM2", "CDKN2A", "PTEN", "RB1", "CREBBP", "WWOX", "ANKRD11", "MAP2K4", "NF1",
#                "AURKAIP1", "MYCL1", "MCL1", "XPR1", "PPP1CB", "PAX8", "CD47", "MECOM", "TACC3", "ANKRD17", "TERT", "SKP2", "ID4", "IKBKB", "SOX17", "DEPTOR", "MYC", "ALG8", "SC5DL", "KRAS", "ERBB3", "FGS2", "METTL17", "AKT1", "ERBB2", "HOXB", "TAF4B", "CCNE1", "BCL2L1", "ZMYND8", "MTMR3")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("CDK1")
# pho_genes <- c("RB1")
# rsds <- c("T356")
# row_order <- c("CDK1_collapsed_PHO", "RB1_T356")


# # Mitotic G2-G2/M phases for OV--------------------------------------------------
# mut_genes <- SMGs[["OV"]]
# cna_genes <- c("LRP1B", "PRIM2", "CDKN2A", "PTEN", "RB1", "CREBBP", "WWOX", "ANKRD11", "MAP2K4", "NF1",
#                "AURKAIP1", "MYCL1", "MCL1", "XPR1", "PPP1CB", "PAX8", "CD47", "MECOM", "TACC3", "ANKRD17", "TERT", "SKP2", "ID4", "IKBKB", "SOX17", "DEPTOR", "MYC", "ALG8", "SC5DL", "KRAS", "ERBB3", "FGS2", "METTL17", "AKT1", "ERBB2", "HOXB", "TAF4B", "CCNE1", "BCL2L1", "ZMYND8", "MTMR3")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("CDK2", "PRKACA", "CDK1")
# pho_genes <- c("CENPF", "HSP90AA1", "NUMA1")
# rsds <- c("S3023", "S598", "S2063")
# row_order <- c("CDK2_collapsed_PHO", "CENPF_S3023",
#                "PRKACA_collapsed_PHO", "HSP90AA1_S598",
#                "CDK1_collapsed_PHO", "NUMA1_S2063")

# # 187687\tSignalling to ERKs for OV--------------------------------------------------
# mut_genes <- SMGs[["OV"]]
# cna_genes <- c("LRP1B", "PRIM2", "CDKN2A", "PTEN", "RB1", "CREBBP", "WWOX", "ANKRD11", "MAP2K4", "NF1",
#                "AURKAIP1", "MYCL1", "MCL1", "XPR1", "PPP1CB", "PAX8", "CD47", "MECOM", "TACC3", "ANKRD17", "TERT", "SKP2", "ID4", "IKBKB", "SOX17", "DEPTOR", "MYC", "ALG8", "SC5DL", "KRAS", "ERBB3", "FGS2", "METTL17", "AKT1", "ERBB2", "HOXB", "TAF4B", "CCNE1", "BCL2L1", "ZMYND8", "MTMR3")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("SRC", "MAPK3")
# pho_genes <- c("CDK1", "RAF1", "SHC1", "BRAF", "RAF1", "SHC1")
# rsds <- c("Y15", "S29", "Y428", "S365", "S621", "Y428")
# row_order <- c("SRC_collapsed_PHO", "CDK1_Y15", "RAF1_S29", "SHC1_Y428", 
#                "MAPK3_collapsed_PHO", "BRAF_S365", "RAF1_S621", "SHC1_Y428")


# # 449147\tSignaling by Interleukins for OV--------------------------------------------------
# mut_genes <- SMGs[["OV"]]
# cna_genes <- c("LRP1B", "PRIM2", "CDKN2A", "PTEN", "RB1", "CREBBP", "WWOX", "ANKRD11", "MAP2K4", "NF1",
#                "AURKAIP1", "MYCL1", "MCL1", "XPR1", "PPP1CB", "PAX8", "CD47", "MECOM", "TACC3", "ANKRD17", "TERT", "SKP2", "ID4", "IKBKB", "SOX17", "DEPTOR", "MYC", "ALG8", "SC5DL", "KRAS", "ERBB3", "FGS2", "METTL17", "AKT1", "ERBB2", "HOXB", "TAF4B", "CCNE1", "BCL2L1", "ZMYND8", "MTMR3")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("IKBKG")
# pho_genes <- c("NFKB1")
# rsds <- c("S907")
# row_order <- c("IKBKG_collapsed_PHO", "NFKB1_S907")

# # 2219528\tPI3K/AKT Signaling in Cancer for OV--------------------------------------------------
# mut_genes <- SMGs[["OV"]]
# cna_genes <- c("LRP1B", "PRIM2", "CDKN2A", "PTEN", "RB1", "CREBBP", "WWOX", "ANKRD11", "MAP2K4", "NF1",
#                "AURKAIP1", "MYCL1", "MCL1", "XPR1", "PPP1CB", "PAX8", "CD47", "MECOM", "TACC3", "ANKRD17", "TERT", "SKP2", "ID4", "IKBKB", "SOX17", "DEPTOR", "MYC", "ALG8", "SC5DL", "KRAS", "ERBB3", "FGS2", "METTL17", "AKT1", "ERBB2", "HOXB", "TAF4B", "CCNE1", "BCL2L1", "ZMYND8", "MTMR3")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("FYN", "AKT1", "EGFR")
# pho_genes <- c("PTPN11", "IRS1", "GAB1")
# rsds <- c("Y584", "S330", "S713")
# row_order <- c("FYN_collapsed_PHO", "PTPN11_Y584", 
#                "AKT1_collapsed_PHO", "IRS1_S330",
#                "EGFR_collapsed_PHO", "GAB1_S713")

# # 177929        Signaling by EGFR for OV--------------------------------------------------
# mut_genes <- SMGs[["OV"]]
# cna_genes <- c("LRP1B", "PRIM2", "CDKN2A", "PTEN", "RB1", "CREBBP", "WWOX", "ANKRD11", "MAP2K4", "NF1",
#                "AURKAIP1", "MYCL1", "MCL1", "XPR1", "PPP1CB", "PAX8", "CD47", "MECOM", "TACC3", "ANKRD17", "TERT", "SKP2", "ID4", "IKBKB", "SOX17", "DEPTOR", "MYC", "ALG8", "SC5DL", "KRAS", "ERBB3", "FGS2", "METTL17", "AKT1", "ERBB2", "HOXB", "TAF4B", "CCNE1", "BCL2L1", "ZMYND8", "MTMR3")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("FYN", "AKT1", "EGFR")
# pho_genes <- c("PTPN11", "IRS1", "GAB1")
# rsds <- c("Y584", "S330", "S713")
# row_order <- c("FYN_collapsed_PHO", "PTPN11_Y584", 
#                "AKT1", "IRS1_S330",
#                "EGFR", "GAB1_S713")

# 1250347        SHC1 events in ERBB4 signaling for BRCA--------------------------------------------------
# mut_genes <- SMGs[["BRCA"]]
# cna_genes <- c("MLL3", "PTEN", "RB1", "MAP2K4",
#                "PIK3CA", "EGFR", "FOXA1", "ERBB2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("EGFR")
# pho_genes <- c("GAB1")
# rsds <- c("S581")
# row_order <- c("EGFR_collapsed_PHO", "GAB1_S581")

# # 450282\tMAPK targets/ Nuclear events mediated by MAP kinases for BRCA--------------------------------------------------
# mut_genes <- SMGs[["BRCA"]]
# cna_genes <- c("MLL3", "PTEN", "RB1", "MAP2K4",
#                "PIK3CA", "EGFR", "FOXA1", "ERBB2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("MAPK1", "MAPK3", "MAPKAPK2")
# pho_genes <- c("RPS6KA1", "JUN", "CREB1")
# rsds <- c("S230", "S63", "S114")
# row_order <- c("MAPK1_collapsed_PHO", "RPS6KA1_S230", 
#                "MAPK3_collapsed_PHO", "JUN_S63",
#                "MAPKAPK2_collapsed_PHO", "CREB1_S114")

# # 179812\tGRB2 events in EGFR signaling for BRCA--------------------------------------------------
# mut_genes <- SMGs[["BRCA"]]
# cna_genes <- c("MLL3", "PTEN", "RB1", "MAP2K4",
#                "PIK3CA", "EGFR", "FOXA1", "ERBB2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("MAP2K1")
# pho_genes <- c("MAPK1")
# rsds <- c("Y187")
# row_order <- c("MAP2K1_collapsed_PHO", "MAPK1_Y187")

# # FOXM1 signaling for BRCA--------------------------------------------------
# mut_genes <- SMGs[["BRCA"]]
# cna_genes <- c("MLL3", "PTEN", "RB1", "MAP2K4",
#                "PIK3CA", "EGFR", "FOXA1", "ERBB2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("CDK1")
# pho_genes <- c("ANAPC1", "BUB1B", "CDC23", "PTTG1")
# rsds <- c("S51S60", "S676", "T596", "S165")
# row_order <- c("CDK1_collapsed_PHO", "ANAPC1_S51S60", "BUB1B_S676", "CDC23_T596", "PTTG1_S165")

# # MTOR signaling for BRCA--------------------------------------------------
# mut_genes <- SMGs[["BRCA"]]
# cna_genes <- c("MLL3", "PTEN", "RB1", "MAP2K4",
#                "PIK3CA", "EGFR", "FOXA1", "ERBB2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("PRKAB1", "PRKAB2")
# pho_genes <- c("MTOR", "RPTOR", "TSC1")
# rsds <- c("S1261", "S705", "T596")
# row_order <- c("PRKAB1_collapsed_PHO", "MTOR_S1261", 
#                "PRKAB2_collapsed_PHO", "RPTOR_S705", "TSC1_T596")

# 2219530\tConstitutive Signaling by Aberrant PI3K in Cancer for BRCA--------------------------------------------------
# mut_genes <- SMGs[["BRCA"]]
# cna_genes <- c("MLL3", "PTEN", "RB1", "MAP2K4",
#                "PIK3CA", "EGFR", "FOXA1", "ERBB2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("FYN", "PTPN11", "EGFR")
# pho_genes <- c("PTPN11", "PDGFRB", "GAB1")
# rsds <- c("Y584", "S705", "S581")
# row_order <- c("FYN_collapsed_PHO", "PTPN11_Y584",
#                "PTPN11_collapsed_PHO", "PDGFRB_S705",
#                "EGFR_collapsed_PHO", "GAB1_S581")

# 199418        Negative regulation of the PI3K/AKT network for BRCA--------------------------------------------------
# mut_genes <- SMGs[["BRCA"]]
# cna_genes <- c("MLL3", "PTEN", "RB1", "MAP2K4",
#                "PIK3CA", "EGFR", "FOXA1", "ERBB2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("AKT1")
# pho_genes <- c("GSK3B")
# rsds <- c("S9")
# row_order <- c("AKT1_collapsed_PHO", "GSK3B_S9")

# 450302        activated TAK1 mediates p38 MAPK activation network for BRCA--------------------------------------------------
# mut_genes <- SMGs[["BRCA"]]
# cna_genes <- c("MLL3", "PTEN", "RB1", "MAP2K4",
#                "PIK3CA", "EGFR", "FOXA1", "ERBB2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("MAP3K8", "MAP2K3")
# pho_genes <- c("MAP2K3", "MAPK14")
# rsds <- c("S189", "T180Y182")
# row_order <- c("MAP3K8_collapsed_PHO", "MAP2K3_S189",
#                "MAP2K3_collapsed_PHO", "MAPK14_T180Y182")

# 69541        Stabilization of p53 network for BRCA--------------------------------------------------
# mut_genes <- SMGs[["BRCA"]]
# cna_genes <- c("MLL3", "PTEN", "RB1", "MAP2K4",
#                "PIK3CA", "EGFR", "FOXA1", "ERBB2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("ATM")
# pho_genes <- c("TP53BP1", "RAD50", "PRKDC")
# rsds <- c("S1072", "S640", "S3205")
# row_order <- c("ATM_collapsed_PHO", "TP53BP1_S1072", "RAD50_S640", "PRKDC_S3205")

# 174143        APC/C-mediated degradation of cell cycle proteins network for BRCA--------------------------------------------------
# mut_genes <- SMGs[["BRCA"]]
# cna_genes <- c("MLL3", "PTEN", "RB1", "MAP2K4",
#                "PIK3CA", "EGFR", "FOXA1", "ERBB2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("PLK1", "CDK1")
# pho_genes <- c("MYC", "RB1")
# rsds <- c("S77", "S807")
# row_order <- c("PLK1_collapsed_PHO", "MYC_S77")

# 110049        MEK activation network for BRCA--------------------------------------------------
# mut_genes <- SMGs[["BRCA"]]
# cna_genes <- c("MLL3", "PTEN", "RB1", "MAP2K4",
#                "PIK3CA", "EGFR", "FOXA1", "ERBB2")
# rna_genes <- c("")
# pro_genes <- c("")
# phog_genes <- c("MAP2K1")
# pho_genes <- c("MAPK1")
# rsds <- c("Y187")
# row_order <- c("MAP2K1_collapsed_PHO", "MAPK1_Y187")

# PLK1-MYC for BRCA--------------------------------------------------
# mut_genes <- SMGs[["BRCA"]]
# cna_genes <- c("MLL3", "PTEN", "RB1", "MAP2K4",
#                "PIK3CA", "EGFR", "FOXA1", "ERBB2",
#                "MYCN")
# rna_genes <- c("")
# pro_genes <- c("MYC", "PLK1", "FBXW7")
# phog_genes <- c("PLK1")
# pho_genes <- c("MYC")
# rsds <- c("S77")
# row_order <- c("PLK1_collapsed_PHO", "MYC_S77", "MYC_PRO", "PLK1_PRO", "FBXW7_PRO")

# bussiness ------------------------------------------------------------------
geneA <- paste(unique(c(mut_genes, cna_genes)), collapse = "_")
geneB <- paste(unique(c(rna_genes, pro_genes, pho_genes)), collapse = "_")
phosphosite <- paste0(rsds, collapse = "_")


# for (cancer in c("UCEC", "BRCA", "CCRCC", "CO", "OV")) {
for (cancer in c("CCRCC")) {
  subdir1 <- paste0(makeOutDir(resultD = resultD), cancer, "/")
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
    mat_value <- mat_value[intersect(row_order, rownames(mat_value)),]
    fn <- paste0(subdir1, paste(rownames(mat_value), collapse = "_"), "_", cancer, "_withlegend.pdf")
    
    # plotting ----------------------------------------------------------------
    
    ## reformat the outlier status column
    if ("esscore_outlier" %in% colnames(col_anno)) {
      tmp <- vector(mode = "character", length = nrow(col_anno))
      tmp[is.na(col_anno$esscore_outlier)] <- NA
      tmp[!is.na(col_anno$esscore_outlier) & col_anno$esscore_outlier == TRUE] <- "TRUE"
      tmp[!is.na(col_anno$esscore_outlier) & col_anno$esscore_outlier == FALSE] <- "FALSE"
      col_anno$esscore_outlier <- tmp
    }
    
    plotpheatmap(mat_value, color.palette, col_anno, ann_colors, 20, 9, fn = fn)
    plotpheatmap(mat_value, color.palette, col_anno, ann_colors, 20, 7, fn = fn)
    
  }
}


