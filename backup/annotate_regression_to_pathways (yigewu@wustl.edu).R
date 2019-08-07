# Yige Wu @ WashU 2018 Jul
# barplot for the top kinases with high and low kinase-substrate pairs

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R')

library(dplyr)
library(ggrepel)
library(colorspace)
# variables ---------------------------------------------------------------
# reg_nonNA2test <- seq(from = 5, to = 30, by = 5)
reg_nonNAs <- c(15, 25); names(reg_nonNAs) <- c("phosphatase", "kinase")
fdr_thres <- c(0.2, 0.05); names(fdr_thres) <- c("phosphatase", "kinase")


# Input phosphatase -------------------------------------------------------
enzyme_type <- "phosphatase"
reg_nonNA <- reg_nonNAs[enzyme_type]
tab_pp <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                               enzyme_type, "_substrate_regression_cptac2p_3can_tumor",
                               "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
tab_pp <- markSigKS(regression = tab_pp, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)
tab_pp$regulated <- (tab_pp$fdr_sig & tab_pp$coef_sig)
tab_pp_reg <- tab_pp[tab_pp$regulated,] 

tab_pp_common <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                               enzyme_type, "_substrate_regression_cptac2p_3can_tumor",
                               "_reg_nonNA", reg_nonNA, "_commonsite_nonNA", reg_nonNA, ".txt"), data.table = F)
tab_pp_common <- markSigSiteCan(regression = tab_pp_common, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)

# Input kinase -------------------------------------------------------
enzyme_type <- "kinase"
reg_nonNA <- reg_nonNAs[enzyme_type]
tab_pk <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                               enzyme_type, "_substrate_regression_cptac2p_3can_tumor",
                               "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
tab_pk <- markSigKS(regression = tab_pk, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)
tab_pk$regulated <- (tab_pk$fdr_sig & tab_pk$coef_sig)
tab_pk_reg <- tab_pk[tab_pk$regulated,] 

tab_pk_common <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                      enzyme_type, "_substrate_regression_cptac2p_3can_tumor",
                                      "_reg_nonNA", reg_nonNA, "_commonsite_nonNA", reg_nonNA, ".txt"), data.table = F)
tab_pk_common <- markSigSiteCan(regression = tab_pk_common, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)

# Annotate substrates -----------------------------------------------------
sup_cans_tab_en_reg <- rbind(tab_pp_reg, tab_pk_reg)
tab_common <- rbind(tab_pp_common, tab_pk_common)

sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(sup_cans_tab_en_reg$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)
sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                 SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
sub_genes2pathways <- sub_genes2pathways[!(sub_genes2pathways$SUB_GENE == "APC" & sub_genes2pathways$SUB_GENE.path != "WNT") & !(sub_genes2pathways$SUB_GENE == "GSK3B" & sub_genes2pathways$SUB_GENE.path != "PI3K") & !(sub_genes2pathways$SUB_GENE == "IRS1" & sub_genes2pathways$SUB_GENE.path != "RTK RAS") & !(sub_genes2pathways$SUB_GENE == "AKT1" & sub_genes2pathways$SUB_GENE.path != "PI3K"),]

sup_cans_tab_en_reg <- merge(sup_cans_tab_en_reg, sub_genes2pathways, all.x = T)
sup_cans_tab_en_reg <- merge(sup_cans_tab_en_reg, tab_common[, c("pair", "Cancer", 
                                                                    "sig_BRCA", "sig_OV", "sig_CO", 
                                                                    "shared3can", "uniq_BRCA", "uniq_OV", "uniq_CO")], by = c("pair", "Cancer"), all.x = T)

# plot --------------------------------------------------------------------
tab_gene_as_en2p <- NULL
for (gene in c("AKT1", "ERBB2", "MAP2K4", "RB1", "CTNNB1", "PTEN", "PTPN22", "PTPRD")) {
  tab_gene_as_en <- sup_cans_tab_en_reg[sup_cans_tab_en_reg$KINASE == gene,]
  tab_gene_as_en_anno <- tab_gene_as_en[!is.na(tab_gene_as_en$SUB_GENE.path) | (!is.na(tab_gene_as_en$uniq_BRCA) & (tab_gene_as_en$uniq_BRCA | tab_gene_as_en$uniq_OV | tab_gene_as_en$uniq_CO)),]
  tab_gene_as_en_anno <- tab_gene_as_en_anno[!(duplicated(tab_gene_as_en_anno[, c("pair", "FDR_pho_kin")]))]
  tab_gene_as_en2p <- rbind(tab_gene_as_en2p, tab_gene_as_en_anno)
}

for (enzyme_type in c("phosphatase")) {
  enzyme_type <- "kinase"
  subdir1 <- paste0(makeOutDir(resultD = resultD), enzyme_type, "/")
  dir.create(subdir1)
  for (reg_nonNA in reg_nonNA2test) {
    reg_nonNA <- 20
    subdir2 <- paste0(subdir1, "reg_nonNA", reg_nonNA, "/")
    dir.create(subdir2)
    
    sup_cans_tab_en <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                            enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor",
                                            "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
    sup_cans_tab_en <- sup_cans_tab_en[order(sup_cans_tab_en$FDR_pho_kin),]
    for (fdr_thres in fdr_thres2process) {
      fdr_thres <- 0.05
      subdir3 <- paste0(subdir2, "fdr_thres", fdr_thres, "/")
      dir.create(subdir3)
      
      sup_cans_tab_en <- markSigKS(regression = sup_cans_tab_en, sig_thres = fdr_thres, enzyme_type = enzyme_type)
      sup_cans_tab_en$regulated <- (sup_cans_tab_en$fdr_sig & sup_cans_tab_en$coef_sig)
      sup_cans_tab_en_reg <- sup_cans_tab_en[sup_cans_tab_en$regulated,]
    }
  }
}

tab_pi3k_intra <- sup_cans_tab_en_reg[sup_cans_tab_en_reg$KINASE %in% pathway_genes[["PI3K"]] & sup_cans_tab_en_reg$SUBSTRATE %in% pathway_genes[["PI3K"]] & sup_cans_tab_en_reg$Cancer == "BRCA",]
tab_pi3k_downstream <- sup_cans_tab_en_reg[sup_cans_tab_en_reg$KINASE %in% pathway_genes[["PI3K"]] & sup_cans_tab_en_reg$Cancer == "BRCA" & sup_cans_tab_en_reg$SELF == "trans",]
tab_pi3k_upstream <- sup_cans_tab_en_reg[sup_cans_tab_en_reg$SUBSTRATE %in% pathway_genes[["PI3K"]] & sup_cans_tab_en_reg$Cancer == "BRCA" & sup_cans_tab_en_reg$SELF == "trans",]

sup_cans_tab_en_reg[sup_cans_tab_en_reg$KINASE %in% c("PIK3CA", "PTEN", "KRAS", "NRAS", "BRAF", "ACVR2A", "SMAD2", "SMAD4", "CTNNB1", "APC") & sup_cans_tab_en_reg$Cancer == "CO" & sup_cans_tab_en_reg$SELF == "trans",]
sup_cans_tab_en_reg[sup_cans_tab_en_reg$SUBSTRATE %in% genes2path & sup_cans_tab_en_reg$Cancer == "CO" & sup_cans_tab_en_reg$SELF == "trans",]

sup_cans_tab_en_reg[sup_cans_tab_en_reg$SUBSTRATE %in% c("KRAS") & sup_cans_tab_en_reg$Cancer == "BRCA" & sup_cans_tab_en_reg$SELF == "trans",]

