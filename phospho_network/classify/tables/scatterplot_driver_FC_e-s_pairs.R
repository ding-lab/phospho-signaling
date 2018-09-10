# Yige Wu @ WashU 2018 Apr
# showing the fold change for enzyme substrate pairs involving drivers for 3 cancers


# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
library(ggrepel)

# inputs ------------------------------------------------------------------
## input druggable gene list
# drug_genes <- fread()
color_direction <- c(set1[1], set1[2], "black")
names(color_direction) <- c("up", "down", "NA")

# mutation drivers associated with consistent high/low pairs --------------
#for (cancer in c("BRCA")) {
for (cancer in c("BRCA", "OV", "CO")) {
  ## input super table
  ksea_diffexp <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_regression_plus_genoalt/", 
                                       cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  ksea_diffexp$pair <- paste0(ksea_diffexp$GENE, ":", ksea_diffexp$SUB_GENE, ":", ksea_diffexp$SUB_MOD_RSD)
  ## filter by driver genes
  driver_tab <- ksea_diffexp[(!is.na(ksea_diffexp$enzyme.driver_type) & ksea_diffexp$enzyme.driver_type %in% c("oncogene", "tsg")) | (!is.na(ksea_diffexp$substrate.driver_type) & ksea_diffexp$substrate.driver_type %in% c("oncogene", "tsg")),]
  driver_tab$driver_direction <- "NA"
  driver_tab$driver_direction[!is.na(driver_tab$enzyme.driver_type) & (driver_tab$enzyme.driver_type == "oncogene") & (driver_tab$enzyme_direction == "up")] <- "up"
  driver_tab$driver_direction[!is.na(driver_tab$enzyme.driver_type) & (driver_tab$enzyme.driver_type == "tsg") & (driver_tab$enzyme_direction == "down")] <- "down"
  driver_tab$driver_direction[!is.na(driver_tab$substrate.driver_type) & (driver_tab$substrate.driver_type == "oncogene") & (driver_tab$substrate_direction == "up")] <- "up"
  driver_tab$driver_direction[!is.na(driver_tab$substrate.driver_type) & (driver_tab$substrate.driver_type == "tsg") & (driver_tab$substrate_direction == "down")] <- "down"
  driver_tab$mut_direction <- ""
  driver_tab$pair_mut <- driver_tab$pair
  tmp <- (!is.na(driver_tab$Fold_Change) & (driver_tab$Fold_Change < 0) & (!is.na(driver_tab$adjusted_pvalue)) & (driver_tab$adjusted_pvalue < 0.1))
  driver_tab$mut_direction[tmp] <- "down"
  driver_tab$pair_mut[tmp] <- paste0(driver_tab$pair[tmp], "(mut-down)")
  tmp <- (!is.na(driver_tab$Fold_Change) & (driver_tab$Fold_Change > 0) & (!is.na(driver_tab$adjusted_pvalue)) & (driver_tab$adjusted_pvalue < 0.1))
  driver_tab$mut_direction[tmp] <- "up"
  driver_tab$pair_mut[tmp] <- paste0(driver_tab$pair[tmp], "(mut-up)")
  driver_cons_tab <- driver_tab[!is.na(driver_tab$consistent) & driver_tab$consistent,]; driver_cons_tab$pair_type <- "consistent"
  tmp <- data.frame(table(driver_cons_tab$GENE))
  outlier_genes <- as.vector(tmp$Var1[tmp$Freq > 50])
  
  p = ggplot(driver_cons_tab)
  p = p + geom_point(aes(x=KSEA_log2FC, y=substrate_log2FC, color = driver_direction), alpha=0.05, stroke = 0)
  p = p + geom_text_repel(aes(KSEA_log2FC, substrate_log2FC, label= ifelse(!(GENE %in% outlier_genes) & (driver_direction != "NA" | mut_direction != ""), as.character(pair_mut), NA), color = driver_direction, segment.color = driver_direction),
                          force = 1, 
                          segment.size = 0.5, segment.alpha = 0.2, 
                          size=1.5,alpha=0.8)
  p = p + geom_vline(xintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
  p = p + geom_hline(yintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
  p = p + scale_color_manual(values = color_direction)
  p = p + theme_nogrid()
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p = p + theme(title = element_text(size = 8),
                legend.key.size = unit(1, "cm"), legend.text = element_text(size = 5), 
                legend.background = element_rect(fill = NA), legend.box.background = element_rect(fill = NA),
                legend.key = element_rect(fill = NA), legend.key.height = unit(0.5, "cm"),
                legend.title = element_text(size = 5), legend.justification = c(1, 0), legend.position = c(1, 0))
  p
  resultDnow <- makeOutDir()
  fn = paste0(resultDnow, cancer, "_consistent_pairs_involving_drivers_mutimpact", "_fdr",sig,".pdf")
  ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
}


# mutation drivers associated with regulated pairs --------------
for (cancer in c("BRCA")) {
  ## input super table
  ksea_diffexp <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_regression_plus_genoalt/", 
                                       cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  ksea_diffexp$pair <- paste0(ksea_diffexp$GENE, ":", ksea_diffexp$SUB_GENE, ":", ksea_diffexp$SUB_MOD_RSD)
  ## filter by driver genes
  driver_tab <- ksea_diffexp[(!is.na(ksea_diffexp$enzyme.driver_type) & ksea_diffexp$enzyme.driver_type %in% c("oncogene", "tsg")) | (!is.na(ksea_diffexp$substrate.driver_type) & ksea_diffexp$substrate.driver_type %in% c("oncogene", "tsg")),]
  driver_tab$driver_direction <- "NA"
  driver_tab$driver_direction[!is.na(driver_tab$enzyme.driver_type) & (driver_tab$enzyme.driver_type == "oncogene") & (driver_tab$enzyme_direction == "up")] <- "up"
  driver_tab$driver_direction[!is.na(driver_tab$enzyme.driver_type) & (driver_tab$enzyme.driver_type == "tsg") & (driver_tab$enzyme_direction == "down")] <- "down"
  driver_tab$driver_direction[!is.na(driver_tab$substrate.driver_type) & (driver_tab$substrate.driver_type == "oncogene") & (driver_tab$substrate_direction == "up")] <- "up"
  driver_tab$driver_direction[!is.na(driver_tab$substrate.driver_type) & (driver_tab$substrate.driver_type == "tsg") & (driver_tab$substrate_direction == "down")] <- "down"
  driver_tab$mut_direction <- "NA"
  driver_tab$pair_mut <- driver_tab$pair
  tmp <- (!is.na(driver_tab$Fold_Change) & (driver_tab$Fold_Change < 0) & (!is.na(driver_tab$adjusted_pvalue)) & (driver_tab$adjusted_pvalue < 0.1))
  driver_tab$mut_direction[tmp] <- "down"
  driver_tab$pair_mut[tmp] <- paste0(driver_tab$pair[tmp], "(mut-down)")
  tmp <- (!is.na(driver_tab$Fold_Change) & (driver_tab$Fold_Change > 0) & (!is.na(driver_tab$adjusted_pvalue)) & (driver_tab$adjusted_pvalue < 0.1))
  driver_tab$mut_direction[tmp] <- "up"
  driver_tab$pair_mut[tmp] <- paste0(driver_tab$pair[tmp], "(mut-up)")
  
  ## get regulated pairs
  driver_reg_tab <- driver_tab[!is.na(driver_tab$FDR_pho_kin),]
  tmp <- data.frame(table(driver_reg_tab$GENE))
  outlier_genes <- as.vector(tmp$Var1[tmp$Freq > 50])
  
  cutoff <- 2
  p = ggplot(driver_reg_tab[!is.na(remove_outliers(x = driver_reg_tab$coef_pho_kin, out_thres = 1.5)),])
  p = p + geom_point(aes(x=coef_pho_kin, y=-log10(FDR_pho_kin), color = mut_direction), alpha=0.05, stroke = 0)
  p = p + geom_text_repel(aes(coef_pho_kin, -log10(FDR_pho_kin), label= ifelse(-log10(FDR_pho_kin) > cutoff | mut_direction != "NA", as.character(pair_mut), NA), color = mut_direction, segment.color = mut_direction),
                          force = 1, 
                          segment.size = 0.5, segment.alpha = 0.2, 
                          size=1.5,alpha=0.8)
  p = p + geom_vline(xintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
  p = p + geom_hline(yintercept = cutoff, color = 'grey', alpha = 0.5, linetype = 2)
  p = p + scale_color_manual(values = color_direction)
  p = p + theme_nogrid()
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p = p + theme(title = element_text(size = 8),
                legend.key.size = unit(1, "cm"), legend.text = element_text(size = 5), 
                legend.background = element_rect(fill = NA), legend.box.background = element_rect(fill = NA),
                legend.key = element_rect(fill = NA), legend.key.height = unit(0.5, "cm"),
                legend.title = element_text(size = 5), legend.justification = c(1, 0), legend.position = c(1, 0))
  p
  resultDnow <- makeOutDir()
  fn = paste0(resultDnow, cancer, "_regulated_pairs_involving_drivers_mutimpact", "_fdr",sig,".pdf")
  ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
}

for (cancer in c("BRCA")) {
  ## input super table
  ksea_diffexp <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_regression_plus_genoalt/", 
                                       cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  ksea_diffexp$pair <- paste0(ksea_diffexp$GENE, ":", ksea_diffexp$SUB_GENE, ":", ksea_diffexp$SUB_MOD_RSD)
  ## filter by driver genes
  driver_tab <- ksea_diffexp[(!is.na(ksea_diffexp$enzyme.driver_type) & ksea_diffexp$enzyme.driver_type %in% c("oncogene", "tsg")) | (!is.na(ksea_diffexp$substrate.driver_type) & ksea_diffexp$substrate.driver_type %in% c("oncogene", "tsg")),]
  driver_tab$driver_direction <- "NA"
  driver_tab$driver_direction[!is.na(driver_tab$enzyme.driver_type) & (driver_tab$enzyme.driver_type == "oncogene") & (driver_tab$enzyme_direction == "up")] <- "up"
  driver_tab$driver_direction[!is.na(driver_tab$enzyme.driver_type) & (driver_tab$enzyme.driver_type == "tsg") & (driver_tab$enzyme_direction == "down")] <- "down"
  driver_tab$driver_direction[!is.na(driver_tab$substrate.driver_type) & (driver_tab$substrate.driver_type == "oncogene") & (driver_tab$substrate_direction == "up")] <- "up"
  driver_tab$driver_direction[!is.na(driver_tab$substrate.driver_type) & (driver_tab$substrate.driver_type == "tsg") & (driver_tab$substrate_direction == "down")] <- "down"
  driver_tab$mut_direction <- "NA"
  driver_tab$pair_mut <- driver_tab$pair
  tmp <- (!is.na(driver_tab$Fold_Change) & (driver_tab$Fold_Change < 0) & (!is.na(driver_tab$adjusted_pvalue)) & (driver_tab$adjusted_pvalue < 0.1))
  driver_tab$mut_direction[tmp] <- "down"
  driver_tab$pair_mut[tmp] <- paste0(driver_tab$pair[tmp], "(mut-down)")
  tmp <- (!is.na(driver_tab$Fold_Change) & (driver_tab$Fold_Change > 0) & (!is.na(driver_tab$adjusted_pvalue)) & (driver_tab$adjusted_pvalue < 0.1))
  driver_tab$mut_direction[tmp] <- "up"
  driver_tab$pair_mut[tmp] <- paste0(driver_tab$pair[tmp], "(mut-up)")
  
  ## get regulated pairs
  driver_reg_tab <- driver_tab[!is.na(driver_tab$FDR_pho_kin) & ((!is.na(driver_tab$Mutated_Sample_Size) & driver_tab$Mutated_Sample_Size > 4) | (is.na(driver_tab$Mutated_Sample_Size))),]
  tmp <- data.frame(table(driver_reg_tab$GENE))
  outlier_genes <- as.vector(tmp$Var1[tmp$Freq > 50])
  
  cutoff <- 2
  p = ggplot(driver_reg_tab[!is.na(remove_outliers(x = driver_reg_tab$coef_pho_kin, out_thres = 1.5)),])
  p = p + geom_point(aes(x=coef_pho_kin, y=-log10(FDR_pho_kin), color = mut_direction), alpha=0.05, stroke = 0)
  p = p + geom_text_repel(aes(coef_pho_kin, -log10(FDR_pho_kin), label= ifelse((-log10(FDR_pho_kin) > cutoff | mut_direction != "NA"), as.character(pair_mut), NA), color = mut_direction, segment.color = mut_direction),
                          force = 1, 
                          segment.size = 0.5, segment.alpha = 0.2, 
                          size=1.5,alpha=0.8)
  p = p + geom_vline(xintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
  p = p + geom_hline(yintercept = cutoff, color = 'grey', alpha = 0.5, linetype = 2)
  p = p + scale_color_manual(values = color_direction)
  p = p + theme_nogrid()
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p = p + theme(title = element_text(size = 8),
                legend.key.size = unit(1, "cm"), legend.text = element_text(size = 5), 
                legend.background = element_rect(fill = NA), legend.box.background = element_rect(fill = NA),
                legend.key = element_rect(fill = NA), legend.key.height = unit(0.5, "cm"),
                legend.title = element_text(size = 5), legend.justification = c(1, 0), legend.position = c(1, 0))
  p
  resultDnow <- makeOutDir()
  fn = paste0(resultDnow, cancer, "_regulated_pairs_involving_drivers_mutimpact_noAKT1", "_fdr",sig,".pdf")
  ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
  
  cutoff <- 1
  p = ggplot(driver_reg_tab[!is.na(remove_outliers(x = driver_reg_tab$coef_pho_kin, out_thres = 1.5)),])
  p = p + geom_point(aes(x=coef_pho_kin, y=-log10(P_pho_kin), color = mut_direction), alpha=0.05, stroke = 0)
  p = p + geom_text_repel(aes(coef_pho_kin, -log10(P_pho_kin), label= ifelse((-log10(P_pho_kin) > cutoff & mut_direction != "NA"), as.character(pair_mut), NA), color = mut_direction, segment.color = mut_direction),
                          force = 1, 
                          segment.size = 0.5, segment.alpha = 0.2, 
                          size=1.5,alpha=0.8)
  p = p + geom_vline(xintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
  p = p + geom_hline(yintercept = cutoff, color = 'grey', alpha = 0.5, linetype = 2)
  p = p + scale_color_manual(values = color_direction)
  p = p + theme_nogrid()
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p = p + theme(title = element_text(size = 8),
                legend.key.size = unit(1, "cm"), legend.text = element_text(size = 5), 
                legend.background = element_rect(fill = NA), legend.box.background = element_rect(fill = NA),
                legend.key = element_rect(fill = NA), legend.key.height = unit(0.5, "cm"),
                legend.title = element_text(size = 5), legend.justification = c(1, 0), legend.position = c(1, 0))
  p
  resultDnow <- makeOutDir()
  fn = paste0(resultDnow, cancer, "_regulated_pairs_involving_drivers_mutimpact_noAKT1", "_pvalue",sig,".pdf")
  ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
}


for (cancer in c( "OV", "CO")) {
  ## input super table
  ksea_diffexp <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_regression_plus_genoalt/", 
                                       cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  ksea_diffexp$pair <- paste0(ksea_diffexp$GENE, ":", ksea_diffexp$SUB_GENE, ":", ksea_diffexp$SUB_MOD_RSD)
  ## filter by driver genes
  driver_tab <- ksea_diffexp[(!is.na(ksea_diffexp$enzyme.driver_type) & ksea_diffexp$enzyme.driver_type %in% c("oncogene", "tsg")) | (!is.na(ksea_diffexp$substrate.driver_type) & ksea_diffexp$substrate.driver_type %in% c("oncogene", "tsg")),]
  driver_tab$driver_direction <- "NA"
  driver_tab$driver_direction[!is.na(driver_tab$enzyme.driver_type) & (driver_tab$enzyme.driver_type == "oncogene") & (driver_tab$enzyme_direction == "up")] <- "up"
  driver_tab$driver_direction[!is.na(driver_tab$enzyme.driver_type) & (driver_tab$enzyme.driver_type == "tsg") & (driver_tab$enzyme_direction == "down")] <- "down"
  driver_tab$driver_direction[!is.na(driver_tab$substrate.driver_type) & (driver_tab$substrate.driver_type == "oncogene") & (driver_tab$substrate_direction == "up")] <- "up"
  driver_tab$driver_direction[!is.na(driver_tab$substrate.driver_type) & (driver_tab$substrate.driver_type == "tsg") & (driver_tab$substrate_direction == "down")] <- "down"
  driver_tab$mut_direction <- "NA"
  driver_tab$pair_mut <- driver_tab$pair
  tmp <- (!is.na(driver_tab$Fold_Change) & (driver_tab$Fold_Change < 0) & (!is.na(driver_tab$adjusted_pvalue)) & (driver_tab$adjusted_pvalue < 0.1))
  driver_tab$mut_direction[tmp] <- "down"
  driver_tab$pair_mut[tmp] <- paste0(driver_tab$pair[tmp], "(mut-down)")
  tmp <- (!is.na(driver_tab$Fold_Change) & (driver_tab$Fold_Change > 0) & (!is.na(driver_tab$adjusted_pvalue)) & (driver_tab$adjusted_pvalue < 0.1))
  driver_tab$mut_direction[tmp] <- "up"
  driver_tab$pair_mut[tmp] <- paste0(driver_tab$pair[tmp], "(mut-up)")
  
  ## get regulated pairs
  driver_reg_tab <- driver_tab[!is.na(driver_tab$FDR_pho_kin),]
  tmp <- data.frame(table(driver_reg_tab$GENE))
  outlier_genes <- as.vector(tmp$Var1[tmp$Freq > 50])
  
  cutoff <- 1
  p = ggplot(driver_reg_tab[!is.na(remove_outliers(x = driver_reg_tab$coef_pho_kin, out_thres = 1.5)),])
  p = p + geom_point(aes(x=coef_pho_kin, y=-log10(FDR_pho_kin), color = mut_direction), alpha=0.05, stroke = 0)
  p = p + geom_text_repel(aes(coef_pho_kin, -log10(FDR_pho_kin), label= ifelse(-log10(FDR_pho_kin) > cutoff | mut_direction != "NA", as.character(pair_mut), NA), color = mut_direction, segment.color = mut_direction),
                          force = 1, 
                          segment.size = 0.5, segment.alpha = 0.2, 
                          size=1.5,alpha=0.8)
  p = p + geom_vline(xintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
  p = p + geom_hline(yintercept = cutoff, color = 'grey', alpha = 0.5, linetype = 2)
  p = p + scale_color_manual(values = color_direction)
  p = p + theme_nogrid()
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p = p + theme(title = element_text(size = 8),
                legend.key.size = unit(1, "cm"), legend.text = element_text(size = 5), 
                legend.background = element_rect(fill = NA), legend.box.background = element_rect(fill = NA),
                legend.key = element_rect(fill = NA), legend.key.height = unit(0.5, "cm"),
                legend.title = element_text(size = 5), legend.justification = c(1, 0), legend.position = c(1, 0))
  p
  resultDnow <- makeOutDir()
  fn = paste0(resultDnow, cancer, "_regulated_pairs_involving_drivers_mutimpact", "_fdr",sig,".pdf")
  ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
  
  cutoff <- 1
  p = ggplot(driver_reg_tab[!is.na(remove_outliers(x = driver_reg_tab$coef_pho_kin, out_thres = 1.5)),])
  p = p + geom_point(aes(x=coef_pho_kin, y=-log10(P_pho_kin), color = mut_direction), alpha=0.05, stroke = 0)
  p = p + geom_text_repel(aes(coef_pho_kin, -log10(P_pho_kin), label= ifelse((-log10(P_pho_kin) > cutoff & mut_direction != "NA"), as.character(pair_mut), NA), color = mut_direction, segment.color = mut_direction),
                          force = 1, 
                          segment.size = 0.5, segment.alpha = 0.2, 
                          size=1.5,alpha=0.8)
  p = p + geom_vline(xintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
  p = p + geom_hline(yintercept = cutoff, color = 'grey', alpha = 0.5, linetype = 2)
  p = p + scale_color_manual(values = color_direction)
  p = p + theme_nogrid()
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p = p + theme(title = element_text(size = 8),
                legend.key.size = unit(1, "cm"), legend.text = element_text(size = 5), 
                legend.background = element_rect(fill = NA), legend.box.background = element_rect(fill = NA),
                legend.key = element_rect(fill = NA), legend.key.height = unit(0.5, "cm"),
                legend.title = element_text(size = 5), legend.justification = c(1, 0), legend.position = c(1, 0))
  p
  resultDnow <- makeOutDir()
  fn = paste0(resultDnow, cancer, "_regulated_pairs_involving_drivers_mutimpact", "_pvalue",sig,".pdf")
  ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
}

genes <- KEGG[[]]
reg_tab <- ksea_diffexp[!is.na(ksea_diffexp$FDR_pho_kin) & ksea_diffexp$FDR_pho_kin < 0.1,]
reg_tab <- ksea_diffexp
