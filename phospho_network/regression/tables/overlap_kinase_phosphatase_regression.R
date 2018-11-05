# Yige Wu @ WashU 2018 Aug
# pie chart showing proportion of regulated pairs by each substrate by each enzyme, also showing regulated ratio
# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R')
library(ggrepel)

# set variables -----------------------------------------------------------
top_kinase2show <- 10
top_substrate2show <- 2

reg_nonNAs <- c(15, 20); names(reg_nonNAs) <- c("phosphatase", "kinase")
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

# Input kinase -------------------------------------------------------
enzyme_type <- "kinase"
reg_nonNA <- reg_nonNAs[enzyme_type]
tab_pk <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                               enzyme_type, "_substrate_regression_cptac2p_3can_tumor",
                               "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
tab_pk <- markSigKS(regression = tab_pk, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)
tab_pk$regulated <- (tab_pk$fdr_sig & tab_pk$coef_sig)
tab_pk_reg <- tab_pk[tab_pk$regulated,] 


# overlap phosphatase and kinase ------------------------------------------
tab_pk_pp <- merge(tab_pk_reg[tab_pk_reg$SELF == "trans", c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "pair")], 
                   tab_pp_reg[tab_pp_reg$SELF == "trans", c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "pair")], 
                   by = c("SUB_GENE", "SUB_MOD_RSD", "Cancer"), suffixes = c(".pk", ".pp"))
sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(tab_pk_pp$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)

sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                 SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
tab_pk_pp <- merge(tab_pk_pp, sub_genes2pathways, all.x = T)
tab_pk_pp_annotated <- tab_pk_pp[!is.na(tab_pk_pp$SUB_GENE.path),]
overlap_pk_sub <- unique(tab_pk_pp_annotated[,c("Cancer", "pair.pk")])
overlap_pp_sub <- unique(tab_pk_pp_annotated[,c("Cancer", "pair.pp")])
tab_pk_reg2p <- merge(overlap_pk_sub, tab_pk_reg, by.x = c("Cancer", "pair.pk"), by.y = c("Cancer", "pair"), all.x = T); tab_pk_reg2p$pair.pk <- NULL
tab_pp_reg2p <- merge(overlap_pp_sub, tab_pp_reg, by.x = c("Cancer", "pair.pp"), by.y = c("Cancer", "pair"), all.x = T); tab_pp_reg2p$pair.pp <- NULL
tab2p <- rbind(tab_pk_reg2p, tab_pp_reg2p)
tab2p$site <- paste0(tab2p$SUB_GENE, "_", tab2p$SUB_MOD_RSD)
tab2p$log10_FDR <- -log10(tab2p$FDR_pho_kin)
tab2p$Cancer <- factor(tab2p$Cancer, levels = cancers_sort)

p <- ggplot()
p <- p + geom_point(data = tab2p, mapping = aes(x = site, y = coef_pho_kin, size = log10_FDR, color = Cancer), alpha = 0.6)
p <- p + scale_color_manual(values = color_cancers2)
p <- p + facet_grid(.~SUB_GENE, scales = "free_x", space = "free")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + geom_text_repel(data = tab2p, mapping = aes(x = site, y = coef_pho_kin, label = KINASE), 
                         size = 2, force = 2, direction = "both", angle = 90)
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
p <- p + theme(panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 9, face = "bold"))
p <- p + theme(panel.spacing.y = unit(0, "lines"))
p
fn = paste0(makeOutDir(resultD = resultD), 'kinase_phosphatase_overlap_substrates.pdf')
pdf(file = fn, height=5, width = 12, useDingbats = F)
print(p)
dev.off()

p <- ggplot()
p <- p + geom_point(data = tab2p[tab2p$Cancer %in% c("BRCA", "CO"),], mapping = aes(x = site, y = coef_pho_kin, size = log10_FDR, color = Cancer), alpha = 0.6)
p <- p + scale_color_manual(values = color_cancers2)
p <- p + facet_grid(.~SUB_GENE, scales = "free_x", space = "free")
p <- p + ylab("standardized coefficient of\nthe kinase/phosphatase phosphoprotein") + xlab("phosphorylation site")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + geom_text_repel(data = tab2p[tab2p$Cancer %in% c("BRCA", "CO"),], mapping = aes(x = site, y = coef_pho_kin, label = KINASE), 
                         size = 2, force = 2, direction = "both")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
p <- p + theme(panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 9, face = "bold"))
p <- p + theme(panel.spacing.y = unit(0, "lines"))
p
fn = paste0(makeOutDir(resultD = resultD), 'kinase_phosphatase_overlap_substrates_horizontal.pdf')
pdf(file = fn, height=5, width = 12, useDingbats = F)
print(p)
dev.off()
