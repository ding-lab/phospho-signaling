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
reg_nonNA <- 20
fdr_thres <- c(0.05, 0.05); names(fdr_thres) <- c("phosphatase", "kinase")

# inputs ------------------------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                   "regression", "_", "cptac2p_cptac3", "_", "tumor", "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
regression$direct <- (!is.na(regression$Source)) & !(regression$Source %in% c("PhosphoNetworks", "MIMP", "NetKIN"))
# Input phosphatase -------------------------------------------------------
enzyme_type <- "phosphatase"
reg_nonNA <- reg_nonNAs[enzyme_type]
tab_pp <- regression[regression$enzyme_type == enzyme_type & regression$SELF == SELF,]
tab_pp <- markSigKS(regression = tab_pp, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)
tab_pp$regulated <- (tab_pp$fdr_sig & tab_pp$coef_sig)
tab_pp_reg <- tab_pp[tab_pp$regulated,] 

# Input kinase -------------------------------------------------------
enzyme_type <- "kinase"
reg_nonNA <- reg_nonNAs[enzyme_type]
tab_pk <- regression[regression$enzyme_type == enzyme_type & regression$SELF == SELF,]
tab_pk <- markSigKS(regression = tab_pk, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)
tab_pk$regulated <- (tab_pk$fdr_sig & tab_pk$coef_sig)
tab_pk_reg <- tab_pk[tab_pk$regulated,] 


# overlap phosphatase and kinase ------------------------------------------
tab_pk_pp <- merge(tab_pk_reg[tab_pk_reg$SELF == "trans", c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "pair", "Source")], 
                   tab_pp_reg[tab_pp_reg$SELF == "trans", c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "pair", "Source")], 
                   by = c("SUB_GENE", "SUB_MOD_RSD", "Cancer"), suffixes = c(".pk", ".pp"))

# plot the substrates in direct target -------------------------------
## filtering
tab_pk_pp_annotated <- tab_pk_pp[(!is.na(tab_pk_pp$Source.pk) & !(tab_pk_pp$Source.pk %in% c("PhosphoNetworks", "MIMP", "NetKIN"))) | !is.na(tab_pk_pp$Source.pp),]
overlap_pk_sub <- unique(tab_pk_pp_annotated[,c("Cancer", "pair.pk")])
overlap_pp_sub <- unique(tab_pk_pp_annotated[,c("Cancer", "pair.pp")])

tab_pk_reg2p <- merge(overlap_pk_sub, 
                      tab_pk_reg,
                      by.x = c("Cancer", "pair.pk"), by.y = c("Cancer", "pair"), all.x = T); tab_pk_reg2p$pair.pk <- NULL
tab_pp_reg2p <- merge(overlap_pp_sub, 
                      tab_pp_reg,
                      by.x = c("Cancer", "pair.pp"), by.y = c("Cancer", "pair"), all.x = T); tab_pp_reg2p$pair.pp <- NULL

tab2p <- rbind(tab_pk_reg2p, tab_pp_reg2p)
tab2p$site <- paste0(tab2p$SUB_GENE, "_", tab2p$SUB_MOD_RSD)
tab2p$log10_FDR <- -log10(tab2p$FDR_pho_kin)
tab2p$Cancer <- order_cancer(tab2p$Cancer)



p <- ggplot()
p <- p + geom_point(data = tab2p, mapping = aes(x = site, y = coef_pho_kin, size = log10_FDR,  
                                                fill = Cancer, color = ifelse(!is.na(Source), "direct target", "covariation")), alpha = 0.6, shape = 21)
p <- p + scale_fill_manual(values = color_cancers2)
p <- p + scale_color_manual(values = c("direct target" = "black", "covariation" = "white"))
# p <- p + facet_grid(.~SUB_GENE, scales = "free_x", space = "free")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + geom_text_repel(data = tab2p, mapping = aes(x = site, y = coef_pho_kin, label = KINASE), 
                         size = 2, force = 2, direction = "both")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
p <- p + theme(panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 9, face = "bold", angle = 90))
p <- p + theme(panel.spacing.y = unit(0, "lines"))
p
fn = paste0(makeOutDir(resultD = resultD), 'kinase_phosphatase_overlap_substrates.pdf')
pdf(file = fn, height=5, width = 12, useDingbats = F)
print(p)
dev.off()


# plot the substrates in oncogenic pathways -------------------------------

sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(tab_pk_pp$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)

sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                 SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
tab_pk_pp <- merge(tab_pk_pp, sub_genes2pathways, all.x = T)

## filtering
tab_pk_pp_annotated <- tab_pk_pp[!is.na(tab_pk_pp$SUB_GENE.path),]
overlap_pk_sub <- unique(tab_pk_pp_annotated[,c("Cancer", "pair.pk")])
overlap_pp_sub <- unique(tab_pk_pp_annotated[,c("Cancer", "pair.pp")])

tab_pk_reg2p <- merge(overlap_pk_sub, 
                      tab_pk_reg,
                      by.x = c("Cancer", "pair.pk"), by.y = c("Cancer", "pair"), all.x = T); tab_pk_reg2p$pair.pk <- NULL
tab_pp_reg2p <- merge(overlap_pp_sub, 
                      tab_pp_reg,
                      by.x = c("Cancer", "pair.pp"), by.y = c("Cancer", "pair"), all.x = T); tab_pp_reg2p$pair.pp <- NULL

tab2p <- rbind(tab_pk_reg2p, tab_pp_reg2p)
tab2p$site <- paste0(tab2p$SUB_GENE, "_", tab2p$SUB_MOD_RSD)
tab2p$log10_FDR <- -log10(tab2p$FDR_pho_kin)
tab2p$Cancer <- order_cancer(tab2p$Cancer)
tab2p$top <- get_top_by_id(value_vector = tab2p$coef_pho_kin, tab2p$site, num_top = 1)
tab2p$bottom <- get_top_by_id(value_vector = -(tab2p$coef_pho_kin), tab2p$site, num_top = 1)

p <- ggplot()
p <- p + geom_point(data = tab2p, mapping = aes(x = site, y = coef_pho_kin, size = log10_FDR,  
                                                fill = Cancer, color = ifelse(tab2p$direct, "direct target", "covariation")), alpha = 0.6, shape = 21)
p <- p + scale_fill_manual(values = color_cancers2)
# p <- p + facet_grid(.~SUB_GENE, scales = "free_x", space = "free")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + geom_text_repel(data = tab2p[tab2p$top | tab2p$bottom | tab2p$direct,], mapping = aes(x = site, y = coef_pho_kin, label = KINASE, color = enzyme_type), 
                         size = 2, force = 2, direction = "both")
p <- p + scale_color_manual(values = c("direct target" = "black", "covariation" = "white", "kinase" = set1[1], "phosphatase" = set1[2]))
p <- p + xlab("substrate phosphosite") + ylab("correlation coefficient(beta)\nwith kinase/phosphatase\nphosphoprotein")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
p <- p + theme(panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 9, face = "bold", angle = 90))
p <- p + theme(panel.spacing.y = unit(0, "lines"))
p <- p + guides(color = F)
p
fn = paste0(makeOutDir(resultD = resultD), 'kinase_phosphatase_overlap_substrates_oncogenic_pathways.pdf')
pdf(file = fn, height=4, width = 10, useDingbats = F)
print(p)
dev.off()
