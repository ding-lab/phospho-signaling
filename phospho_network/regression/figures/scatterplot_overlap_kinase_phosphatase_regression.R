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
SELF <- "trans"

# input druggable kinases -------------------------------------------------
drug_genes <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/reference_files/gene_drug_list/Premed_raw_databases/drugBank/drug_list.txt", data.table = F, col.names = "gene")


# inputs ------------------------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                       "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                        data.table = F)
regression <- regression %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP"))] | pair_pro %in% psp_tab$pair_pro)

regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression <- annotate_ks_source(regression = regression)

# inputs regulatory sites------------------------------------------------------------------
regulatory_sites <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/Phospho_databases/PhosphositePlus/Mar_04_2019/Regulatory_sites", data.table = F)
regulatory_sites <- regulatory_sites %>%
  mutate(SUB_MOD_RSD = str_split_fixed(string = MOD_RSD, pattern = "-", n = 2)[,1])
regulatory_sites2merge <- regulatory_sites %>%
  select(GENE, SUB_MOD_RSD, ON_FUNCTION, ON_PROCESS, ON_PROT_INTERACT, ON_OTHER_INTERACT) %>%
  unique
regression <- merge(regression, regulatory_sites2merge, by.x = c("SUBSTRATE", "SUB_MOD_RSD"), by.y = c("GENE", "SUB_MOD_RSD"), all.x = T)
regression <- regression %>%
  mutate(is.functional = (!is.na(ON_FUNCTION) | !is.na(ON_PROCESS) | !is.na(ON_PROT_INTERACT) | !is.na(ON_OTHER_INTERACT)))

# Input phosphatase -------------------------------------------------------
enzyme_type <- "phosphatase"
tab_pp <- regression[regression$enzyme_type == enzyme_type & regression$SELF == SELF,]
tab_pp_reg <- tab_pp[tab_pp$regulated,] 

# Input kinase -------------------------------------------------------
enzyme_type <- "kinase"
tab_pk <- regression[regression$enzyme_type == enzyme_type & regression$SELF == SELF,]
tab_pk_reg <- tab_pk[tab_pk$regulated,] 

# overlap phosphatase and kinase ------------------------------------------
tab_pk_pp <- merge(tab_pk_reg[tab_pk_reg$SELF == "trans", c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "pair", "Source")], 
                   tab_pp_reg[tab_pp_reg$SELF == "trans", c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "pair", "Source")], 
                   by = c("SUB_GENE", "SUB_MOD_RSD", "Cancer"), suffixes = c(".pk", ".pp"))

# plot the substrates phosphosites with functional annotation -------------------------------
tab_pk_pp <- merge(tab_pk_reg[tab_pk_reg$is.functional, c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "pair", "Source")], 
                   tab_pp_reg[tab_pp_reg$is.functional, c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "pair", "Source")], 
                   by = c("SUB_GENE", "SUB_MOD_RSD", "Cancer"), suffixes = c(".pk", ".pp"))
# tab_pk_pp_annotated <- tab_pk_pp[(!is.na(tab_pk_pp$Source.pk) & !(tab_pk_pp$Source.pk %in% c("PhosphoNetworks", "MIMP", "NetKIN"))) | !is.na(tab_pk_pp$Source.pp),]
overlap_pk_sub <- unique(tab_pk_pp[,c("Cancer", "pair.pk")])
overlap_pp_sub <- unique(tab_pk_pp[,c("Cancer", "pair.pp")])

tab_pk_reg2p <- merge(overlap_pk_sub, 
                      unique(tab_pk_reg[,c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "pair", "coef_pho_kin", "FDR_pho_kin", "is.direct")]),
                      by.x = c("Cancer", "pair.pk"), by.y = c("Cancer", "pair"), all.x = T); tab_pk_reg2p$pair.pk <- NULL
tab_pp_reg2p <- merge(overlap_pp_sub, 
                      unique(tab_pp_reg[,c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "pair", "coef_pho_kin", "FDR_pho_kin", "is.direct")]),
                      by.x = c("Cancer", "pair.pp"), by.y = c("Cancer", "pair"), all.x = T); tab_pp_reg2p$pair.pp <- NULL

tab2p <- rbind(tab_pk_reg2p, tab_pp_reg2p)
tab2p$site <- paste0(tab2p$SUB_GENE, "_", tab2p$SUB_MOD_RSD)
tab2p$log10_FDR <- -log10(tab2p$FDR_pho_kin)
tab2p$Cancer <- order_cancer(tab2p$Cancer)

pos <- position_jitter(width = 0.1, seed = 1)
p <- ggplot()
p <- p + geom_point(data = tab2p, mapping = aes(x = SUB_MOD_RSD, y = coef_pho_kin, size = log10_FDR,  
                                                fill = Cancer, color = is.direct), alpha = 0.6, shape = 21, position = pos)
p <- p + scale_size_continuous(range = c(3,6))
p <- p + scale_fill_manual(values = color_5cancers)
p <- p + scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"))
p <- p + facet_grid(.~SUB_GENE, scales = "free_x", space = "free")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + geom_text_repel(data = tab2p, mapping = aes(x = SUB_MOD_RSD, y = coef_pho_kin, label = GENE), 
                         size = 3.5, force = 2, position = pos, direction = "both", segment.alpha = 0.5)
p <- p + xlab("substrate_phosphosite") + ylab("beta coeffcient for\nkinase/phosphatase phosphoprotein level")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, face = "bold", size = 10))
p <- p + theme(panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 12, face = "bold", angle = 90))
p <- p + theme(panel.spacing.y = unit(0, "lines"))
p <- p + theme(strip.background = element_rect(fill = "white", color = "white"))
p
fn = paste0(makeOutDir(resultD = resultD), 'kinase_phosphatase_overlap_functional_substrate_phosphosite.pdf')
pdf(file = fn, height=5, width = 12, useDingbats = F)
print(p)
dev.off()

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

pos <- position_jitter(width = 0.1, seed = 1)
p <- ggplot()
p <- p + geom_point(data = tab2p, mapping = aes(x = SUB_MOD_RSD, y = coef_pho_kin, size = log10_FDR,  
                                                fill = Cancer, color = is.direct), alpha = 0.6, shape = 21, position = pos)
p <- p + scale_size_continuous(range = c(3,6))
p <- p + scale_fill_manual(values = color_5cancers)
p <- p + scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"))
p <- p + facet_grid(.~SUB_GENE, scales = "free_x", space = "free")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + geom_text_repel(data = tab2p, mapping = aes(x = SUB_MOD_RSD, y = coef_pho_kin, label = KINASE), 
                         size = 3.5, force = 2, position = pos, direction = "both", segment.alpha = 0.5)
p <- p + xlab("substrate_phosphosite") + ylab("beta coeffcient for\nkinase/phosphatase phosphoprotein level")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, face = "bold", size = 10))
p <- p + theme(panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 12, face = "bold", angle = 90))
p <- p + theme(panel.spacing.y = unit(0, "lines"))
p <- p + theme(strip.background = element_rect(fill = "white", color = "white"))
p
fn = paste0(makeOutDir(resultD = resultD), 'kinase_phosphatase_overlap_substrates.pdf')
pdf(file = fn, height=5, width = 12, useDingbats = F)
print(p)
dev.off()

# plot the substrates in oncogenic pathways -------------------------------
tab_pk_pp <- merge(tab_pk_reg[tab_pk_reg$SELF == "trans", c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "pair", "Source")], 
                   tab_pp_reg[tab_pp_reg$SELF == "trans", c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "pair", "Source")], 
                   by = c("SUB_GENE", "SUB_MOD_RSD", "Cancer"), suffixes = c(".pk", ".pp"))
sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(tab_pk_pp$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)

sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                 SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
tab_pk_pp <- merge(tab_pk_pp, sub_genes2pathways, all.x = T)

## filtering
tab_pk_pp_annotated <- tab_pk_pp[!is.na(tab_pk_pp$SUB_GENE.path),]
overlap_pk_sub <- unique(tab_pk_pp_annotated[,c("Cancer", "pair.pk")])
overlap_pp_sub <- unique(tab_pk_pp_annotated[,c("Cancer", "pair.pp")])

tab_pk_reg2p <- merge(overlap_pk_sub, 
                      unique(tab_pk_reg[,c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "pair", "coef_pho_kin", "FDR_pho_kin", "is.direct")]),
                      by.x = c("Cancer", "pair.pk"), by.y = c("Cancer", "pair"), all.x = T); tab_pk_reg2p$pair.pk <- NULL
tab_pp_reg2p <- merge(overlap_pp_sub, 
                      unique(tab_pp_reg[,c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "pair", "coef_pho_kin", "FDR_pho_kin", "is.direct")]),
                      by.x = c("Cancer", "pair.pp"), by.y = c("Cancer", "pair"), all.x = T); tab_pp_reg2p$pair.pp <- NULL

tab2p <- rbind(tab_pk_reg2p, tab_pp_reg2p)
tab2p$site <- paste0(tab2p$SUB_GENE, "_", tab2p$SUB_MOD_RSD)
tab2p$log10_FDR <- -log10(tab2p$FDR_pho_kin)
tab2p$Cancer <- order_cancer(tab2p$Cancer)
tab2p$top <- get_top_by_id(value_vector = tab2p$coef_pho_kin, tab2p$site, num_top = 1)
tab2p$bottom <- get_top_by_id(value_vector = -(tab2p$coef_pho_kin), tab2p$site, num_top = 1)
tab2p$GENE_is.smg <- get_SMG_by_cancer(gene_vector = tab2p$GENE, cancer_vector = tab2p$Cancer)
tab2p_text_highlight <- tab2p[tab2p$GENE %in% drug_genes$gene,]
tab2p_text_blur <- tab2p[!(tab2p$GENE %in% drug_genes$gene),]


pos <- position_jitter(width = 0.1, seed = 1)
p <- ggplot()
p <- p + geom_point(data = tab2p, mapping = aes(x = SUB_MOD_RSD, y = coef_pho_kin, size = log10_FDR,  
                                                fill = Cancer, color = is.direct), alpha = 0.6, shape = 21, position = pos)
p <- p + scale_size_continuous(range = c(3,6))
p <- p + scale_fill_manual(values = color_5cancers)
p <- p + scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"))
p <- p + facet_grid(.~SUB_GENE, scales = "free_x", space = "free")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + geom_text_repel(data = tab2p_text_highlight, mapping = aes(x = SUB_MOD_RSD, y = coef_pho_kin, label = GENE), 
                         size = 3.5, force = 2, position = pos, direction = "both", segment.alpha = 0.5, angle = 90)
p <- p + geom_text_repel(data = tab2p_text_blur, mapping = aes(x = SUB_MOD_RSD, y = coef_pho_kin, label = GENE), 
                         size = 2, force = 2, position = pos, direction = "both", segment.alpha = 0.5, color = "grey50", angle = 90)
p <- p + xlab("substrate_phosphosite") + ylab("beta coeffcient for\nkinase/phosphatase phosphoprotein level")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, face = "bold", size = 10))
p <- p + theme(panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 12, face = "bold", angle = 90))
p <- p + theme(panel.spacing.y = unit(0, "lines"))
p <- p + theme(strip.background = element_rect(fill = "white", color = "white"))
p
fn = paste0(makeOutDir(resultD = resultD), 'kinase_phosphatase_overlap_substrates_oncogenic_pathways.pdf')
pdf(file = fn, height=6, width = 15, useDingbats = F)
print(p)
dev.off()


# spotcheck ---------------------------------------------------------------
regression %>%
  filter(GENE == "PTEN") %>%
  filter(regulated == T) %>%
  filter(SUB_GENE == "SHC1") %>%
  filter(SUB_MOD_RSD == "S139")

pho_tab <- loadParseProteomicsData(cancer = "UCEC", expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
pho_tab <- loadParseProteomicsData(cancer = "UCEC", expression_type  = "PHO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
pho_tab %>%
  filter(Gene == "SHC1")

regression %>%
  filter(GENE == "MET") %>%
  filter(regulated == T) %>%
  filter(SUB_GENE == "PTK2") %>%
  filter(Cancer == "CCRCC")

regression %>%
  filter(GENE == "PTPN11") %>%
  filter(regulated == T) %>%
  filter(SUB_GENE == "EGFR") %>%
  filter(Cancer == "CCRCC")

regression %>%
  filter(GENE == "MAPK14") %>%
  filter(regulated == T) %>%
  filter(SUB_GENE == "EGFR") %>%
  filter(Cancer == "CCRCC")
