# Yige Wu @ WashU 2018 Jan
# draw a grid showing the distribution of correlated kinase-substrate pairs are distributed across oncogenic pathways

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
setwd(baseD)
source('./cptac2p_analysis/preprocess_files/preprocess_files_shared.R')
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
library(pheatmap)

# set variables -----------------------------------------------------------
reg_nonNA <- 20
num_top <- 10
outlier_sd <- 1.5
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")
enzyme_type <- "kinase"
reg_sig <- c(0.05, 0.05); names(reg_sig) <- c("kinase", "phosphatase")

# input druggable gene list -----------------------------------------------
drug_genes <- read.delim(file = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/reference_files/gene_drug_list/Premed_raw_databases/drugBank/drug_list.txt", header = F, col.names = "gene")
drug_genes <- as.vector(drug_genes$gene)

# input regression --------------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                   "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                    data.table = F)
regression <- regression %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP"))] | pair_pro %in% psp_tab$pair_pro)

regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression <- annotate_ks_source(regression = regression)


# calculate # ks pairs -------------------------------------------------------
tab2p <- regression %>%
  filter(GENE %in% drug_genes) %>%
  filter(regulated == T) %>%
  select(GENE, pair, Cancer) %>%
  unique() %>%
  select(GENE, Cancer) %>%
  table() %>%
  data.frame()

tab2p_text <- tab2p %>%
  filter(GENE %in% c("AKT1", "MET", "MTOR", "ERBB2", "EGFR", "MAP2K1"))

p <- ggplot()
p <- p + geom_boxplot(data = tab2p, mapping = aes(x = Cancer, y = Freq, fill = Cancer), alpha = 0.5)
p <- p + geom_point(data = tab2p, mapping = aes(x = Cancer, y = Freq), alpha = 0.8)
p <- p + geom_text_repel(data = tab2p_text, mapping = aes(x = Cancer, y = Freq, label = GENE))
p <- p + geom_line(data = tab2p_text, mapping = aes(x = Cancer, y = Freq, group = GENE), alpha = 0.6, color = "grey50", linetype = 2)
p <- p + theme_bw()
p 
ggsave(filename = paste0(makeOutDir(resultD = resultD), "number_of_regulated_druggable_kinase_substrate_pairs_per_cancer.pdf"), width = 8, height = 6)



