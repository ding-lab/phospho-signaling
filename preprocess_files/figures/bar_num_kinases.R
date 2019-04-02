# Yige Wu @ WashU Mar. 2019
## plot the number of kinases within the phosphorysite data

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


# input filtered phosphosites ---------------------------------------------
pho_head_cancers <- fread(input = "./cptac2p/analysis_results/preprocess_files/tables/filter_phosphosites_for_ks/pho_head_cancers_nonNA201IQR_below_Median.txt", data.table = F)

# input the list of kinases -----------------------------------------------
kinases_all <- fread(input = "./cptac2p/resources/gene_lists/All_Kinase_cleaned.txt", data.table = F)
kinases_all %>%
  nrow()


# input list of phosphatases ----------------------------------------------
phosphatases_all <- readxl::read_xlsx("./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/Phospho_databases/DEPOD/DEPOD_201410_human_phosphatases.xlsx")
phosphatases_all %>%
  nrow()

# input regression --------------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                   "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                    data.table = F)
regression <- regression %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP"))] | pair_pro %in% psp_tab$pair_pro)

regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression <- annotate_ks_source(regression = regression)


# construct the data frame for plotting -----------------------------------
tab2p <- unique(c(kinases_all$gene, phosphatases_all$`Gene symbol`))
tab2p <- data.frame(gene_symbol = tab2p, enzyme_type = ifelse(tab2p %in% kinases_all$gene, "kinase", "phosphatase"))
tab2p$enzyme_in_phospho <- (tab2p$gene %in% pho_head_cancers$Gene)
tab2p$enzyme_substrate_in_phosho <- (tab2p$gene %in% regression$GENE)
tab2p$bar_fill <- "enzyme"
tab2p$bar_fill[tab2p$enzyme_in_phospho] <- "enzyme_in_phospho"
tab2p$bar_fill[tab2p$enzyme_substrate_in_phosho] <- "enzyme_substrate_in_phosho"

p <- ggplot()
p <- p + geom_bar(data = tab2p, mapping = aes(x = enzyme_type, fill = bar_fill), stat = "count", position = "stack")
p <- p + scale_fill_manual(values = c("enzyme" = "grey50", "enzyme_in_phospho" = "#addd8e", "enzyme_substrate_in_phosho" = "#238443"))
p <- p + ylab('number of kinases/phosphatases')
p <- p + theme_bw() + theme(axis.title.x = element_blank())
p
ggsave(filename = paste0(makeOutDir(resultD = resultD), "bar_num_kinases_phosphatases.pdf"), width = 5, height = 4)


# get the numbers in the plot ---------------------------------------------
tab2p %>%
  filter(enzyme_type == "kinase") %>%
  filter(enzyme_in_phospho == T) %>%
  nrow()

tab2p %>%
  filter(enzyme_type == "phosphatase") %>%
  filter(enzyme_in_phospho == T) %>%
  nrow()

tab2p %>%
  filter(enzyme_type == "kinase") %>%
  filter(enzyme_substrate_in_phosho == T) %>%
  nrow()

tab2p %>%
  filter(enzyme_type == "phosphatase") %>%
  filter(enzyme_substrate_in_phosho == T) %>%
  nrow()

regression %>%
  filter(enzyme_type == "kinase") %>%
  select(SUB_GENE) %>%
  unique %>%
  nrow()

regression %>%
  filter(enzyme_type == "phosphatase") %>%
  select(SUB_GENE) %>%
  unique %>%
  nrow()

