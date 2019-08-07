# Yige Wu @ WashU Mar 2019
## calculate the number of proteins invovled in regulated kinase-substrate pairs

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
setwd(baseD)
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")

# set variables -----------------------------------------------------------
data2process <- matrix(data = c("BRCA", "CDAP", "tumor", "scaled", "cptac2p",
                                "CO", "CDAP", "tumor", "scaled", "cptac2p",
                                "UCEC", "PGDAC", "tumor", "median_polishing", "cptac3",
                                "OV", "CDAP", "tumor", "scaled", "cptac2p",
                                "CCRCC", "PGDAC", "tumor", "MD_MAD", "cptac3"), ncol = 5, byrow = T)
reg_sig <- c(0.05, 0.05); names(reg_sig) <- c("kinase", "phosphatase")
omnipath_tab <- annotate_ks_source(omnipath_tab)

# input KEGG --------------------------------------------------------------
load(paste0(pan3can_shared_dataD, "Gene_family/2015-08-01_Gene_Set.RData"))
all_keggs <- names(KEGG)
all_keggs

# input regression --------------------------------------------------------
for (iteration in 1) {
  regression <- NULL
  for (i in 1:nrow(data2process)) {
    cancer <- data2process[i,1]
    pipeline_type <- data2process[i,2]
    sample_type <- data2process[i,3]
    norm_type <- data2process[i,4]
    cptac_phase <- data2process[i,5]
    
    tn = paste0(ppnD, "regression/tables/regression_simulate_substrate_only/", 
                "/iteration1/", cancer, "/", "regression_", cptac_phase , "_", cancer, "_", sample_type, "_", pipeline_type, "_", norm_type, "_iteration1.txt")
    tab_tmp <- fread(input = tn, data.table = F)
    tab_tmp <- tab_tmp %>%
      filter(!is.na(KINASE)) %>%
      mutate(SUB_GENE = SUBSTRATE)
    regression <- rbind(regression, tab_tmp)
    
    
  }
  regression %>% nrow()
  regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
  
}



# test BRCA ---------------------------------------------------------------
cancer_tmp <- "BRCA"
path_tmp <- ""
key <- "PI3K"
kegg_name_tmp <- all_keggs[grepl(key, all_keggs, ignore.case = T)]
kegg_name_tmp
regression_tmp <- regression %>%
  filter(regulated == T) %>%
  filter(Cancer == cancer_tmp) %>%
  filter(SELF == "trans") %>%
  filter(GENE %in% KEGG[[kegg_name_tmp]]) %>%
  filter(SUB_GENE %in% KEGG[[kegg_name_tmp]])

nrow(regression_tmp)
