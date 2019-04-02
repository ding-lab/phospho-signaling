# Yige Wu @ March 2019 WashU
# show the number of associated substrate phosphosites per kinase with functional annotation, especially those that haven't been reported before


# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
setwd(baseD)
source('./cptac2p_analysis/preprocess_files/preprocess_files_shared.R')
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")

# input regression --------------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                   "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                    data.table = F)
reg_sig <- c(0.05, 0.05); names(reg_sig) <- c("kinase", "phosphatase")
regression %>% nrow()
regression <- regression %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP"))] | pair_pro %in% psp_tab$pair_pro)
regression %>% nrow()

regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression <- annotate_ks_source(regression = regression)



# set variables -----------------------------------------------------------
# genes2process <- c("MET")
# genes2process <- c("BRAF")
genes2process <- c("RAF1")
# genes2process <- c("PTK2")

# cancers2process <- c("CCRCC")
cancers2process <- unique(regression$Cancer)

cancer2ProteinPaintColor <- function(vector_cancer_type) {
  vector_color_string <- vector(mode = "character", length = length(vector_cancer_type))
  vector_color_string[vector_cancer_type == "CCRCC"] <- "M"
  vector_color_string[vector_cancer_type == "UCEC"] <- "P"
  vector_color_string[vector_cancer_type == "CO"] <- "S"
  vector_color_string[vector_cancer_type == "OV"] <- "F"
  vector_color_string[vector_cancer_type == "BRCA"] <- "deletion"
  return(vector_color_string)
}

# Write table -------------------------------------------------------------
for (gene_tmp in genes2process) {
  for (cancer_tmp in cancers2process) {
    regression_tmp <- regression %>%
      # filter(SELF == "cis") %>%
      filter(regulated == T) %>%
      filter(SUB_GENE %in% gene_tmp) %>%
      filter(Cancer %in% cancer_tmp) %>%
      mutate(p_coord = str_split_fixed(string = SUB_MOD_RSD, pattern = "[STY]", 3)[,2]) %>%
      mutate(is_single = (str_split_fixed(string = SUB_MOD_RSD, pattern = "[STY]", 3)[,3] == "")) %>%
      filter(is_single == T)
    regression_tmp$color <- cancer2ProteinPaintColor(regression_tmp$Cancer)
    
    table2w <- regression_tmp %>%
      select(SUB_MOD_RSD, p_coord, color) %>%
      unique()
    
    write.table(x = table2w, file = paste0(makeOutDir(resultD = resultD), cancer_tmp, "_", gene_tmp, ".txt"), sep = ";", quote = F, row.names = F, col.names = F)
  }
  regression_tmp <- regression %>%
    # filter(SELF == "cis") %>%
    filter(regulated == T) %>%
    filter(SUB_GENE %in% gene_tmp) %>%
    mutate(p_coord = str_split_fixed(string = SUB_MOD_RSD, pattern = "[STY]", 3)[,2]) %>%
    mutate(is_single = (str_split_fixed(string = SUB_MOD_RSD, pattern = "[STY]", 3)[,3] == "")) %>%
    filter(is_single == T)
  regression_tmp$color <- cancer2ProteinPaintColor(regression_tmp$Cancer)
  
  table2w <- regression_tmp %>%
    select(SUB_MOD_RSD, p_coord, color) %>%
    unique()
  
  write.table(x = table2w, file = paste0(makeOutDir(resultD = resultD), gene_tmp, ".txt"), sep = ";", quote = F, row.names = F, col.names = F)
}
