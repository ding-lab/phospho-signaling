# Yige Wu @ March 2019 WashU
# show the number of associated substrate phosphosites per kinase with functional annotation, especially those that haven't been reported before


# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
setwd(baseD)
source('./cptac2p_analysis/preprocess_files/preprocess_files_shared.R')
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")


# inputs regulatory sites------------------------------------------------------------------
regulatory_sites <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/Phospho_databases/PhosphositePlus/Mar_04_2019/Regulatory_sites", data.table = F)
regulatory_sites <- regulatory_sites %>%
  mutate(SUB_MOD_RSD = str_split_fixed(string = MOD_RSD, pattern = "-", n = 2)[,1])
regulatory_sites2merge <- regulatory_sites %>%
  select(GENE, SUB_MOD_RSD, ON_FUNCTION, ON_PROCESS, ON_PROT_INTERACT, ON_OTHER_INTERACT) %>%
  unique

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
table(regression$Cancer)


# annotate the regulatory sites -------------------------------------------
regression <- merge(regression, regulatory_sites2merge, by.x = c("SUBSTRATE", "SUB_MOD_RSD"), by.y = c("GENE", "SUB_MOD_RSD"), all.x = T)


# get all the terms in ON_PROCESS -----------------------------------------
processes_terms <- regression %>%
  filter(regulated == T) %>%
  select(ON_PROCESS) %>%
  unique() %>%
  unlist()

processes_terms <-  str_split(string = processes_terms, pattern = "; ") %>%
  unlist() %>%
  unique %>%
  sort()

processes_terms_estimate <- regression %>%
  filter(regulated == T) %>%
  select(ON_PROCESS) %>%
  unlist()

processes_terms_estimate <-  str_split(string = processes_terms_estimate, pattern = "; ") %>%
  unlist() %>%
  table()

processes_terms_estimate <- data.frame(ON_PROCESS_term = names(processes_terms_estimate), Freq = as.vector(processes_terms_estimate))
processes_terms_estimate <- processes_terms_estimate %>%
  filter(ON_PROCESS_term != "") %>%
  mutate(ON_PROCESS_term_parent = str_split_fixed(string = ON_PROCESS_term, pattern = ",", n = 2)[,1])
  arrange(Freq)
processes_terms_estimate

terms2process <- c("cytoskeletal reorganization", "chromatin organization",
                   "endocytosis", "DNA repair", "translation", "cell adhesion", "exocytosis", "autophagy", "cell differentiation", "carcinogenesis", "cell growth", "transcription", "apoptosis", "cell motility")
# terms2process <- c("apoptosis")
terms2process <- unique(processes_terms_estimate$ON_PROCESS_term_parent)
for (ON_PROCESS_term in terms2process) {
  # get table and make figure about ON_PROCESS_tmp -----------------------------------------------
  regression[, ON_PROCESS_term] <- sapply(1:nrow(regression), FUN = function(i, vector_function, ON_PROCESS_term) {
    if (!is.na(vector_function[i]) & grepl(x = vector_function[i], pattern = ON_PROCESS_term)) {
      tmp <- "altered"
      if (!grepl(x = vector_function[i], pattern = paste0(ON_PROCESS_term, ", inhibited")) & grepl(x = vector_function[i], pattern = paste0(ON_PROCESS_term, ", induced"))) {
        tmp <- "induced"
      }
      if (grepl(x = vector_function[i], pattern = paste0(ON_PROCESS_term, ", inhibited")) & !grepl(x = vector_function[i], pattern = paste0(ON_PROCESS_term, ", induced"))) {
        tmp <- "inhibited"
      }
    } else {
      tmp <- NA
    }
    return(tmp)
  }, vector_function = regression$ON_PROCESS, ON_PROCESS_term = ON_PROCESS_term)
}


tab2p <- regression %>%
  filter(regulated == T) %>%
  filter(SELF == "trans") %>%
  select(GENE, pair_cancer, terms2process)

tab2p <- melt(tab2p, id.vars = c("GENE", "pair_cancer"))

tab2p <- tab2p %>%
  filter(!is.na(value)) %>%
  select(GENE, pair_cancer, variable) %>%
  unique %>%
  select(GENE, variable) %>%
  table %>%
  data.frame()

# tab2p <- tab2p %>%
#   filter(Freq > 0)

mat_value <- dcast(data = tab2p, formula = variable ~ GENE, value.var = "Freq")
rownames(mat_value) <- mat_value$variable
mat_value[,1] <- NULL
mat_value <- mat_value[,colSums(mat_value, na.rm = T) > 0]
my_heatmap <- pheatmap(mat_value, 
                       color = c("white", col_paletteR(20)),
                       # breaks = breaks,
                       na_col = "white", 
                       cluster_cols=T, show_colnames = T, 
                       cluster_rows=F, show_rownames = T)
fn <- paste0(makeOutDir(resultD = resultD), "num_pair_cancer_w.functionally_annotated_sub_phosphosite.pdf")
save_pheatmap_pdf(x = my_heatmap, 
                  filename = fn, 
                  width = 18, height = 4)


tab2p <- regression %>%
  filter(regulated == T) %>%
  select(GENE, pair_cancer, terms2process)

tab2p <- melt(tab2p, id.vars = c("GENE", "pair_cancer"))

tab2p %>%
  filter(!is.na(value)) %>%
  filter(GENE == "MAPKAPK2") %>%
  filter(Cancer == "CCRCC")

kinase_dark_tab <- fread("../Downloads/All_Kinase_cleaned.txt", data.table = F)

tab2p_dark <- tab2p %>%
  filter(!is.na(value)) %>%
  filter(GENE %in% kinase_dark_tab$gene[kinase_dark_tab$is_dark_kinase])
