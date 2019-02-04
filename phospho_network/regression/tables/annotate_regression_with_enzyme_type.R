# Yige Wu @ WashU 2019 Jan
## annotate the regresion results with supporting evidence

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source('cptac2p_analysis/phospho_network/phospho_network_shared.R')

# input regression table --------------------------------------------------
if (!exists("table2annotate")) {
  stop("input regression table named table2annotate first!")
}
table2annotate$enzyme_type <- NULL

# gather merged source table ----------------------------------------------
table2annotate$enzyme_type <- ""
table2annotate$enzyme_type[table2annotate$GENE %in% kinases] <- "kinase"
table2annotate$enzyme_type[table2annotate$GENE %in% phosphatases] <- "phosphatase"
