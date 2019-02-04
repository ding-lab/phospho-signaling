# Yige Wu @ WashU 2019 Jan
## annotate the regresion results with supporting evidence

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source('./cptac2p_analysis/phospho_network/phospho_network_shared.R')

# input regression table --------------------------------------------------
if (!exists("table2annotate")) {
  stop("input regression table first!")
}
table2annotate$Source <- NULL

# input kinase-substrate sources ------------------------------------------
omnipath_tab <- load_omnipath()
psp_tab <- load_psp()

# gather merged source table ----------------------------------------------
source_tab <- data.frame(pair = psp_tab$pair[!(psp_tab$pair %in% omnipath_tab$pair)], Source = "PhosphoSite")
source_tab <- rbind(source_tab, omnipath_tab[, c("pair", "Source")])

table2annotate <- merge(table2annotate, source_tab, all.x = T)
