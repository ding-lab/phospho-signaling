# Yige Wu @ WashU 2018 Oct
## parse outputs from multiple runs of mutational impact

# source ------------------------------------------------------------------
source("Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")

library(dplyr)


# inputs ------------------------------------------------------------------
mut_pho <- fread(input = "./cptac2p/analysis_results/phospho_network/genoalt/tables/parse_mut_impact/cptac2p_mut_cnv.txt", data.table = F)


# mutations in TP53 associated with upstream -----------------------------------------------------------
tab2w <- mut_pho
tab2w <- tab2w[tab2w$GENE == "TP53" & tab2w$SUB_GENE != "TP53",]
# tab2w <- tab2w[tab2w$p < 0.05,]
tab2w <- tab2w[order(tab2w$fdr),]

# mutations in TP53 upstream -----------------------------------------------------------
tab2w <- mut_pho
tab2w <- tab2w[tab2w$SUB_GENE == "TP53" & tab2w$GENE != "TP53",]
## no significant events associated with mutations in enzymes regulating TP53

# anything with significant FDR -----------------------------------------------------------
tab2w <- mut_pho
tab2w <- tab2w[tab2w$fdr < 0.1 & tab2w$fdr > 0,]

# mutations in AKT1 -----------------------------------------------------------
tab2w <- mut_pho
tab2w <- tab2w[tab2w$GENE == "AKT1" & tab2w$SUB_GENE != "AKT1",]
tab2w <- tab2w[order(tab2w$fdr),]
