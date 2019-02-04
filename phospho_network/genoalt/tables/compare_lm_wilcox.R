# Yige at WashU @ Nov 2018
## for deciding which stats should be used for mutational impact section

# source ------------------------------------------------------------------
source("Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
library(dplyr)
library(eulerr)

mut_pho_wilcox <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_cna_impact_cptac2p/mut_cnv_sig_cans.txt"), data.table = F)
mut_pho_wilcox$id <- paste0(mut_pho_wilcox$GENE, ":", mut_pho_wilcox$SUB_GENE, ":", mut_pho_wilcox$SUB_MOD_RSD, ":", mut_pho_wilcox$cancer)
mut_pho_lm <- fread(input = "./cptac2p/analysis_results/phospho_network/regression/tables/regression_enzyme_mutation/kinase_substrate_regression_cptac2p_3can_tumor.txt", data.table = F)
mut_pho_lm$id <- paste0(mut_pho_lm$KINASE, ":", mut_pho_lm$SUBSTRATE, ":", mut_pho_lm$SUB_MOD_RSD, ":", mut_pho_lm$Cancer)

nrow(mut_pho_wilcox[mut_pho_wilcox$p > 0 & mut_pho_wilcox$p <  0.05,])
nrow(mu_pho_lm[mu_pho_lm$P_mut_kin > 0 & mu_pho_lm$P_mut_kin <  0.05,])

## venn diagram for phosphosites
dat <- data.frame(id = unique(c(as.vector(mut_pho_wilcox$id), as.vector(mu_pho_lm$id))))
dat$wilcox <- (dat$id %in% mut_pho_wilcox$id[mut_pho_wilcox$p > 0 & mut_pho_wilcox$p <  0.05])
dat$lm <- (dat$id %in% mut_pho_lm$id[mu_pho_lm$P_mut_kin > 0 & mu_pho_lm$P_mut_kin <  0.05])

fn = paste(makeOutDir(resultD = resultD), 'test', '.pdf',sep ="")
grid.newpage()
pdf(fn, height = 12, width = 12, useDingbats = FALSE)
fit <- euler(combinations = dat, input = "disjoint", shape = 'circle')

p <-plot(fit, quantities = list(fontsize = 50), labels = F, legend = list(fontsize = 50))
grid.draw(p)
dev.off()

mut_pho_lm[mut_pho_lm$KINASE == "PIK3CA" & mut_pho_lm$SUBSTRATE == "AKT1",]
