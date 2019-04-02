# Yige Wu @ WashU 2018 July
## show the number of phosphosites detected in prospective collection overlapping with phosphoELM 

# source ------------------------------------------------------------------
source("Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R")

library(readr)
library(eulerr)
library(ggrepel)

# set variables -----------------------------------------------------------


# input CPTAC prospective phosphosite genomic coordinates(hg38) -----------------------------------
cptac_phosphosites <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid_hg38.txt"), data.table = F)
## create the column that contains the formatted genomic coordinates
cptac_phosphosites <- cptac_phosphosites %>%
  mutate(g_coord_hg38 = paste0(seq_region_name, ":g.", start, "_", end))
cptac_phosphosites %>%
  head()


# input the actual CPTAC phosphosites we are using ------------------------
cptac_phosphosites_used <- fread(input = "./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers.txt", data.table = F)
cptac_phosphosites <- cptac_phosphosites %>%
  filter(id %in% cptac_phosphosites_used$id)

# input PSP phosphosites --------------------------------------------------
psp_phosphosites <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_psp_phosphosite2genomic_coord/psp_sites_single_site_wpid_hg38.txt"), data.table = F)
psp_phosphosites <- psp_phosphosites %>%
  mutate(g_coord_hg38 = paste0(seq_region_name, ":g.", start, "_", end))


# create data frame for venn diagram --------------------------------------
df_phosphosites <- unique(c(cptac_phosphosites$g_coord_hg38, psp_phosphosites$g_coord_hg38))
df_phosphosites %>%
  head()
df_phosphosites <- data.frame(g_coord_hg38 = df_phosphosites)
df_phosphosites$CPTAC <- (df_phosphosites$g_coord_hg38 %in% cptac_phosphosites$g_coord_hg38)
df_phosphosites$PSP <- (df_phosphosites$g_coord_hg38 %in% psp_phosphosites$g_coord_hg38)


# plot venn diagram -------------------------------------------------------
tab2p <- df_phosphosites

fn = paste(makeOutDir(resultD = resultD),'venn_CPTAC_PSP.pdf',sep ="")
grid.newpage()
pdf(fn, height = 5, width = 5, useDingbats = FALSE)
fit <- euler(combinations = tab2p[, c("CPTAC", "PSP")], input = "disjoint", shape = 'ellipse')
p <-plot(fit, quantities = list(fontsize = 20), labels = list(fontsize = 20))
grid.draw(p)
dev.off()
