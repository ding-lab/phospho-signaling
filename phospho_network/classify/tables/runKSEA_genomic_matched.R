# Yige Wu @ WashU 2018 Mar
# running KSEA


# source ------------------------------------------------------------------
library(KSEAapp)
source('/Users/yigewu/Box Sync/cptac2p_analysis/cptac2p_analysis_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
makeOutDir()

# loop around cancer types using kinase-substrate data from OmniPath + NetworKIN------------------------------------------------
## input OmniPath
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin/omnipath_networkin_enzyme_substrate_site_level_union.csv"))
for (m_cutoff in 1:3) {
  for (NetworKIN.cutoff in 2:5) {
    outputDtmp = paste0("/Users/yigewu/Box\ Sync/cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA_genomic_matched/OmniPath_NetworKIN_", 
                        "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m_cutoff, "/")
    dir.create(outputDtmp)
    for (cancer in c("BRCA", "OV", "CO", "UCEC")) {
      df <- fread(input = paste0(ppnD, "classify/tables/generate_table4KSEA/", cancer, "_phosphposite_FC_for_KSEA.csv"),
                  data.table = F)
      outputDtmp = paste0(ppnD, "classify/tables/runKSEA_genomic_matched/OmniPath_NetworKIN_", 
                          "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m_cutoff, "/", cancer, "/")
      dir.create(outputDtmp)
      setwd(outputDtmp)
      KSEA.Complete(KSData =  ptms_site_pairs_sup , 
                    PX = df, NetworKIN = T, NetworKIN.cutoff = NetworKIN.cutoff, m.cutoff = m_cutoff, p.cutoff = 0.05)
    }
  }
}