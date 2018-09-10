# Yige Wu @ WashU 2018 Mar
# running KSEA


# source ------------------------------------------------------------------
library(KSEAapp)
source('/Users/yigewu/Box Sync/cptac2p_analysis/cptac2p_analysis_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# inputs ------------------------------------------------------------------
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin/omnipath_networkin_enzyme_substrate_site_level_union.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[order(ptms_site_pairs_sup$GENE, ptms_site_pairs_sup$SUB_GENE, ptms_site_pairs_sup$SUB_MOD_RSD, ptms_site_pairs_sup$networkin_score, str_length(ptms_site_pairs_sup$Source), decreasing = T),]
ptms_site_pairs_sup$pair <- paste0(ptms_site_pairs_sup$GENE, ":", ptms_site_pairs_sup$SUB_GENE, ":", ptms_site_pairs_sup$SUB_MOD_RSD)
ptms_site_pairs_sup <- ptms_site_pairs_sup[!duplicated(ptms_site_pairs_sup$pair),]
ptms_site_pairs_sup$pair <- NULL
ptms_site_pairs_sup$Source[ptms_site_pairs_sup$Source != "NetKIN"] <- "Signor"
# ptms_site_pairs_sup[ptms_site_pairs_sup$GENE == "PTPN11" & ptms_site_pairs_sup$SUB_GENE == "MPZL1" & ptms_site_pairs_sup$SUB_MOD_RSD == "Y263", "Source"] <- "Signor"
outputD1 <- makeOutDir(resultD = resultD)


# loop around cancer types using kinase-substrate data from OmniPath + NetworKIN------------------------------------------------
## input OmniPath

for (m_cutoff in 2:2) {
  for (NetworKIN.cutoff in 5:5) {
    outputD2 = paste0(outputD1,
                        "OmniPath_NetworKIN_", "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m_cutoff, "/")
    dir.create(outputD2)
    for (cancer in c("BRCA", "OV", "CO", "UCEC")) {
      df <- fread(input = paste0(ppnD, "classify/tables/generate_table4KSEA/", cancer, "_phosphposite_FC_for_KSEA.csv"),
                  data.table = F)
      df <- df[order(df$Gene, df$Residue.Both),]
      outputD3 = paste0(outputD2, cancer, "/")
      dir.create(outputD3)
      setwd(outputD3)
      KSEA.Complete(KSData =  ptms_site_pairs_sup , 
                    PX = df, NetworKIN = T, NetworKIN.cutoff = NetworKIN.cutoff, m.cutoff = m_cutoff, p.cutoff = 0.05)
      source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
    }
  }
}