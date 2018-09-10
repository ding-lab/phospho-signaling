# Yige Wu @ WashU 2018 Mar
# running KSEA


# source ------------------------------------------------------------------
library(KSEAapp)


# inputs ------------------------------------------------------------------
PSP_NetworKIN_Kinase_Substrate_Dataset_July2016 <- read_csv("~/Box Sync/cptac2p/KSEA-master/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv")
PSP_NetworKIN_Kinase_Substrate_Dataset_July2016 <- data.frame(PSP_NetworKIN_Kinase_Substrate_Dataset_July2016)
makeOutDir()


# loop around cancer types using kinase-substrate data from OmniPath------------------------------------------------
## input OmniPath
ptms_site_pairs_sup <- read_csv("~/Box Sync/cptac2p/cptac_shared/analysis_results/phospho_network/compile_enzyme_substrate/compile_omnipath/omnipath_phosphositeplus_depod_enzyme_substrate_site_level_union.csv")


outputDtmp = paste0("/Users/yigewu/Box\ Sync/cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA/OmniPath", "/")
dir.create(outputDtmp)
for (cancer in c("BRCA", "OV", "CO")) {
  df <- fread(input = paste0("/Users/yigewu/Box\ Sync/", "cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/generate_table4KSEA/", cancer, "_phosphposite_FC_for_KSEA.csv"),
              data.table = F)
  outputDtmp = paste0("/Users/yigewu/Box\ Sync/cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA/OmniPath/", cancer, "/")
  dir.create(outputDtmp)
  setwd(outputDtmp)
  KSEA.Complete(KSData =  ptms_site_pairs_sup , 
                PX = df, NetworKIN = T, NetworKIN.cutoff = 2, m.cutoff = 1, p.cutoff = 0.05)
}

# loop around cancer types using kinase-substrate data from OmniPath + NetworKIN------------------------------------------------
## input OmniPath
ptms_site_pairs_sup <- read_csv("~/Box Sync/cptac2p/cptac_shared/analysis_results/phospho_network/compile_enzyme_substrate/compile_omnipath/omnipath_networkin_enzyme_substrate_site_level_union.csv")
# op_site_pairs_sup <- read_csv("~/Box Sync/cptac2p/cptac_shared/analysis_results/phospho_network/compile_enzyme_substrate/compile_omnipath/omnipath_phosphositeplus_depod_enzyme_substrate_site_level_union.csv")
# op_site_pairs_sup <- data.frame(op_site_pairs_sup)
# ptms_site_pairs_sup <- rbind(PSP_NetworKIN_Kinase_Substrate_Dataset_July2016[PSP_NetworKIN_Kinase_Substrate_Dataset_July2016$Source == "NetworKIN",], op_site_pairs_sup)
for (m_cutoff in 1:3) {
  for (NetworKIN.cutoff in 2:5) {
    outputDtmp = paste0("/Users/yigewu/Box\ Sync/cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA/OmniPath_NetworKIN_", 
                        "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m_cutoff, "/")
    dir.create(outputDtmp)
    for (cancer in c("BRCA", "OV", "CO")) {
      df <- fread(input = paste0("/Users/yigewu/Box\ Sync/", "cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/generate_table4KSEA/", cancer, "_phosphposite_FC_for_KSEA.csv"),
                  data.table = F)
      outputDtmp = paste0("/Users/yigewu/Box\ Sync/cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA/OmniPath_NetworKIN_", 
                          "NetworKIN.cutoff", NetworKIN.cutoff, "_m.cutoff", m_cutoff, "/", cancer, "/")
      dir.create(outputDtmp)
      setwd(outputDtmp)
      KSEA.Complete(KSData =  ptms_site_pairs_sup , 
                    PX = df, NetworKIN = T, NetworKIN.cutoff = NetworKIN.cutoff, m.cutoff = m_cutoff, p.cutoff = 0.05)
    }
  }
}


# loop around cancer types using kinase-substrate data from PhosphositePlus + NetworKIN------------------------------------------------
ptms_site_pairs_sup <- PSP_NetworKIN_Kinase_Substrate_Dataset_July2016
for (NetworKIN.cutoff in 2:5) {
  outputDtmp = paste0("/Users/yigewu/Box\ Sync/cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA/PSP_NetworKIN", NetworKIN.cutoff, "cutoff", "/")
  dir.create(outputDtmp)
  for (cancer in c("BRCA", "OV", "CO")) {
    df <- fread(input = paste0("/Users/yigewu/Box\ Sync/", "cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/generate_table4KSEA/", cancer, "_phosphposite_FC_for_KSEA.csv"),
                data.table = F)
    outputDtmp = paste0("/Users/yigewu/Box\ Sync/cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA/PSP_NetworKIN", NetworKIN.cutoff, "cutoff", "/", cancer, "/")
    dir.create(outputDtmp)
    setwd(outputDtmp)
    KSEA.Complete(KSData =  ptms_site_pairs_sup ,
                  PX = df, NetworKIN = T, NetworKIN.cutoff = NetworKIN.cutoff, m.cutoff = 1, p.cutoff = 0.05)
  }
}





# loop around cancer types using kinase-substrate data from PhosphositePlus------------------------------------------------
ptms_site_pairs_sup <- PSP_NetworKIN_Kinase_Substrate_Dataset_July2016
outputDtmp = paste0("/Users/yigewu/Box\ Sync/cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA/PSP/")
dir.create(outputDtmp)
for (cancer in c("BRCA", "OV", "CO")) {
  df <- fread(input = paste0("/Users/yigewu/Box\ Sync/", "cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/generate_table4KSEA/", cancer, "_phosphposite_FC_for_KSEA.csv"),
              data.table = F)
  outputDtmp = paste0("/Users/yigewu/Box\ Sync/cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/runKSEA/PSP/", cancer, "/")
  dir.create(outputDtmp)
  setwd(outputDtmp)
  KSEA.Complete(KSData =  ptms_site_pairs_sup , 
                PX = df, NetworKIN = F, m.cutoff = 1, p.cutoff = 0.05)
}


