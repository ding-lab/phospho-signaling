# Yige Wu @ WashU Aug 2018
## correct the case ID swap at RNA-seq level for PAM50 calls

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# set variables -----------------------------------------------------------
swaped_caseIDs1 <- c("11BR025", "11BR019")

# inputs ------------------------------------------------------------------
old_fpkm <- fread(input = paste0(cptac2p_genomicD, "FPKM/BRCA/11062017/fpkm_prospective_breastcancer_110617.csv"), data.table = F)
colnames(old_fpkm)
ncol(old_fpkm)
nrow(old_fpkm)
# processing --------------------------------------------------------------
partIDs <- colnames(old_fpkm)
pos1 <- which(partIDs == swaped_caseIDs1[1])
pos2 <- which(partIDs == swaped_caseIDs1[2])
partIDs[pos1] <- swaped_caseIDs1[2]
partIDs[pos2] <- swaped_caseIDs1[1]
new_fpkm <- old_fpkm
colnames(new_fpkm) <- partIDs

write.table(x = new_fpkm, file = paste0(cptac2p_genomicD, "FPKM/BRCA/fpkm_prospective_breastcancer_110617_swap_corrected.csv"),
            row.names = F, quote = F, sep = ",")



# see how many proteomics sample have fpkm ----------------------------------
pro <- loadProteinNormalizedTumor("BRCA")
sampIDmap <- loadSampMap()
partIDs <- sampID2partID(sampleID_vector = colnames(pro), sample_map = sampIDmap)
partIDs[!(partIDs %in% colnames(old_fpkm))]
## three proteomics cases don't have RNA-seq data