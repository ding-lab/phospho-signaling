# Yige Wu @ WashU Aug 2018
## correct the case ID swap at RNA-seq level for PAM50 calls

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# set variables -----------------------------------------------------------
swaped_caseIDs1 <- c("11BR025", "11BR019")

# inputs ------------------------------------------------------------------
old_pam50_map <- fread(input = paste0(cptac2p_genomicD, "FPKM/BRCA/11062017/fpkm_pam50_brca_110617_pam50scores.txt"), data.table = F)

# processing --------------------------------------------------------------
partIDs <- as.vector(old_pam50_map$V1)
pos1 <- which(partIDs == swaped_caseIDs1[1])
pos2 <- which(partIDs == swaped_caseIDs1[2])
partIDs[pos1] <- swaped_caseIDs1[2]
partIDs[pos2] <- swaped_caseIDs1[1]
new_pam50_map <- old_pam50_map
new_pam50_map$V1 <- partIDs

write.table(x = new_pam50_map, file = paste0(cptac_sharedD, "BRCA/fpkm_pam50_brca_080718_pam50scores_swap_corrected.txt"),
            row.names = F, quote = F, sep = "\t")
