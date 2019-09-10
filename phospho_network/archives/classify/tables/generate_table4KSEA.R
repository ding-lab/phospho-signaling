# Yige Wu @ WashU 2018 Mar
# generate the table for running KSEA


# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/get_tumor_normal_pair.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/cptac2p_analysis_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
resultDnow <- makeOutDir()


# generate phosphosite fold change btw tumor and normal ------------------------------------------------
# for (isimputed in c("", "_imputed")) {
for (isimputed in c("")) {
  for (cancer in c("UCEC")) {
  # for (cancer in c("OV", "BRCA", "CO")) {
    pho <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PHO", "_formatted_normalized_replicate_averaged", isimputed, ".txt"),
                 data.table = F)
    df_tn_match <- tumor_normal_match[[cancer]]
    tumorSampIDs <- as.vector(df_tn_match$tumorSampIDs[!is.na(df_tn_match$normalSampIDs)])
    normalSampIDs <- as.vector(df_tn_match[tumorSampIDs, "normalSampIDs"])
    pho_head <- formatPhosphosite(phosphosite_vector = as.vector(pho$Phosphosite), gene_vector = as.vector(pho$Gene))
    
    pho.t <- pho[,tumorSampIDs]; pho.n <- pho[, normalSampIDs]
    fold_change <- 2^(rowMeans(pho.t, na.rm = T) - rowMeans(pho.n, na.rm = T))
    
    df <- data.frame(Protein = rep("NULL", times = c(nrow(pho_head))),
                     Gene = pho_head$SUBSTRATE,
                     Peptide = rep("NULL", times = c(nrow(pho_head))),
                     Residue.Both = pho_head$SUB_MOD_RSD,
                     p = rep("NULL", times = c(nrow(pho_head))),
                     FC = fold_change)
    write.table(x = df, file = paste0(makeOutDir(), cancer, "_phosphposite_FC_for_KSEA", isimputed, ".csv"),
                row.names = F, quote = F, sep = ",")
  }
}

# generate global phosphorylation (phosphoprotein level) fold change btw tumor and normal ------------------------------------------------
for (isimputed in c("", "_imputed")) {
  for (cancer in c("OV", "BRCA", "CO")) {
    phog <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "collapsed_PHO", "_formatted_normalized_replicate_averaged", isimputed, ".txt"),
                 data.table = F)
    df_tn_match <- tumor_normal_match[[cancer]]
    tumorSampIDs <- as.vector(df_tn_match$tumorSampIDs[!is.na(df_tn_match$normalSampIDs)])
    normalSampIDs <- as.vector(df_tn_match[tumorSampIDs, "normalSampIDs"])
    
    pho.t <- phog[,tumorSampIDs]; pho.n <- phog[, normalSampIDs]
    fold_change <- 2^(rowMeans(pho.t, na.rm = T) - rowMeans(pho.n, na.rm = T))
    
    df <- data.frame(Protein = rep("NULL", times = c(nrow(phog))),
                     Gene = phog$Gene,
                     Peptide = rep("NULL", times = c(nrow(phog))),
                     Residue.Both = rep("NULL", times = c(nrow(phog))),
                     p = rep("NULL", times = c(nrow(phog))),
                     FC = fold_change)
    write.table(x = df, file = paste0(makeOutDir(), cancer, "_collapsed_PHO_FC_for_KSEA", isimputed, ".csv"),
                row.names = F, quote = F, sep = ",")
  }
}
