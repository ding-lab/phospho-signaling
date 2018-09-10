# Yige Wu @ WashU March 2018
# intergrate MSI score, gene/protein abundance and quantile for MMR genes

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# input -------------------------------------------------------------------
## input MMR genes
mmr_gene <- read_delim("~/Box Sync/MSI_CPTAC/Data/mmr_gene.txt","\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
mmr_gene <- mmr_gene$X1

## input MSI scores
msi_score <- read_delim("~/Box Sync/MSI_CPTAC/Data/MSIsensor_Score_qgao/CPTAC.MSI.score.tsv",
                        "\t", escape_double = FALSE, trim_ws = TRUE)


## input clinical info
clinical <- read.delim(paste0(inputD, "Specimen_Data_20161005_Yige_20180307.txt"), sep = "\t")
clinical <- data.frame(clinical)

## initiate list to record Mann-Whitney test results
msi_exps_list <- vector("list")
msi_expq_list <- vector("list")

## cptac collection name
cptac_collection <- c("cptac2p", "cptac3")
names(cptac_collection) <- c("CO", "UCEC")


# save expression score ---------------------------------------------------

for (cancer in c("CO", "UCEC")) {
  msi_exps_list[[cancer]] <- vector("list")

  ## input expression score
  exps <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/get_expression_quantile/", cptac_collection[cancer], "_protein_score.RDS"))
  exps <- exps[[cancer]]
  exps <- exps[exps$Gene %in% mmr_gene,]
  
  ## transform sample IDs to patient IDs
  exps$Participant.ID <- sampID2partID(sampleID_vector = exps$Specimen.ID, sample_map = clinical)
  msi_exps <- merge(msi_score, exps, by.x = c("Sample"), by.y = c("Participant.ID"))
  msi_exps$MSI_status <- ifelse(msi_exps$Score >= 3.5, paste0("MSI-H_tumors"), paste0("other_tumors"))
  msi_exps_list[[cancer]][["protein"]] <- msi_exps
}

for (cancer in c("CO", "UCEC")) {
  exps <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/get_mRNA_quantile/", cptac_collection[cancer], "_mRNA_score.RDS"))
  exps <- exps[[cancer]]
  exps <- exps[exps$Gene %in% mmr_gene,]
  ## merge exps with MSI score
  msi_exps <- merge(msi_score, exps, by.x = c("Sample"), by.y = c("Participant.ID"))
  msi_exps$MSI_status <- ifelse(msi_exps$Score >= 3.5, paste0("MSI-H_tumors"), paste0("other_tumors"))
  msi_exps$log2RPKM <- log2(msi_exps$score+1)
  msi_exps_list[[cancer]][["mRNA"]] <- msi_exps
}

for (cancer in c("CO", "UCEC")) {
  ## input expression score
  exps <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/get_phosphoprotein_quantile/", cptac_collection[cancer], "_phosphoprotein_score.RDS"))
  exps <- exps[[cancer]]
  exps <- exps[exps$Gene %in% mmr_gene,]
  
  ## transform sample IDs to patient IDs
  exps$Participant.ID <- sampID2partID(sampleID_vector = exps$Specimen.ID, sample_map = clinical)
  msi_exps <- merge(msi_score, exps, by.x = c("Sample"), by.y = c("Participant.ID"))
  msi_exps$MSI_status <- ifelse(msi_exps$Score >= 3.5, paste0("MSI-H_tumors"), paste0("other_tumors"))
  msi_exps_list[[cancer]][["phosphoprotein"]] <- msi_exps
}

resultDnow <- makeOutDir()
saveRDS(object = msi_exps_list, file = paste0(resultDnow, "mmr_gene_msi_exps_list.RDS"))


# save expression quantile ---------------------------------------------------
for (cancer in c("CO", "UCEC")) {
  msi_expq_list[[cancer]] <- vector("list")
  
  ## input expression score
  expq <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/get_expression_quantile/", cptac_collection[cancer], "_protein_quantile.RDS"))
  expq <- expq[[cancer]]
  expq <- expq[expq$Gene %in% mmr_gene,]
  
  ## transform sample IDs to patient IDs
  expq$Participant.ID <- sampID2partID(sampleID_vector = expq$Specimen.ID, sample_map = clinical)
  msi_expq <- merge(msi_score, expq, by.x = c("Sample"), by.y = c("Participant.ID"))
  msi_expq$MSI_status <- ifelse(msi_expq$Score >= 3.5, paste0("MSI-H_tumors"), paste0("other_tumors"))
  msi_expq_list[[cancer]][["protein"]] <- msi_expq
}

for (cancer in c("CO", "UCEC")) {
  expq <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/get_mRNA_quantile/", cptac_collection[cancer], "_mRNA_quantile.RDS"))
  expq <- expq[[cancer]]
  expq <- expq[expq$Gene %in% mmr_gene,]
  ## merge expq with MSI score
  msi_expq <- merge(msi_score, expq, by.x = c("Sample"), by.y = c("Participant.ID"))
  msi_expq$MSI_status <- ifelse(msi_expq$Score >= 3.5, paste0("MSI-H_tumors"), paste0("other_tumors"))
  msi_expq_list[[cancer]][["mRNA"]] <- msi_expq
}

for (cancer in c("CO", "UCEC")) {
  expq <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/get_phosphoprotein_quantile/", cptac_collection[cancer], "_phosphoprotein_quantile.RDS"))
  expq <- expq[[cancer]]
  expq <- expq[expq$Gene %in% mmr_gene,]
  ## merge expq with MSI score
  expq$Participant.ID <- sampID2partID(sampleID_vector = expq$Specimen.ID, sample_map = clinical)
  
  msi_expq <- merge(msi_score, expq, by.x = c("Sample"), by.y = c("Participant.ID"))
  msi_expq$MSI_status <- ifelse(msi_expq$Score >= 3.5, paste0("MSI-H_tumors"), paste0("other_tumors"))
  msi_expq_list[[cancer]][["phosphoprotein"]] <- msi_expq
}

resultDnow <- makeOutDir()
saveRDS(object = msi_expq_list, file = paste0(resultDnow, "mmr_gene_msi_expq_list.RDS"))



