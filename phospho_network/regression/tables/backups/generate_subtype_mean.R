# Yige Wu @ WashU 2018 Jan
# Usage: generate subtype mean protein/phosphoprotien abundance

# fixed input -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
sample_map <- read.delim(paste0(inputD, "Specimen_Data_20161005_Yige_20180127.txt"))

# BRCA pam50 --------------------------------------------------------------
cancer <- "BRCA"
subtype_map <- read_excel("~/Box Sync/cptac2p/cptac_shared/5_CPTAC2_Breast_Prospective_Collection_BI/proteome-data-v1.01011/data/20171106_CPTAC2_ProspectiveBreastCancer_Sample_Annotations_Table_v79.xlsx")
pho_t <- fread(input = paste0(inputD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt"), data.table = F)
pho_n <- fread(input = paste0(inputD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_Control.txt"), data.table = F)
sampleIDs <- colnames(pho_t)[!(colnames(pho_t) %in% c("Gene", "Phosphosite"))]
pam50 <- sampID2pam50(sampleIDs, subtype_map)
cohorts <- c("LumA","LumB","Her2","Basal")

## generate phosphosite mean
pho_subtype_mean <- data.frame(matrix(rep(vector(mode = "numeric", length = nrow(pho_t)), length(cohorts)), ncol=length(cohorts), byrow=T))
colnames(pho_subtype_mean) <- cohorts
for (cohort in cohorts) {
  subtype_sample <- sampleIDs[pam50 == cohort]
  pho_subtype_mean[,cohort] <- rowMeans(pho_t[,c(subtype_sample)], na.rm = TRUE)
}
pho_subtype_mean[,'Normal'] <- rowMeans(pho_n[,!(colnames(pho_n) %in% c("Gene", "Phosphosite"))], na.rm = T)
pho_subtype_mean <- cbind(pho_t[,c("Gene", "Phosphosite")], pho_subtype_mean)
write.table(pho_subtype_mean, row.names = F, quote=F, sep = '\t', file=paste(inputD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_pam50_mean.txt",sep=""))

## generate protein mean
pro_t <- fread(input = paste0(inputD, cancer,"/",prefix[cancer], "_PRO_formatted_normalized_noControl.txt"), data.table = F)
pro_n <- fread(input = paste0(inputD, cancer,"/",prefix[cancer], "_PRO_formatted_normalized_Control.txt"), data.table = F)
pro_subtype_mean <- data.frame(matrix(rep(vector(mode = "numeric", length = nrow(pro_t)),4), ncol=4, byrow=T))
colnames(pro_subtype_mean) <- cohorts
for (cohort in cohorts) {
  subtype_sample <- sampleIDs[pam50 == cohort]
  pro_subtype_mean[,cohort] <- rowMeans(pro_t[,c(subtype_sample)], na.rm = TRUE)
}
pro_subtype_mean[,'Normal'] <- rowMeans(pro_n[,!(colnames(pro_n) %in% c("Gene", "Phosphosite"))], na.rm = T)
pro_subtype_mean <- cbind(pro_t$Gene, pro_subtype_mean); colnames(pro_subtype_mean)[1] <- "Gene"
write.table(pro_subtype_mean, row.names = F, quote=F, sep = '\t', file=paste(inputD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_pam50_mean.txt",sep=""))

## generate phosphoprotein mean
phog_t <- fread(input = paste0(baseD,"cptac2p/cptac_shared/",cancer,"/",prefix[cancer], "_PHO_normalized_collapsed_noControl.txt"), data.table = F)
phog_n <- fread(input = paste0(baseD,"cptac2p/cptac_shared/",cancer,"/",prefix[cancer], "_PHO_normalized_collapsed_Control.txt"), data.table = F)
phog_subtype_mean <- data.frame(matrix(rep(vector(mode = "numeric", length = nrow(phog_t)), length(cohorts)), ncol=length(cohorts), byrow=T))
colnames(phog_subtype_mean) <- cohorts
for (cohort in cohorts) {
  subtype_sample <- sampleIDs[pam50 == cohort]
  phog_subtype_mean[,cohort] <- rowMeans(phog_t[,c(subtype_sample)], na.rm = TRUE)
}
phog_subtype_mean[,'Normal'] <- rowMeans(phog_n[,!(colnames(phog_n) %in% c("Gene", "Phosphosite"))], na.rm = T)
phog_subtype_mean <- cbind(phog_t$Gene, phog_subtype_mean); colnames(phog_subtype_mean)[1] <- "Gene"
write.table(phog_subtype_mean, row.names = F, quote=F, sep = '\t', file=paste(inputD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_collapsed_pam50_mean.txt",sep=""))

# CO MSI/MSS --------------------------------------------------------------
cancer <- "CO"
pho_t <- fread(input = paste0(inputD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt"), data.table = F)
pro_t <- fread(input = paste0(inputD, cancer,"/",prefix[cancer], "_PRO_formatted_normalized_noControl.txt"), data.table = F)
pho_n <- fread(input = paste0(inputD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_Control.txt"), data.table = F)
pro_n <- fread(input = paste0(inputD, cancer,"/",prefix[cancer], "_PRO_formatted_normalized_Control.txt"), data.table = F)
subtype_map <- read_excel("~/Box Sync/cptac2p/cptac_shared/CPTAC_Biospecimens_Clinical_Data/17_September_2016/CPTAC Colon Clinical Data_20160915.xls")
sampleIDs <- colnames(pho_t)[!(colnames(pho_t) %in% c("Gene", "Phosphosite"))]
msi <- sampID2msi(sampleIDs, sample_map, subtype_map)
cohorts <- c("MSI", "MSS", "other")

## generate phosphosite mean
pho_subtype_mean <- data.frame(matrix(rep(vector(mode = "numeric", length = nrow(pho_t)), length(cohorts)), ncol=length(cohorts), byrow=T))
colnames(pho_subtype_mean) <- cohorts
for (cohort in cohorts) {
  subtype_sample <- sampleIDs[msi == cohort]
  pho_subtype_mean[,cohort] <- rowMeans(pho_t[,c(subtype_sample)], na.rm = TRUE)
}
pho_subtype_mean[,'Normal'] <- rowMeans(pho_n[,!(colnames(pho_n) %in% c("Gene", "Phosphosite"))], na.rm = T)
pho_subtype_mean <- cbind(pho_t[,c("Gene", "Phosphosite")], pho_subtype_mean)
write.table(pho_subtype_mean, row.names = F, quote=F, sep = '\t', file=paste(inputD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_MSI_type_mean.txt",sep=""))

## generate protein mean
pro_subtype_mean <- data.frame(matrix(rep(vector(mode = "numeric", length = nrow(pro_t)), length(cohorts)), ncol=length(cohorts), byrow=T))
colnames(pro_subtype_mean) <- cohorts
for (cohort in cohorts) {
  subtype_sample <- sampleIDs[msi == cohort]
  pro_subtype_mean[,cohort] <- rowMeans(pro_t[,c(subtype_sample)], na.rm = TRUE)
}
pro_subtype_mean[,'Normal'] <- rowMeans(pro_n[,!(colnames(pro_n) %in% c("Gene", "Phosphosite"))], na.rm = T)
pro_subtype_mean <- cbind(pro_t$Gene, pro_subtype_mean); colnames(pro_subtype_mean)[1] <- "Gene"
write.table(pro_subtype_mean, row.names = F, quote=F, sep = '\t', file=paste(inputD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_MSI_type_mean.txt",sep=""))

## generate phosphoprotein mean
phog_t <- fread(input = paste0(baseD,"cptac2p/cptac_shared/",cancer,"/",prefix[cancer], "_PHO_normalized_collapsed_noControl.txt"), data.table = F)
phog_n <- fread(input = paste0(baseD,"cptac2p/cptac_shared/",cancer,"/",prefix[cancer], "_PHO_normalized_collapsed_Control.txt"), data.table = F)

phog_subtype_mean <- data.frame(matrix(rep(vector(mode = "numeric", length = nrow(phog_t)), length(cohorts)), ncol=length(cohorts), byrow=T))
colnames(phog_subtype_mean) <- cohorts
for (cohort in cohorts) {
  subtype_sample <- sampleIDs[msi == cohort]
  phog_subtype_mean[,cohort] <- rowMeans(phog_t[,c(subtype_sample)], na.rm = TRUE)
}
phog_subtype_mean[,'Normal'] <- rowMeans(phog_n[,!(colnames(phog_n) %in% c("Gene", "Phosphosite"))], na.rm = T)
phog_subtype_mean <- cbind(phog_t$Gene, phog_subtype_mean); colnames(phog_subtype_mean)[1] <- "Gene"
write.table(phog_subtype_mean, row.names = F, quote=F, sep = '\t', file=paste(inputD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_collapsed_MSI_type_mean.txt",sep=""))


# OV vs Basal --------------------------------------------------------------
cancer <- "OV"
cohorts <- c("Tumor")
## generate phosphosite mean
pho_t <- fread(input = paste0(inputD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt"), data.table = F)
sampleIDs <- colnames(pho_t)[!(colnames(pho_t) %in% c("Gene", "Phosphosite"))]
pho_n <- fread(input = paste0(inputD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_Control.txt"), data.table = F)
pho_brca <- fread(input = paste0(inputD, "BRCA","/",prefix["BRCA"], "_PHO_formatted_normalized_pam50_mean.txt"), data.table = F)
pho_subtype_mean <- data.frame(matrix(rep(vector(mode = "numeric", length = nrow(pho_t)), length(cohorts)), ncol=length(cohorts), byrow=T))
colnames(pho_subtype_mean) <- cohorts
for (cohort in cohorts) {
  subtype_sample <- sampleIDs
  pho_subtype_mean[,cohort] <- rowMeans(pho_t[,c(subtype_sample)], na.rm = TRUE)
}
pho_subtype_mean <- cbind(pho_t[,c("Gene", "Phosphosite")], pho_subtype_mean)
pho_subtype_mean <- merge(pho_subtype_mean, 
                          pho_brca[,c("Gene", "Phosphosite", "Basal", "Her2", "LumA", "LumB")],
                          by = c("Gene", "Phosphosite"), all.x = T, sort = F)
pho_subtype_mean[,'Normal'] <- rowMeans(pho_n[,!(colnames(pho_n) %in% c("Gene", "Phosphosite"))], na.rm = T)
write.table(pho_subtype_mean, row.names = F, quote=F, sep = '\t', file=paste(inputD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_tumor_normal_WBasal_mean.txt",sep=""))

## generate protein mean
pro_t <- fread(input = paste0(inputD, cancer,"/",prefix[cancer], "_PRO_formatted_normalized_noControl.txt"), data.table = F)
pro_n <- fread(input = paste0(inputD, cancer,"/",prefix[cancer], "_PRO_formatted_normalized_Control.txt"), data.table = F)
pro_brca <- fread(input = paste(inputD, "BRCA", "/", prefix["BRCA"], "_PRO_formatted_normalized_pam50_mean.txt",sep=""), data.table = F)
pro_subtype_mean <- data.frame(matrix(rep(vector(mode = "numeric", length = nrow(pro_t)), length(cohorts)), ncol=length(cohorts), byrow=T))
colnames(pro_subtype_mean) <- cohorts
for (cohort in cohorts) {
  subtype_sample <- sampleIDs
  pro_subtype_mean[,cohort] <- rowMeans(pro_t[,c(subtype_sample)], na.rm = TRUE)
}
pro_subtype_mean <- cbind(pro_t$Gene, pro_subtype_mean); colnames(pro_subtype_mean)[1] <- "Gene"
pro_subtype_mean <- merge(pro_subtype_mean, pro_brca[,c("Gene", "Basal", "Her2", "LumA", "LumB")], by = c("Gene"), all.x = T, sort = F)
pro_subtype_mean[,'Normal'] <- rowMeans(pro_n[,!(colnames(pro_n) %in% c("Gene", "Phosphosite"))], na.rm = T)
write.table(pro_subtype_mean, row.names = F, quote=F, sep = '\t', file=paste(inputD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_tumor_normal_WBasal_mean.txt",sep=""))

## generate phosphoprotein mean
phog_t <- fread(input = paste0(baseD,"cptac2p/cptac_shared/",cancer,"/",prefix[cancer], "_PHO_normalized_collapsed_noControl.txt"), data.table = F)
phog_n <- fread(input = paste0(baseD,"cptac2p/cptac_shared/",cancer,"/",prefix[cancer], "_PHO_normalized_collapsed_Control.txt"), data.table = F)
phog_brca <- fread(input = paste(inputD, "BRCA", "/", prefix["BRCA"], "_PHO_formatted_normalized_collapsed_pam50_mean.txt",sep=""), data.table = F)
phog_subtype_mean <- data.frame(matrix(rep(vector(mode = "numeric", length = nrow(phog_t)), length(cohorts)), ncol=length(cohorts), byrow=T))
colnames(phog_subtype_mean) <- cohorts
for (cohort in cohorts) {
  subtype_sample <- sampleIDs
  phog_subtype_mean[,cohort] <- rowMeans(phog_t[,c(subtype_sample)], na.rm = TRUE)
}
phog_subtype_mean <- cbind(phog_t$Gene, phog_subtype_mean); colnames(phog_subtype_mean)[1] <- "Gene"
phog_subtype_mean <- merge(phog_subtype_mean, phog_brca[,c("Gene","Basal", "Her2", "LumA", "LumB")], by = c("Gene"), all.x = T, sort = F)
phog_subtype_mean[,'Normal'] <- rowMeans(phog_n[,!(colnames(phog_n) %in% c("Gene", "Phosphosite"))], na.rm = T)
write.table(phog_subtype_mean, row.names = F, quote=F, sep = '\t', file=paste(inputD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_collapsed_tumor_normal_WBasal_mean.txt",sep=""))



