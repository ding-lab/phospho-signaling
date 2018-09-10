# Yige Wu @ WashU 2018 Mar
## annotate BRCA PRM raw results with kinase familly info, LOD, LOQ, etc

# souce -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# inputs ------------------------------------------------------------------
cancer <- "BRCA"
prm_raw.f <- fread(input = paste0("cptac2p/cptac_shared/analysis_results/preprocess_files/format_RPM_peptide/prm_genename_annotated.txt"),
                   data.table = F)
prm_raw.f.cv <- prm_raw.f[,c(colnames(prm_raw.f)[grepl(pattern = "CV", x = colnames(prm_raw.f))])]
cv_nonzeroorna <- rowSums(!is.na(prm_raw.f.cv) & (prm_raw.f.cv != 0 ))
cv_median <- sapply(1:nrow(prm_raw.f), FUN = function(n, m) median(as.numeric(as.vector(m[n,])), na.rm = T), m = prm_raw.f.cv)
prm_raw.f$num_nonzeroorna <- cv_nonzeroorna
prm_raw.f$cv_median <- cv_median
num_peptide <- data.frame(table(prm_raw.f[,c("Gene_Name")]))
prm_raw.f <- merge(prm_raw.f, num_peptide, by.x = c("Gene_Name"), by.y = c("Var1"), all.x = T)
prm_loq <- fread(input = "cptac2p/reid/Merged_Kinase_Met_PCS_equations_TD_LLOQ_101917.csv", data.table = F)
prm_raw.f <- merge(prm_raw.f, prm_loq[,c("peptide", "phospho", "rlm_slope", "rlm_intercept", "LOD", "CPTAC LLOQ")], by.x = c("Peptide.Modified.Sequence"), by.y = c("peptide"), all.x = T)

## input kinase family info
kinase_family_per_synonym <- fread(input = "./cptac2p/cptac_shared/analysis_results/validation/tables/cor_prm_ms_match_specimen/kinase_family_per_synonym,txt", data.table = F)

## input metabolome annotation info
metabolism_group <- read_excel(path = "./cptac2p/cptac_shared/analysis_results/validation/tables/annotate_prm_raw/proteins_in_prm_notkinase_DAVID_functional_annotation.xlsx")
metabolism_group2m <- metabolism_group[, c("ID","panel", "GroupName")]
colnames(metabolism_group2m) <- c("Gene_Name", "panel", "GroupName")

## merge with gene annotation
prm_raw.f <- merge(prm_raw.f, kinase_family_per_synonym[, c("synonym", "GroupName", "ProteinNames")], by.x = c("Gene_Name"), by.y = c("synonym"), all.x = T)
prm_raw.f$is_kinase <- (grepl(pattern = "kinase", x = as.vector(prm_raw.f$ProteinNames)) | !is.na(prm_raw.f$ProteinNames))
prm_raw.k <- prm_raw.f[prm_raw.f$is_kinase,]
prm_raw.k$panel <- "Kinome_panel"
prm_raw.k$panel[prm_raw.k$Gene_Name == "PIK3R4"] <- "Metabolome_panel"
prm_raw.m <- prm_raw.f[!(prm_raw.f$is_kinase),]
prm_raw.m$GroupName <- NULL
prm_raw.m <- merge(prm_raw.m, metabolism_group2m, by = c("Gene_Name"), all.x = T)
prm_raw.f <- rbind(prm_raw.k, prm_raw.m[, colnames(prm_raw.k)])
write.table(prm_raw.f, row.names = F, quote=F, sep = '\t', file=paste(makeOutDir(), cancer, "_prm_genename_annotated_specimen_lod_loq_labeled.txt",sep=""))
