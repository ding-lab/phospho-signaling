# Yige Wu @ WashU Mar 2018 
# preprocess the proteome and phosphoproteome data
# remember to check no duplicated samples from the same participant

# source -------------------------------------------------------------------
source('~/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files_shared.R')

# add cancer type column clinical file --------------------------------------------------
clinical <- fread(input = paste0(cptac_sharedD, "Specimen_Data_20161005_Yige_20180307.txt"), data.table = F)
if (!("cancer" %in% colnames(clinical))) {
  tmp <- vector(mode = "character", length = nrow(clinical))
  partIDs <- as.vector(clinical$Participant.ID)
  tmp[grepl(pattern = "BR", x = partIDs)] <- "BRCA"
  tmp[grepl(pattern = "OV", x = partIDs)] <- "OV"
  tmp[grepl(pattern = "CO", x = partIDs)] <- "CO"
  which(tmp == "")
  clinical$cancer <- tmp
  clinical <- clinical[!(clinical$cancer == ""),]
} 

# initiate summary table --------------------------------------------------
## need to edit the prefix vector first in the preprocess_files.R
cancers <- c("UCEC", "CCRCC", "LUAD")
sum_table <- matrix(NA, nrow = length(cancers), ncol = 12, 
                    dimnames = list(row=prefix[cancers], 
                                    col=c("samples_all_pro","samples_tumor_pro", "samples_normal_pro",
                                          "genes_pro","genes_10NonNA_pro",
                                          "samples_all_pho","samples_tumor_pho", "samples_normal_pho",
                                          "genes_pho","genes_10NonNA_pho",
                                          "phosphosites","phosphosites_10NonNA")))


# process proteome and phosphoproteome data for UCEC---------------------------------------------
raw["UCEC","protein"] <- paste(base_ucec_pnnl,"/UCEC_batch1-4_SummaryReports/CPTAC3_Uterine_Corpus_Endometrial_Carcinoma_Proteome.tmt10.tsv", sep = "")
raw["UCEC", "phosphoprotein"] <- paste(base_ucec_pnnl,"/UCEC_batch1-4_SummaryReports/CPTAC3_Uterine_Corpus_Endometrial_Carcinoma_Phosphoproteome.phosphosite.tmt10.tsv", sep = "")

for (cancer in c("UCEC")) {
  sink(paste0(cptac_sharedD, cancer, "/", 'preprocessing_log.txt'), append=FALSE, split=FALSE)
  
  ## process proteome
  print(paste0('-Start processing ', cancer, ' proteome data at:'))
  print(raw[cancer,"protein"])
  
  ## input proteome data
  Pro <- fread(raw[cancer, "protein"], data.table=FALSE)
  print(paste0('proteome data has ', ncol(Pro), ' columns'))
  
  ## input ID mapping data from supporting documents
  # id_map <- read_excel("~/Box Sync/MSI_CPTAC/Data/clinical_info/UCEC_TMT_10282017_batches.xlsx")
  id_map <- read_excel(paste0(cptac_sharedD, "CPTAC_Biospecimens_Clinical_Data/UCEC_TMT_Sample_Mapping_AllBatches.xlsx"))
  id_map <- id_map[!is.na(id_map$`Aliquot ID`),]
  ## format the Id map
  id_map$Specimen.Label <- gsub("\\s+","",id_map$`Aliquot ID`)
  
  ## add ID mapping info to the clinical file
  tmp <- data.frame(Participant.ID = id_map$ParticipantID, Specimen.Label = id_map$Specimen.Label, tumor_normal = id_map$Group,
                    protein = NA, phosphoprotein = NA, 
                    Parent.Specimen.Label = id_map$Specimen.Label, Parent.Specimen.Title = id_map$Group,
                    Specimen.Comments = NA, Pooled. = NA, Specimen.Class = NA, Specimen.Type = NA, cancer = cancer)
  clinical <- rbind(tmp, clinical[!(clinical$Specimen.Label %in% tmp$Specimen.Label),])
  clinical <- unique(clinical)
  rm(tmp)
  
  # keep shared peptide log ratio
  Pro.f = unshared_pro(Pro)
  Pro.f2p <- data.frame(Gene = rownames(Pro.f))
  Pro.f2p <- cbind(Pro.f2p, Pro.f)

  # normalize by sample
  Pro.n = normalize_by_sample(Pro.f)
  Pro.n2p = data.frame(Gene = rownames(Pro.n))
  Pro.n2p <- cbind(Pro.n2p, Pro.n)
  
  ## annotate clinical file with whether seen in proteomics data
  clinical <- anno_clinical(Pro.n, clinical, "protein")
  
  # average the biological replicates
  Pro.n_aveRep <- aveRep(Pro.n, clinical)
  Pro.n_aveRep2p <- data.frame(Gene = rownames(Pro.n_aveRep))
  Pro.n_aveRep2p <- cbind(Pro.n_aveRep2p, Pro.n_aveRep)
  
  ## annotate clinical file with whether seen in proteomics data
  clinical <- anno_clinical(Pro.n_aveRep2p, clinical, "protein")
  
  # keep tumor only samples and exclude duplicates
  Pro.n_tumor <- get_tumor(Pro.n_aveRep, clinical)
  Pro.n_tumor2p <- print_format_pro(Pro.n_tumor)

  Pro.n_normal <- get_normal(Pro.n_aveRep, clinical)
  Pro.n_normal2p <- print_format_pro(Pro.n_normal)

  # exclude genes with < 10 non-NA samples
  Pro.n_u_max10NA = Pro.n_tumor[rowSums(!is.na(Pro.n_tumor))>=10,]
  Pro.n_u_max10NA2p <- print_format_pro(Pro.n_u_max10NA)

  ## process phosphoproteome
  print(paste0('Start processing ', cancer, ' phosphoproteome data'))
  Pho <- fread(raw[cancer, "phosphoprotein"], data.table=FALSE)
  Pho.f2p = unshared_pho(Pho, Pho[,c("Gene","Phosphosite")])
  Pho.n = normalize_by_sample(Pho.f2p[,!(colnames(Pho.f2p) %in% c("Gene", "Phosphosite"))])
  Pho.n2p <- Pho[,c("Gene","Phosphosite")]
  Pho.n2p <- cbind(Pho.n2p, Pho.n)

  # annotate clinical data
  clinical <- anno_clinical(Pho.n, clinical, "phosphoprotein")

  # average the biological replicates
  Pho.n_aveRep <- aveRep(Pho.n, clinical)
  Pho.n_aveRep2p <- Pho[,c("Gene","Phosphosite")]; Pho.n_aveRep2p <- cbind(Pho.n_aveRep2p, Pho.n_aveRep)
  
  # annotate clinical data
  clinical <- anno_clinical(Pho.n_aveRep, clinical, "phosphoprotein")

  # keep tumor only samples
  Pho.n_tumor <- get_tumor(Pho.n_aveRep, clinical)
  Pho.n_tumor2p <- print_format_pho(Pho.n_tumor, Pho[,c("Gene","Phosphosite")])

  Pho.n_normal <- get_normal(Pho.n_aveRep, clinical)
  Pho.n_normal2p <- print_format_pho(Pho.n_normal, Pho[,c("Gene","Phosphosite")])

  max10NA <- rowSums(!is.na(Pho.n_tumor))>=10
  Pho.n_u_max10NA <- Pho.n_tumor[max10NA,]
  Pho.n_u_max10NA2p <- print_format_pho(Pho.n_u_max10NA, Pho[max10NA,c("Gene","Phosphosite")])

  # collect statistics
  sum_table[prefix[cancer],] <- c(ncol(Pro.f), ncol(Pro.n_tumor), ncol(Pro.n_normal),
                                  nrow(Pro.n), nrow(Pro.n_u_max10NA),
                                  ncol(Pho.n), ncol(Pho.n_tumor), ncol(Pho.n_normal),
                                  length(unique(Pho$Gene)), length(unique(Pho.n_u_max10NA2p$Gene)),
                                  nrow(Pho.n), nrow(Pho.n_u_max10NA))
  
  ## write out tables
  write.table(Pro.f2p, row.names = T, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PRO_formatted.txt",sep=""))
  write.table(Pro.n2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized.txt",sep=""))
  write.table(Pro.n_aveRep2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_replicate_averaged.txt",sep=""))
  write.table(Pro.n_tumor2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_noControl.txt",sep=""))
  write.table(Pro.n_normal2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_Control.txt",sep=""))
  write.table(Pro.n_u_max10NA2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_noControl_max10NA.txt",sep=""))
  write.table(Pho.f2p, quote=F, sep = '\t', row.names = F, file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted.txt", sep=""))
  write.table(Pho.n2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized.txt",sep=""))
  write.table(Pho.n_aveRep2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_replicate_averaged.txt",sep=""))
  write.table(Pho.n_tumor2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_noControl.txt",sep=""))
  write.table(Pho.n_normal2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_Control.txt",sep=""))
  write.table(Pho.n_u_max10NA2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_noControl_max10NA.txt",sep=""))
  sink()
  closeAllConnections()
}
clinical$Specimen.Comments <- NULL
write.table(x = clinical, file = paste0(cptac_sharedD, "Specimen_Data_20161005_Yige_20180815.txt"), quote = F, sep = "\t", row.names = F)
write.table(x = sum_table, file = paste0(cptac_sharedD, "sample_genes_phosphosites_summary_CPTAC3.csv"), quote = F, sep = ",", row.names = T)


# process proteome and phosphoproteome data for CCRCC---------------------------------------------
base_CCRCC_jhu <- paste(cptac_sharedD, "CCRCC", sep = "")
raw["CCRCC","protein"] <- paste(base_CCRCC_jhu,"/CCRCC_batch1-2_SummaryReports/CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.tmt10.tsv", sep = "")
raw["CCRCC", "phosphoprotein"] <- paste(base_CCRCC_jhu,"/CCRCC_batch1-2_SummaryReports/CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Phosphoproteome.phosphosite.tmt10.tsv", sep = "")

## input ID mapping data from supporting documents
id_map_raw <- read_excel("~/Box Sync/cptac2p/cptac_shared/6_CPTAC3_Clear_Cell_Renal_Cell_Carcinoma/CCRCC_batch1_SupportingDocuments/CCRCC_TMT10_Label_to_Sample_Mapping_File_JHU_20171130_edkak.xlsx")
## reshape the id map
channels <- unique(str_split_fixed(colnames(id_map_raw)[grepl("TMT10", colnames(id_map_raw))], pattern = " ", 2)[,1])
partIDs <- sapply(channels, FUN = function(c, map) as.vector(map[, paste0(c, " Participant ID")]), map = id_map_raw)
sampIDs <- sapply(channels, FUN = function(c, map) as.vector(map[, paste0(c, " Specimen Label")]), map = id_map_raw)
id_map <- data.frame(ParticipantID = unlist(partIDs), Specimen.Label = unlist(sampIDs))
id_map <- unique(id_map[!is.na(id_map$ParticipantID),])
id_map$Group <- "other"
id_map$Group[grepl("CPT", id_map$Specimen.Label)] <- ifelse(substr(x = id_map$Specimen.Label[grepl("CPT", id_map$Specimen.Label)], start = 12, stop = 13) == "01", "Normal", "Tumor")

for (cancer in c("CCRCC")) {
  sink(paste0(cptac_sharedD, cancer, "/", 'preprocessing_log.txt'), append=FALSE, split=FALSE)
  
  ## process proteome
  print(paste0('-Start processing ', cancer, ' proteome data'))
  
  ## input proteome data
  Pro <- fread(raw[cancer, "protein"], data.table=FALSE)
  
  ## add ID mapping info to the clinical file
  # tmp <- data.frame(Participant.ID = id_map$ParticipantID, Specimen.Label = id_map$Specimen.Label, tumor_normal = id_map$Group,
  #                   protein = NA, phosphoprotein = NA, 
  #                   Parent.Specimen.Label = id_map$Specimen.Label, Parent.Specimen.Title = id_map$Group,
  #                   Specimen.Comments = NA, Pooled. = NA, Specimen.Class = NA, Specimen.Type = NA)
  # clinical <- rbind(clinical, tmp)
  # rm(tmp)
  
  # keep shared peptide log ratio
  Pro.f = unshared_pro(Pro)
  Pro.f2p <- data.frame(Gene = rownames(Pro.f))
  Pro.f2p <- cbind(Pro.f2p, Pro.f)
  
  # normalize by sample
  Pro.n = normalize_by_sample(Pro.f)
  Pro.n2p = data.frame(Gene = rownames(Pro.n))
  Pro.n2p <- cbind(Pro.n2p, Pro.n)
  
  ## annotate clinical file with whether seen in proteomics data
  # clinical <- anno_clinical(Pro.n, clinical, "protein")
  
  # average the biological replicates
  # Pro.n_aveRep <- aveRep(Pro.n, clinical)
  Pro.n_aveRep <- Pro.n
  Pro.n_aveRep2p <- data.frame(Gene = rownames(Pro.n_aveRep))
  Pro.n_aveRep2p <- cbind(Pro.n_aveRep2p, Pro.n_aveRep)
  
  # ## annotate clinical file with whether seen in proteomics data
  # clinical <- anno_clinical(Pro.n_aveRep2p, clinical, "protein")
  
  # keep tumor only samples and exclude duplicates
  # Pro.n_tumor <- get_tumor(Pro.n_aveRep, clinical)
  Pro.n_tumor <- Pro.n_aveRep[, substr(colnames(Pro.n_aveRep), start = 12, stop = 12) == "0" & substr(colnames(Pro.n_aveRep), start = 13, stop = 13) != "1"]
  Pro.n_tumor2p <- print_format_pro(Pro.n_tumor)
  
  # Pro.n_normal <- get_normal(Pro.n_aveRep, clinical)
  Pro.n_normal <- Pro.n_aveRep[, substr(colnames(Pro.n_aveRep), start = 12, stop = 12) == "0" & substr(colnames(Pro.n_aveRep), start = 13, stop = 13) == "1"]
  Pro.n_normal2p <- print_format_pro(Pro.n_normal)
  
  # exclude genes with < 10 non-NA samples
  Pro.n_u_max10NA = Pro.n_tumor[rowSums(!is.na(Pro.n_tumor))>=10,]
  Pro.n_u_max10NA2p <- print_format_pro(Pro.n_u_max10NA)
  
  ## process phosphoproteome
  print(paste0('Start processing ', cancer, ' phosphoproteome data'))
  Pho <- fread(raw[cancer, "phosphoprotein"], data.table=FALSE)
  Pho.f2p = unshared_pho(Pho, Pho[,c("Gene","Phosphosite")])
  Pho.n = normalize_by_sample(Pho.f2p[,!(colnames(Pho.f2p) %in% c("Gene", "Phosphosite"))])
  Pho.n2p <- Pho[,c("Gene","Phosphosite")]
  Pho.n2p <- cbind(Pho.n2p, Pho.n)
  
  # annotate clinical data
  # clinical <- anno_clinical(Pho.n, clinical, "phosphoprotein")
  
  # average the biological replicates
  # Pho.n_aveRep <- aveRep(Pho.n, clinical)
  Pho.n_aveRep <- Pho.n
  Pho.n_aveRep2p <- Pho[,c("Gene","Phosphosite")]; Pho.n_aveRep2p <- cbind(Pho.n_aveRep2p, Pho.n_aveRep)
  
  # annotate clinical data
  # clinical <- anno_clinical(Pho.n_aveRep, clinical, "phosphoprotein")
  
  # keep tumor only samples
  # Pho.n_tumor <- get_tumor(Pho.n_aveRep, clinical)
  Pho.n_tumor <- Pho.n_aveRep[, substr(colnames(Pho.n_aveRep), start = 12, stop = 12) == "0" & substr(colnames(Pho.n_aveRep), start = 13, stop = 13) != "1"]
  Pho.n_tumor2p <- print_format_pho(Pho.n_tumor, Pho[,c("Gene","Phosphosite")])
  
  # Pho.n_normal <- get_normal(Pho.n_aveRep, clinical)
  Pho.n_normal <- Pho.n_aveRep[, substr(colnames(Pho.n_aveRep), start = 12, stop = 12) == "0" & substr(colnames(Pho.n_aveRep), start = 13, stop = 13) == "1"]
  Pho.n_normal2p <- print_format_pho(Pho.n_normal, Pho[,c("Gene","Phosphosite")])
  
  max10NA <- (rowSums(!is.na(Pho.n_tumor))>=10)
  Pho.n_u_max10NA <- Pho.n_tumor[max10NA,]
  Pho.n_u_max10NA2p <- print_format_pho(Pho.n_u_max10NA, Pho[max10NA,c("Gene","Phosphosite")])
  
  # collect statistics
  sum_table[prefix[cancer],] <- c(ncol(Pro.f), ncol(Pro.n_tumor), ncol(Pro.n_normal),
                                  nrow(Pro.n), nrow(Pro.n_u_max10NA),
                                  ncol(Pho.n), ncol(Pho.n_tumor), ncol(Pho.n_normal),
                                  length(unique(Pho$Gene)), length(unique(Pho.n_u_max10NA2p$Gene)),
                                  nrow(Pho.n), nrow(Pho.n_u_max10NA))
  
  ## write out tables
  write.table(Pro.f2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PRO_formatted.txt",sep=""))
  write.table(Pro.n2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized.txt",sep=""))
  write.table(Pro.n_aveRep2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_replicate_averaged.txt",sep=""))
  write.table(Pro.n_tumor2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_noControl.txt",sep=""))
  write.table(Pro.n_normal2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_Control.txt",sep=""))
  write.table(Pro.n_u_max10NA2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_noControl_max10NA.txt",sep=""))
  write.table(Pho.f2p, quote=F, sep = '\t', row.names = F, file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted.txt", sep=""))
  write.table(Pho.n2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized.txt",sep=""))
  write.table(Pho.n_aveRep2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_replicate_averaged.txt",sep=""))
  write.table(Pho.n_tumor2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_noControl.txt",sep=""))
  write.table(Pho.n_normal2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_Control.txt",sep=""))
  write.table(Pho.n_u_max10NA2p, row.names = F, quote=F, sep = '\t', file=paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_noControl_max10NA.txt",sep=""))
  sink()
  closeAllConnections()
}

# write out summary table -------------------------------------------------
# write.csv(sum_table, paste0(cptac_sharedD, "sample_genes_phosphosites_summary.csv"))
write_delim(clinical, delim = "\t", paste0(cptac_sharedD, "Specimen_Data_20161005_Yige_20180702.txt"))


