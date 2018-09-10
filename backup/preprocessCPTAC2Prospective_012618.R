### preprocessCPTAC2Prospective.R ###
# Yige Wu @ WashU Dec 2017 
# preprocess the proteome and phosphoproteome data

# input -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
clinical <- read.delim(paste0(baseD, "Specimen_Data_20161005_Yige_20180126.txt"), sep = "\t")
clinical <- data.frame(clinical)
clinical_anno <- clinical

# initiate summary table --------------------------------------------------
sum_table <- matrix(NA, nrow = 3, ncol = 12, 
                    dimnames = list(row=prefix, 
                                    col=c("samples_all_pro","samples_tumor_pro", "samples_normal_pro",
                                          "genes_pro","genes_10NonNA_pro",
                                          "samples_all_pho","samples_tumor_pho", "samples_normal_pho",
                                          "genes_pho","genes_10NonNA_pho",
                                          "phosphosites","phosphosites_10NonNA")))


# process proteome and phosphoproteome data ---------------------------------------------
for (cancer in c("BRCA", "OV", "CO")) {
  ## process proteome
  Pro <- fread(raw[cancer, "protein"], data.table=FALSE)

  # keep shared peptide log ratio
  Pro.f = unshared_pro(Pro)
  
  # annotate clinical data
  clinical <- anno_clinical(Pro.f, clinical, "protein")

  # normalize by sample
  Pro.n = normalize_by_sample(Pro.f)
  
  # average the biological replicates
  Pro.n_aveRep <- aveRep(Pro.n, clinical)
  
  # keep tumor only samples and exclude duplicates
  Pro.n_tumor <- get_tumor(Pro.n_aveRep, clinical)
  Pro.n_tumor2p <- print_format_pro(Pro.n_tumor)
  
  Pro.n_normal <- get_normal(Pro.n_aveRep, clinical)
  Pro.n_normal2p <- print_format_pro(Pro.n_normal)
  
  # exclude genes with > 10 NA samples
  Pro.n_u_max10NA = Pro.n_tumor[rowSums(!is.na(Pro.n_tumor))>=10,]
  Pro.n_u_max10NA2p <- print_format_pro(Pro.n_u_max10NA)
  
  ## process phosphoproteome
  Pho <- fread(raw[cancer, "phosphoprotein"], data.table=FALSE)
  Pho.f2p = unshared_pho(Pho, Pho[,c("Gene","Phosphosite")])
  Pho.n = normalize_by_sample(Pho.f2p[,!(colnames(Pho.f2p) %in% c("Gene", "Phosphosite"))])
  Pho.n2p <- Pho[,c("Gene","Phosphosite")]
  Pho.n2p <- cbind(Pho.n2p, Pho.n)
  
  # annotate clinical data
  clinical <- anno_clinical(Pho.n, clinical, "phosphoprotein")
  
  # keep tumor only samples and average the duplicates
  Pho.n_tumor <- tumor_aveRep(Pho.n, clinical)
  Pho.n_tumor2p <- print_format_pho(Pho.n_tumor, Pho[,c("Gene","Phosphosite")])
  
  Pho.n_normal <- normal_aveRep(Pho.n, clinical)
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
  write.table(Pro.f, row.names = T, quote=F, sep = '\t', file=paste(baseD, cancer, "/", prefix[cancer], "_PRO_formatted.txt",sep=""))
  write.table(Pro.n, row.names = T, quote=F, sep = '\t', file=paste(baseD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized.txt",sep=""))
  write.table(Pro.n_tumor2p, row.names = F, quote=F, sep = '\t', file=paste(baseD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_noControl.txt",sep=""))
  write.table(Pro.n_normal2p, row.names = F, quote=F, sep = '\t', file=paste(baseD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_Control.txt",sep=""))
  write.table(Pro.n_u_max10NA2p, row.names = F, quote=F, sep = '\t', file=paste(baseD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_noControl_max10NA.txt",sep=""))
  write.table(Pho.f2p, quote=F, sep = '\t', row.names = F, file=paste(baseD, cancer, "/", prefix[cancer], "_PHO_formatted.txt", sep=""))
  write.table(Pho.n2p, row.names = F, quote=F, sep = '\t', file=paste(baseD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized.txt",sep=""))
  write.table(Pho.n_tumor2p, row.names = F, quote=F, sep = '\t', file=paste(baseD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_noControl.txt",sep=""))
  write.table(Pho.n_normal2p, row.names = F, quote=F, sep = '\t', file=paste(baseD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_Control.txt",sep=""))
  write.table(Pho.n_u_max10NA2p, row.names = F, quote=F, sep = '\t', file=paste(baseD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_noControl_max10NA.txt",sep=""))

}

write.csv(sum_table, paste0(baseD, "sample_genes_phosphosites_summary.csv"))
write_delim(clinical, delim = "\t", paste0(baseD, "Specimen_Data_20161005_Yige_20171214.xls"))
write_delim(clinical[!is.na(clinical$protein) & clinical$protein & clinical$tumor_normal=="Tumor",], delim = "\t", paste0(baseD, "sample_mapping"))


# remember to check no duplicated samples from the same participant

