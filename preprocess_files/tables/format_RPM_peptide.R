# Yige Wu @ WashU 2018 Mar
## format the prospective BRCA PRM peptide-level abundance

# souce -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
resultDnow <- makeOutDir()

# inputs ------------------------------------------------------------------
## checked the proteins by HUMAN
prm_raw <- read_excel("~/Box Sync/cptac2p/reid/CPTAC 007 PRM - Peptide Level Results - fmol Kinase.xlsx")
partID_map <- read_excel("~/Box Sync/cptac2p/reid/Sample list_CPTAC_007.xlsx", skip = 1)

## take only mean value columns
prm_means <- prm_raw[, grep("Average", colnames(prm_raw), ignore.case = T)]
prm_means <- as.matrix(prm_means)

## transform triplicate ID to sample ID
triplicate_names <- colnames(prm_means)
triplicate_start <- str_split_fixed(triplicate_names, pattern = " |-", 4)[,3]
partID_map_start <- partID_map[partID_map$`AS No.` %in% triplicate_start,]

## annotate gene names from the head columns
prm_head <- prm_raw[,c("Protein Name","Peptide Modified Sequence","#")]
prm_head <- data.frame(prm_head)
prm_head$index <- 1:nrow(prm_raw)
tmp <- str_split_fixed(string = prm_head$Protein.Name, pattern = "\\|", 3)
tmp1 <- tmp[,1]
tmp2 <- tmp[,3]
### add in the gene names already there
tmp3 <- str_split_fixed(string = tmp2, pattern = " ", 2)
tmp4 <- str_split_fixed(string = tmp3[,2], pattern = "GN=", 2)
tmp5 <- str_split_fixed(string = tmp4[,2], pattern = " PE", 2)

### transfrom the rest
tmp11 <- tmp[tmp[,3] == "" & grepl(pattern = "HUMAN", x = tmp[,1]),1]
tmp12 <- tmp[tmp3[,2] == "" & grepl(pattern = "HUMAN", x = tmp[,3]),3]
genenames_orig <- data.frame(From = unique(c(tmp11, tmp12)))
tmp13 <- tmp[tmp[,3] == "" & !grepl(pattern = "HUMAN", x = tmp[,1]),1]

# resultDnow <- makeOutDir()
# write.table(genenames_orig, file = paste0(resultDnow, "uniprot2genename.txt"), quote = F, col.names = F, row.names = F)

uniprot2genename_mapped <- read_delim("~/Box Sync/cptac2p/cptac_shared/analysis_results/preprocess_files/format_RPM_peptide/uniprot2genename_mapped.txt",
                                      "\t", escape_double = FALSE, trim_ws = TRUE)
uniprot2genename_mapped <- data.frame(uniprot2genename_mapped)
uniprot2genename_notmapped <- read_delim("~/Box Sync/cptac2p/cptac_shared/analysis_results/preprocess_files/format_RPM_peptide/uniprot2genename_notmapped.txt",
                                      "\t", escape_double = FALSE, trim_ws = TRUE)
uniprot2genename_notmapped <- data.frame(From = uniprot2genename_notmapped$`not mapped`, To = c("PFKM", "PSPH"))
proteinname2genename <- data.frame(From = unique(tmp13), To = c("HIST1H2AG", "HIST1H2BA", "HIST1H3A", "HIST1H4A"))
genenames <- rbind(uniprot2genename_mapped, uniprot2genename_notmapped, proteinname2genename)
genenames_orig <- merge(genenames_orig, genenames, all.x = T, sort = F)

prm_head$Gene_Name[tmp5[,1] != ""] = tmp5[tmp5[,1] != "",1]
prm_head$From <- NA
prm_head$From[tmp[,3] == "" & grepl(pattern = "HUMAN", x = tmp[,1])] <- tmp11
prm_head$From[tmp3[,2] == "" & grepl(pattern = "HUMAN", x = tmp[,3])] <- tmp12
prm_head$From[tmp[,3] == "" & !grepl(pattern = "HUMAN", x = tmp[,1])] <- tmp13
prm_head <- merge(prm_head, genenames, all.x = T, sort = F)
prm_head$Gene_Name[is.na(prm_head$Gene_Name)] = prm_head$To[is.na(prm_head$Gene_Name)]
prm_head <- prm_head[order(prm_head$index),]
prm_raw_nohead <- data.frame(prm_raw)
prm_raw_nohead <- data.frame(prm_raw_nohead[, !(colnames(prm_raw_nohead) %in% colnames(prm_head))])
prm_raw.f <- cbind(prm_head[, c("Protein.Name","Peptide.Modified.Sequence", "X.", "Gene_Name")], prm_raw_nohead)
prm_raw.f <- data.frame(prm_raw.f)
write_delim(prm_raw.f, path = paste0(makeOutDir(), "prm_genename_annotated.txt"), delim = "\t")

for (datatype in c("Average", "Stdev", "CV")) {
  prm_raw.f_specimen_labeled <- prm_raw.f[,c("Protein.Name","Peptide.Modified.Sequence", "X.", "Gene_Name", colnames(prm_raw.f)[grepl(pattern = datatype, x = colnames(prm_raw.f))])]
  colnames(prm_raw.f_specimen_labeled) <- c("Protein.Name","Peptide.Modified.Sequence", "X.", "Gene_Name", as.vector(partID_map_start$Label))
  write_delim(prm_raw.f_specimen_labeled, path = paste0(makeOutDir(), "prm_", datatype,"_genename_annotated_specimen_labeled.txt"), delim = "\t")
  
  prm_raw.f_PPI_labeled <- prm_raw.f[,c("Protein.Name","Peptide.Modified.Sequence", "X.", "Gene_Name", colnames(prm_raw.f)[grepl(pattern = datatype, x = colnames(prm_raw.f))])]
  colnames(prm_raw.f_PPI_labeled) <- c("Protein.Name","Peptide.Modified.Sequence", "X.", "Gene_Name", as.vector(partID_map_start$PPI))
  write_delim(prm_raw.f_PPI_labeled, path = paste0(makeOutDir(), "prm_", datatype,"_genename_annotated_PPI_labeled.txt"), delim = "\t")
  
}

## take the median value of mean values of all peptides for each protein
prm_mean_median <- matrix(data = NA, nrow = length(unique(prm_head$Gene_Name)), ncol = ncol(prm_means))
rownames(prm_mean_median) <- unique(prm_head$Gene_Name)
for (g in unique(prm_head$Gene_Name)) {
  for (c in 1:ncol(prm_means)) {
    prm_mean_median[g, c] <- median(as.numeric(as.vector(prm_means[prm_head$Gene_Name == g,c])), na.rm = T)
  }
}

prm_mean_median <- data.frame(prm_mean_median)
colnames(prm_mean_median) <- partID_map_start$PPI
prm_mean_median$Gene_Name <- rownames(prm_mean_median)
prm_mean_median <- prm_mean_median[, c("Gene_Name", partID_map_start$PPI)]
write.table(prm_mean_median, file = paste0(resultDnow, "prm_kinase_median_per_triplicate_PPI_labeled.txt"), quote = F, col.names = T, row.names = F)

colnames(prm_mean_median) <- c("Gene_Name", partID_map_start$Label)
write.table(prm_mean_median, file = paste0(resultDnow, "prm_kinase_median_per_triplicate_specimen_labeled.txt"), quote = F, col.names = T, row.names = F)

