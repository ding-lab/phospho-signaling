# Yige Wu @ WashU 2018 Aug
## For formating LUAD discovery manuscript data in data freeze to be in concord with CDAP proteomics data and genomics data

# source ------------------------------------------------------------------
code_top_dir <- "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/"
path2phospho_network_shared <- paste0(code_top_dir, "phospho_network/phospho_network_shared.R")
source(path2phospho_network_shared)


# Install pkgs ------------------------------------------------------------
# for (pkg_name_tmp in c("phantasus")) {
#   if (!(pkg_name_tmp %in% installed.packages()[,1])) {
#     if (!requireNamespace("BiocManager", quietly = TRUE))
#       install.packages("BiocManager")
#     
#     BiocManager::install("phantasus")
#   }
# }
# library(phantasus)

# set variables -----------------------------------------------------------
cancer <- "LUAD"
pipeline_type <- "PGDAC"

# rewrite protein data (Umich) ------------------------------------------------------
pro_data <- fread(input = "Ding_Lab/Projects_Current/CPTAC/PGDAC/CPTAC3.0_LUAD_Data/formatted/20190628/prot.LUAD.allgene.tsv", data.table = F)
write.table(x = pro_data, file = paste0(makeOutDir(), "LUAD_PRO_Tumor_PGDAC_SD_partID.txt"), quote = F, row.names = F, sep = "\t")

# 
# # Rewrite phosphorylation gene level --------------------------------------
# phog_data <- fread(input = "Ding_Lab/Projects_Current/CPTAC/PGDAC/JHU_PCC_shared_data/LUAD_analysis/proteomics/phospho/Michigan/P_Batch1-5/6_CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_P_abundance_groupby=0_protNorm=2.tsv", data.table = F)
# phog_data <- data.frame(phog_data)
# ncol(phog_data)
# head(phog_data)
# samples <- colnames(phog_data)[!(colnames(phog_data) %in% c("Index", "NumberPSM", "Proteins", "ReferenceIntensity"))]
# pro_values.a <- phog_data[, samples]
# pro_values.r <- phog_data[, "ReferenceIntensity"]
# pro_values.a_r <- pro_values.a - pro_values.r
# phog_data <- cbind(data.frame(Gene = phog_data$Index), pro_values.a_r)
# write.table(x = phog_data, file = paste0(makeOutDir(resultD = resultD), "6_CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_P_abundance_groupby=0_protNorm=2.tsv.formatted.txt"), quote = F, row.names = F, sep = "\t")

# rewrite phosphorylation data ----------------------------------------------
pho_data <- fread(input = "Ding_Lab/Projects_Current/CPTAC/PGDAC/CPTAC3.0_LUAD_Data/formatted/20190628/phos.LUAD.allgene.tsv", data.table = F)
pho_data <- data.frame(pho_data)
samples <- colnames(pho_data)[!(colnames(pho_data) %in% c("Gene"))]
pho_values <- pho_data[, samples]
pho_data <- cbind(data.frame(Gene = str_split_fixed(string = pho_data$Gene, pattern = "[-stym]", n = 3)[,1], 
                             Phosphosite = str_split_fixed(string = pho_data$Gene, pattern = "[-stym]", n = 3)[,2]), pho_values)
pho_data <- pho_data[pho_data$Phosphosite != "",]
write.table(x = pho_data, file = paste0(makeOutDir(), "LUAD_PHO_Tumor_PGDAC_SD_partID.txt"), quote = F, row.names = F, sep = "\t")

# input sample mapping file -----------------------------------------------
# sample_map <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/PGDAC/A6_PGDAC_Data_Clear_Cell_Renal_Cell_Carcinoma/MountSinai/CPTAC3_LUAD/Batch1-4/CPTAC-3 LUAD samples tumor-normal info.csv", data.table = F)
# sample_map <- unique(sample_map[, c("Aliquot ID", "Case ID", "Tissue Type")])
# rownames(sample_map) <- sample_map$`Aliquot ID`
# sample_map <- sample_map[!(sample_map$`Case ID` %in% c("C3N-01175", "C3N-01180", "C3N-00313", "C3N-00435", "C3N-00832")),]
# 
# # rewrite proteomics/phosphoproteomics data split by tumor and normal ----------------------------------
# for (sample_type in c("tumor", "normal")) {
#   sampIDs <- sample_map$`Aliquot ID`[grepl(x = sample_map$`Tissue Type`, pattern= sample_type, ignore.case = T)]
#   sampIDs <- intersect(sampIDs, samples)
#   partIDs <- sample_map[sampIDs, "Case ID"]
#   
#   ## write protein data
#   pro2w <- pro_data[, c("Gene", sampIDs)]
#   pro2w.partID <- pro2w
#   colnames(pro2w.partID) <- c("Gene", partIDs)
#   expresson_type <- "PRO"
#   write.table(x = pro2w.partID, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "MD_MAD", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
#   
#   ## write gene-level phosphorylation data
#   phog2w <- phog_data[, c("Gene", sampIDs)]
#   phog2w.partID <- phog2w
#   colnames(phog2w.partID) <- c("Gene", partIDs)
#   expresson_type <- "collapsed_PHO"
#   write.table(x = phog2w.partID, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "MD_MAD", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
#   
#   
#   ## write phosphorylatin data
#   pho2w <- pho_data[, c("Gene", "Phosphosite", sampIDs)]
#   pho2w.partID <- pho2w
#   colnames(pho2w.partID) <- c("Gene", "Phosphosite", partIDs)
#   expresson_type <- "PHO"
#   write.table(x = pho2w.partID, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "MD_MAD", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
# }
# 
# # Parse Michigan’s CNV ------------------------------------------------------
# cna <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/PGDAC/JHU_PCC_shared_data/LUAD_analysis/genomics/WXS_CNV/michigan/kirc_cnv_genelr.csv", data.table = F)
# cna <- cna[, c("gene_name", colnames(cna)[grepl(pattern = "lr_", x = colnames(cna))])]
# colnames(cna) <- str_split_fixed(string = colnames(cna), pattern = "lr_", 2)[,2]
# colnames(cna)[1] <- "gene"
# cna.val <- as.matrix(cna[,-1])
# cna.val[is.infinite(cna.val)] <- NA
# colnames_tmp <- colnames(cna)
# cna <- cbind(data.frame(gene = cna$gene), data.frame(cna.val))
# colnames(cna) <- colnames_tmp
# cancer <- "LUAD"
# write.table(x = cna, file = paste0(makeOutDir(resultD = resultD), "somatic_CNA", ".", cancer, ".partID.txt"), quote = F, row.names = F, sep = "\t")
# 
# 
