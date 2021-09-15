# Yige Wu @ WashU 2018 Jan
## shared parameters and functions for cptac2p_analysis

#  library -----------------------------------------------------------------
library(data.table)
library(readxl)

# set variables -----------------------------------------------------------
cancers_sort <- c("BRCA", "OV", "CO")
prefix <- NULL
prefix["BRCA"] <- "BRCA_BI"
prefix["OV"] <- "OV_PNNL"
prefix["CO"] <- "CO_PNNL"
prefix["UCEC"] <- "UCEC_PNNL"
prefix["CCRCC"] <- "CCRCC_JHU"
prefix["LUAD"] <- "LUAD_BI"

datatypes <- c("PRO", "PHO", "collapsed_PHO")

cancer_full_names <- c("Breast Cancer", "Ovarian Cancer", "Colorectal Cancer"); names(cancer_full_names) <- cancers_sort

# heavily-used data paths ---------------------------------------------------------
## base directory at Box Sync
baseD = "~/Box/"
setwd(baseD)
dir2dinglab_projects <- paste0(baseD, "Ding_Lab/Projects_Current/")

## CPTAC2 prospective and CPTAC proteomics data
dir2cptac_shared <- paste0(dir2dinglab_projects, "CPTAC/cptac_shared/")
dir2cptac_pgdac <- paste0(dir2dinglab_projects, "CPTAC/PGDAC/")

## CPTAC2 prospective genomics data
dir2cptac2_prospective <- paste0(dir2dinglab_projects, "CPTAC/CPTAC_Prospective_Samples/")

## CPTAC2 retrospective proteomics/genomics data
dir2cptac2_retrospective <- paste0(dir2dinglab_projects, "CPTAC/pan3can_shared_data/")

## processed data directory
dir2phospho_signaling <- paste0(dir2dinglab_projects, "PanCan_Phospho-signaling/")
resultD <- paste0(dir2phospho_signaling, "analysis_results/")
preprocess_filesD <- paste0(resultD, "preprocess_files/")

## the path to scripts relating to phospho-signaling
ppnD <- paste0(resultD, "phospho_network/")

maf_files <- list.files(path = paste0(dir2cptac2_prospective, "Somatic/"))

#list.files(path = paste0(dir2cptac2_prospective, "FPKM/08092018/"), recursive = T)
fpkm_files <- paste0(dir2cptac2_prospective, "FPKM/08092018/", 
                     c("fpkm_prospective_breastcancer_110617_swap_corrected.csv",
                       "fpkm_ov_allgenes.011718.csv",
                       "fpkm_crc_allgenes.v1.0.10232017.csv"))
names(fpkm_files) <- cancers_sort



# gene lists --------------------------------------------------------------
SWI_genes <- c("PBRM1", "ARID1A", "SMARCA4", "ARID1B", "ARID4A", "SMARCA5", "VHL")

## SMGs
SMGs <- list()
SMGs[["BRCA"]] <- c("PIK3CA", "TP53", "MAP3K1", "MAP2K4", "GATA3", "MLL3", "CDH1", "PTEN", "PIK3R1", "AKT1", "RUNX1", "CBFB", "TBX3", "NCOR1", "CTCF", "FOXA1", "SF3B1", "CDKN1B", "RB1", "AFF2", "NF1", "PTPN22", "PTPRD")
SMGs[["OV"]] <- c("TP53", "NF1", "BRCA1", "BRCA2", "RB1", "CDK12")
SMGs[["CO"]] <- unique(c("APC", "TP53", "KRAS", "PIK3CA", "FBXW7", "SMAD4", "TCF7L2", "NRAS", "BRAF", "CTNNB1", "SMAD2", "FAM123B", "SOX9", "ATM", "ARID1A", "ACVR2A", "APC", "TGFBR2", "MSH3", "MSH6", "SLC9A9", "TCF7L2"))
SMGs[["UCEC"]] <- c("FLNA",
  "MAP3K4",
  "HUWE1",
  "NSD1",
  "JAK1",
  "RPL22",
  "SCAF4",
  "FBXW7",
  "KMT2D",
  "INPPL1",
  "ZFHX3",
  "KMT2B",
  "TP53",
  "CTCF",
  "CTNNB1",
  "KRAS",
  "PIK3R1",
  "ARID1A",
  "PIK3CA",
  "PTEN")
SMGs[["CCRCC"]] <- c("VHL", "PBRM1", "SETD2", "KDM5C", "PTEN", "BAP1", "MTOR", "TP53")
SMGs[["LIHC"]] <- c("TP53", "AXIN1", "RB1", "CTNNB1", "ARID1A", "ARID2", "BAP1", "NFE2L2", "KEAP1", "ALB", "APOB")
SMGs[["LUAD"]] <- c("KRAS", "EGFR", "BRAF", "ROS1", "ALK", "RET", "MAP2K1", "HRAS", "NRAS", "MET", "ERBB2", "RIT1", "NF1")
SMGs[["GBM"]] <- c(
  # RTK genes
  'EGFR', # 'ERBB2', 'ERBB3', 'ERBB4',
  'PDGFRA', 'PDGFRB',
  'MET',
  # 'FGFR1', 'FGFR2', 'FGFR3',
  # PI3K genes
  'PIK3CA', 'PIK3R1',
  # 'PIK3CG', 'PIK3C2G',
  # 'PIK3CB', 'PIK3C2B', 'PIK3C2A', 'PIK3R2',
  'PTEN',
  # MAPK genes
  'NF1',
  'BRAF',
  # TP53 genes
  'TP53',
  'MDM2', 'MDM4', 'MDM1', 'PPM1D',
  'CDKN2A',
  # RB1 genes
  'RB1',
  'CDKN2C', 'CDKN1A', 'CDKN2B',
  'CDK4',
  'CDK6',
  # Chromosome modification genes
  'IDH1', 'IDH2',
  'ATRX', 'SMARCAL1'
)

## significant SCNA
CNGs <- list()
CNGs[["GBM"]] <- c("EGFR", "PRDM2", "MDM4", "AKT3", "MYCN", "SOX2", "FGFR3", "PDGFRA", "CDK6", "MET", "MYC", "CCND2", "CDK4", "MDM2", "IRS2", "AKT1", "HYDIN", "GRB2", "CCNE1",
                   "CDKN2C", "LSAMP", "QKI", "CDKN2A", "CDKN2B", "PTEN", "RB1", "NPAS3", "TP53", "NF1")


# oncogenes and TSGs ------------------------------------------------------
driver_genes <- read_excel("~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/Resources/Knowledge/Gene_Lists/mmc1.xlsx", 
                           sheet = "Table S1", skip = 3)
driver_genes <- data.frame(driver_genes)
oncogenes <- driver_genes$Gene[grepl(x = driver_genes$Tumor.suppressor.or.oncogene.prediction..by.20.20.., pattern = "oncogene")]
tsgs <- driver_genes$Gene[grepl(x = driver_genes$Tumor.suppressor.or.oncogene.prediction..by.20.20.., pattern = "tsg")]
oncogenes <- c(oncogenes, "CTNND1", "PIK3R1", "INSR", "PPM1D")
tsgs <- tsgs[!(tsgs %in% c("CDH1", "CTNND1", "PIK3R1", "INSR", "PPM1D"))]
oncogenes <- unique(oncogenes)
tsgs <- unique(tsgs)

# functions ---------------------------------------------------------------
loadMSIMap <- function() {
  msi_score <- fread(input = "./Ding_Lab/Projects_Current/MSI/MSI_shared_data/MSI_score/CPTAC.MSI.score.tsv", data.table = F)
  msi_score$MSI_type <- ""
  msi_score$MSI_type[msi_score$Score >= 3.5] <- "MSI-H"
  msi_score$MSI_type[msi_score$Score < 3.5 & msi_score$Score >= 1.0] <- "MSI-L"
  msi_score$MSI_type[msi_score$Score < 1.0] <- "MSS"
  return(msi_score)
}

loadPAM50Map <- function() {
  file <- fread(paste0(dir2cptac_shared, "BRCA/fpkm_pam50_brca_080718_pam50scores_swap_corrected.txt"), data.table = F)
  return(file)
}

partID2pam50 = function(patientID_vector, pam50_map) {
  partIDs <- patientID_vector
  subtypes <- NULL
  for (partID in partIDs) {
    subtype <- unique(as.character(pam50_map$Call[pam50_map$V1 == partID]))
    subtype <- subtype[!is.na(subtype)]
    if (length(subtype) > 0) {
      subtypes <- c(subtypes, subtype)
    } else {
      subtypes <- c(subtypes, "no_RNA-seq")
    }
  }
  names(subtypes) <- partIDs
  return(subtypes)
}

partID2MSI = function(patientID_vector, subtype_map) {
  partIDs <- patientID_vector
  subtypes <- NULL
  for (partID in partIDs) {
    subtype <- unique(as.character(subtype_map$MSI_type[subtype_map$Sample == partID]))
    subtype <- subtype[!is.na(subtype)]
    if (length(subtype) > 0) {
      subtypes <- c(subtypes, subtype)
    } else {
      subtypes <- c(subtypes, "no_genomic_data")
    }
  }
  names(subtypes) <- partIDs
  
  return(subtypes)
}

partID2UCECsubtype = function(patientID_vector) {
  subtype_map <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/supporting_documents/UCEC_CPTAC3_meta_table_subtyped_v2.0_filtered_FIGOed.txt", data.table = F)
  partIDs <- patientID_vector
  subtypes <- NULL
  for (partID in partIDs) {
    subtype <- unique(as.character(subtype_map$Genomic_subtype[subtype_map$Proteomics_Participant_ID == partID]))
    subtype <- subtype[!is.na(subtype)]
    if (length(subtype) > 0) {
      subtypes <- c(subtypes, subtype)
    } else {
      subtypes <- c(subtypes, "no_subtype_info")
    }
  }
  names(subtypes) <- partIDs
  
  return(subtypes)
}

partID2CCRCCImmune = function(patientID_vector) {
  subtype_map <- fread(input = paste0(resultD, "phospho_network/druggability/figures/heatmap_mutation_CNA_RNA_PRO_PHO_cascade_for_ccRCC/CCRCC/CPTAC3_ccRCC_discovery_set_sample_annotation.csv"), data.table = F, sep = ",")
  partIDs <- patientID_vector
  subtypes <- NULL
  for (partID in partIDs) {
    subtype <- unique(as.character(subtype_map$immune_group[subtype_map$partID == partID]))
    subtype <- subtype[!is.na(subtype)]
    if (length(subtype) > 0) {
      subtypes <- c(subtypes, subtype)
    } else {
      subtypes <- c(subtypes, "no_subtype_info")
    }
  }
  names(subtypes) <- partIDs
  
  return(subtypes)
}

# loadpam50fromclinTab <- function() {
#   file <- read_excel(paste0(dir2cptac_shared, "5_CPTAC2_Breast_Prospective_Collection_BI/proteome-data-v1.01011/data/20171106_CPTAC2_ProspectiveBreastCancer_Sample_Annotations_Table_v79.xlsx"))
#   return(file)
# }


# load data functions -----------------------------------------------------
loadSampMap <- function() {
  file <- fread(input = paste0(dir2cptac_shared, "Specimen_Data_20161005_Yige_20180815.txt"), data.table = F, sep = "\t")
  return(file)
}

loadMaf <- function(cancer, maf_files) {
  if (cancer %in% c("BRCA", "OV", "CO")) {
    maf <- fread(input = paste0(dir2cptac2_prospective, "Somatic/", maf_files[grepl(x = maf_files, pattern = cancer)]), data.table = F, fill=TRUE)
  } else if (cancer %in% c("UCEC")) {
    maf <- fread(input = paste0(dir2cptac_pgdac, "Endometrium_CPTAC3/01_Data_tables/Baylor_DataFreeze_V2/UCEC_somatic_mutation_site_level_V2.0.maf"), data.table = F, fill=TRUE)
  } else if (cancer %in% c("CCRCC")) {
    maf <- fread(input = paste0(dir2cptac_pgdac, "ccRCC_discovery_manuscript/ccRCC_expression_matrices/Somatic_Variants/ccrcc.somatic.consensus.gdc.umichigan.wu.112918.maf"), data.table = F, fill=TRUE) 
  } else if (cancer  == "LIHC") {
    maf <- fread(input = paste0(resultD, "preprocess_files/tables/parse_China_Liver/LIHC_AllSamples.filtered.mutect2.partID.maf"), data.table = F, fill=TRUE) 
  } else if (cancer == "GBM") {
    maf <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC3-GBM/Data_Freeze/cptac3_gbm_data_freeze_v2.1.20190927/somatic_mutation/tindaisy/tindaisy_all_cases_filtered.v2.1.20190927.maf.gz"), data.table = F)
  }
  maf <- data.frame(maf)
  print(paste0("MAF has ", nrow(maf), " lines\n"))
  
  return(maf)
}

# Load cBioPortal data ----------------------------------------------------
cBioPortalD <- "./Ding_Lab/Projects_Current/cBioPortal/"
loadTCGAMaf <- function(cancer) {
  files_tmp <- list.files(path = paste0(cBioPortalD, cancer, "/"), recursive = T)
  files_tmp <- files_tmp[grepl(x = files_tmp, pattern = "data_mutations_extended.txt") & grepl(x = files_tmp, pattern = "tcga", ignore.case = T)]
  files_tmp
  if (length(files_tmp) > 1) {
    stop("multiple maf files!")
  }
  if (cancer != "UCEC") {
    maf <- read_delim(file = paste0(cBioPortalD,  cancer, "/", files_tmp), 
                      "\t", escape_double = FALSE, trim_ws = TRUE, col_names = T)
  } else {
    maf <- read_delim(file = paste0(cBioPortalD,  cancer, "/", files_tmp), 
                      "\t", escape_double = FALSE, trim_ws = TRUE, col_names = T, skip = 1)
  }

  maf <- data.frame(maf)
  print(paste0("MAF has ", nrow(maf), " lines\n"))
  return(maf)
}

loadTCGACNA <- function(cancer) {
  files_tmp <- list.files(path = paste0(cBioPortalD, cancer, "/"), recursive = T)
  files_tmp <- files_tmp[grepl(x = files_tmp, pattern = "data_CNA") & grepl(x = files_tmp, pattern = "tcga", ignore.case = T) & !grepl(x = files_tmp, pattern = "microRNA", ignore.case = T) ]
  if (length(files_tmp) > 1) {
    stop("multiple CNA files!")
  }
  
  data_CNA <- read_delim(file = paste0(cBioPortalD,  cancer, "/", files_tmp), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
  colnames(data_CNA)
  partIDs <- colnames(data_CNA)[!(colnames(data_CNA) %in% c("Hugo_Symbol", "Entrez_Gene_Id"))]
  
  data_CNA_mat <- as.matrix(data_CNA[, colnames(data_CNA)[!(colnames(data_CNA) %in% c("Hugo_Symbol", "Entrez_Gene_Id"))]])
  cna_status_mat <- matrix(data = "neutral", nrow = nrow(data_CNA_mat), ncol = ncol(data_CNA_mat))
  cna_status_mat[which(data_CNA_mat <= -1, arr.ind = T)] <- "deletion"
  cna_status_mat[which(data_CNA_mat >= 1, arr.ind = T)] <- "amplification"
  cna_status_df <- data.frame(gene = data_CNA$Hugo_Symbol)
  cna_status_df <- cbind(cna_status_df, data.frame(cna_status_mat))
  colnames(cna_status_df) <- c("gene", colnames(data_CNA)[!(colnames(data_CNA) %in% c("Hugo_Symbol", "Entrez_Gene_Id"))])
  print("WARNING: CNA file might have multiple entries for one gene!")
  return(cna_status_df)
}

loadTCGARNA <- function(cancer) {
  files_tmp <- list.files(path = paste0(cBioPortalD, cancer, "/"), recursive = T)
  files_tmp
  files_tmp <- files_tmp[grepl(x = files_tmp, pattern = "median_Zscores") & grepl(x = files_tmp, pattern = "mRNA") & !grepl(x = files_tmp, pattern = "meta") & grepl(x = files_tmp, pattern = "tcga", ignore.case = T)]
  files_tmp
  if (length(files_tmp) > 1) {
    files_tmp <- files_tmp[grepl(x = files_tmp, pattern = "v2")]
  }
  if (length(files_tmp) > 1) {
    stop("multiple RNA files!")
  }
  if (length(files_tmp) == 0) {
    files_tmp <- list.files(path = paste0(cBioPortalD, cancer, "/"), recursive = T)
    
    files_tmp <- files_tmp[grepl(x = files_tmp, pattern = "median_Zscores") & grepl(x = files_tmp, pattern = "data_RNA") & grepl(x = files_tmp, pattern = "tcga", ignore.case = T)]
    
  }
  if (cancer == "OV") {
    files_tmp <- list.files(path = paste0(cBioPortalD, cancer, "/"), recursive = T)
    
    files_tmp <- files_tmp[grepl(x = files_tmp, pattern = "median_Zscores") & grepl(x = files_tmp, pattern = "expression") & !grepl(x = files_tmp, pattern = "meta") & grepl(x = files_tmp, pattern = "tcga", ignore.case = T)]
  }
  files_tmp
  data_RNA <- read_delim(file = paste0(cBioPortalD,  cancer, "/", files_tmp), 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
  colnames(data_RNA)
  head(data_RNA)
  partIDs <- colnames(data_RNA)[!(colnames(data_RNA) %in% c("Hugo_Symbol", "Entrez_Gene_Id"))]
  RNA_status_df <- data.frame(gene = data_RNA$Hugo_Symbol)
  RNA_status_df <- cbind(RNA_status_df, data_RNA[, partIDs])
  colnames(RNA_status_df) <- c("gene", partIDs)
  print("WARNING: RNA file might have multiple entries for one gene!")
  head(RNA_status_df)
  return(RNA_status_df)
}

loadTCGARPPA <- function(cancer, expression_type) {
  file_tmp <- paste0(resultD, "preprocess_files/tables/parse_TCGA_data/", cancer, "_", expression_type, "_tumor_TCGA_RPPA_partID.txt")
  tab_tmp <- read_delim(file = file_tmp, 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
  return(tab_tmp)
}

# generate mutation matrix ------------------------------------------------
generate_somatic_mutation_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  maf <- maf[maf$Hugo_Symbol %in% genes4mat,]
  nrow(maf)
  maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    variant_class <- paste0(unique(x), collapse = ",")
    return(variant_class)
  }, value.var = "Variant_Classification", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  return(mut_mat)
}

# functions relating to CNA -----------------------------------------------

loadCNA <- function(cancer) {
  ## input CNA values
  if (cancer %in% c("BRCA", "OV", "CO")) {
    ### if it is breast cancer data, irregular columns names actually won't overlap with proteomics data
    cna <- fread(input = paste0(dir2cptac2_prospective, "copy_number/gatk/v1.3.swap_contamination_fixed/prospective_somatic/gene_level/v1.3.CPTAC2_prospective.2018-03-19/", toupper(substr(cancer, start = 1, stop = 2)), "/gene_level_CNV.", substr(cancer, start = 1, stop = 2), ".v1.3.2018-03-19.tsv"), data.table = F)
  }
  if (cancer %in% c("UCEC", "CCRCC")) {
    cna <- fread(input = paste0(resultD, "preprocess_files/tables/parse_", cancer, "_data_freeze/somatic_CNA.", cancer, ".partID.txt"), data.table = F)
  }
  if (cancer %in% c("LIHC")) {
    cna <- fread(input = paste0(resultD, "preprocess_files/tables/parse_", "China_Liver", "/somatic_CNA.", cancer, ".partID.txt"), data.table = F)
  }
  if (cancer == "GBM") {
    cna <- fread(input = paste0(resultD, "preprocess_files/tables/parse_", cancer, "_data_freeze_v2_0/somatic_CNA.", cancer, ".partID.txt"), data.table = F)
  }
  return(cna)
}

amp_thres_cans <- c(log2(1.1), log2(1.1), log2(1.1), 0.2, 0.1, 0.8, 0.2); names(amp_thres_cans) <- c(cancers_sort, "UCEC", "CCRCC", "LIHC", "GBM")
del_thres_cans <- c(log2(0.9), log2(0.9), log2(0.9), -0.2, -0.1, -0.8, -0.2); names(del_thres_cans) <- c(cancers_sort, "UCEC", "CCRCC", "LIHC", "GBM")

loadCNAstatus <- function(cancer) {
  ## input CNA values
  if (cancer %in% c("BRCA", "OV", "CO")) {
    ### if it is breast cancer data, irregular columns names actually won't overlap with proteomics data
    cna <- fread(input = paste0(dir2cptac2_prospective, "copy_number/gatk/v1.3.swap_contamination_fixed/prospective_somatic/gene_level/v1.3.CPTAC2_prospective.2018-03-19/", toupper(substr(cancer, start = 1, stop = 2)), "/gene_level_CNV.", substr(cancer, start = 1, stop = 2), ".v1.3.2018-03-19.tsv"), data.table = F)
  }
  if (cancer %in% c("UCEC", "CCRCC")) {
    cna <- fread(input = paste0(resultD, "preprocess_files/tables/parse_", cancer, "_data_freeze/somatic_CNA.", cancer, ".partID.txt"), data.table = F)
  }
  if (cancer %in% c("LIHC")) {
    cna <- fread(input = paste0(resultD, "preprocess_files/tables/parse_", "China_Liver", "/somatic_CNA.", cancer, ".partID.txt"), data.table = F)
  }
  if (cancer == "GBM") {
    cna <- fread(input = paste0(resultD, "preprocess_files/tables/parse_", cancer, "_data_freeze_v2_0/somatic_CNA.", cancer, ".partID.txt"), data.table = F)
  }
  cna_head <- cna$gene
  cna_mat <- cna[, colnames(cna)[!(colnames(cna) %in% "gene")]]
  cna_status <- matrix(data = "neutral", nrow = nrow(cna_mat), ncol = ncol(cna_mat))
  cna_status[cna_mat > amp_thres_cans[cancer]] <- "gain"
  cna_status[cna_mat < del_thres_cans[cancer]] <- "loss"
  cna_status <- data.frame(cbind(cna$gene, cna_status))
  colnames(cna_status) <- colnames(cna)
  return(cna_status)
}

loadRNA <- function(cancer) {
  if (cancer == "BRCA") {
    rna_fn <- "fpkm_prospective_breastcancer_110617_swap_corrected.csv"
  }
  if (cancer == "OV") {
    rna_fn <- "fpkm_ov_allgenes.011718_formatted.csv"
  }
  if (cancer == "CO") {
    rna_fn <- "fpkm_crc_allgenes.v1.0.10232017_formatted.csv"
  }
  if (cancer  == "UCEC") {
    rna_fn <- "UCEC_RNAseq_linear_gene_RSEM_UQ_log2(+1)_V2.0_tumor96_PartID.tsv"
  }
  if (cancer  ==  "CCRCC") {
    rna_fn <- "ccRcc_RNA_rpkm_Mich_formatted_tumor.csv"
  }
  if (cancer == "LIHC") {
    rna_fn <- "LIHC_RNA_tumor_PGDAC_RSEM_partID.txt"
  }
  rna_tab <- fread(input = paste0(dir2dinglab_projects, "TP53_shared_data/resources/rna/", rna_fn), data.table = F)
  
  colnames(rna_tab)[1] <- "gene"
  if (cancer %in% c("BRCA", "OV", "CO", "CCRCC", "LIHC")) {
    rna_mat <- as.matrix(rna_tab[,-1])
    rna_mat_log2 <- log2(rna_mat+1)
    rna_tab <- data.frame(gene = rna_tab$gene)
    rna_tab <- cbind(rna_tab, as.data.frame(rna_mat_log2))
  }
  return(rna_tab) 
}

loadParseProteomicsData <- function(cancer, expression_type, sample_type, pipeline_type, norm_type) {
  ## expresson_type: PRO or PHO (phosphosite level) or PHO_collapsed (protein level)
  ## sample_type: tumor or normal
  ## pipeline_type: CDAP or PGDAC
  ## norm_type: unnormalized or scaled
  if (cancer %in% c("BRCA", "OV", "CO")) {
    dir1 <- paste0(resultD, "preprocess_files/tables/parse_cptac2p_3can_data/")
  } else if (pipeline_type == "PGDAC") {
    if (cancer == "LIHC") {
      if (expression_type == "collapsed_PHO") {
        dir1 <- paste0(resultD, "preprocess_files/tables/parse_China_Liver", "/")
      } else {
        dir1 <- paste0(resultD, "preprocess_files/tables/parse_China_Liver_log2_noimputation_NA50_median", "/")
      }
    } else if (cancer == "LUAD") {
      dir1 <- paste0(resultD, "preprocess_files/tables/parse_", cancer, "_data_freeze" , "/")
    } else if (cancer == "UCEC") {
      dir1 <- paste0(resultD, "preprocess_files/tables/parse_", cancer, "_data_freeze", "_v2" , "/")
    } else if (cancer == "GBM") {
      dir1 <- paste0(resultD, "preprocess_files/tables/parse_", cancer, "_data_freeze", "_v2_1" , "/")
    } else {
      dir1 <- paste0(resultD, "preprocess_files/tables/parse_", cancer, "_data_freeze", "" , "/")
    }
  } else if (pipeline_type == "CDAP") {
    dir1 <- paste0(resultD, "preprocess_files/tables/parse_", cancer, "_CDAP_data/")
  }
  exp_data <- fread(input = paste0(dir1, cancer, "_", expression_type, "_", sample_type, "_", pipeline_type, "_", norm_type, "_", "partID", ".txt"), data.table = F)
  return(exp_data)
}

loadGeneList = function(gene_type, cancer, is.soft.limit) {
  if (gene_type == "RTK") {
    RTK_file = read.table(paste(dir2cptac2_retrospective,"reference_files/RTKs_list.txt",sep = ""))
    genelist <- as.vector(t(RTK_file))
  }
  driver_table <- read_excel("~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/Resources/Knowledge/Gene_Lists/mmc1.xlsx", 
                                sheet = "Table S1", skip = 3)
  
  driver_table <- data.frame(driver_table)
  if (cancer != "CO") {
    driver_cancer = driver_table[driver_table$Cancer == cancer, ]
    ## oncogenes for OV: RAB25
  } else if (cancer == "CO") {
    driver_cancer = driver_table[driver_table$Cancer == "COADREAD", ]
  }
  if (gene_type == "driver") {
    genelist <- as.vector(driver_cancer$Gene)
  }
  if (gene_type == "tsg" || gene_type == "oncogene") {
    if (is.soft.limit == "soft") {
      genelist <- as.vector(driver_cancer$Gene)[grepl(pattern = gene_type, x = as.vector(driver_cancer$Tumor.suppressor.or.oncogene.prediction..by.20.20..))]
    } else {
      genelist <- as.vector(driver_cancer$Gene)[driver_cancer$Tumor.suppressor.or.oncogene.prediction..by.20.20.. == gene_type]
    }
  }
  
  return(genelist)
}


# ordering ----------------------------------------------------------------
order_exp_type <- function(x) {
  y <- factor(x, levels = c("RNA", "PRO", "PHO"))
  return(y)
}

order_variant_class <- function(x) {
  y <- factor(x, levels = c("missense", "truncation", "not_silent"))
  return(y)
}


order_cancer <- function(x) {
  y <- factor(x, levels = c("BRCA", "OV", "CO", "UCEC", "CCRCC", "LIHC"))
  return(y)
}

order_cancer_rev <- function(x) {
  y <- factor(x, levels = rev(c("BRCA", "OV", "CO", "UCEC", "CCRCC", "LIHC")))
  return(y)
}

order_pathway <- function(x) {
  pathways_x <- unique(x)
  pathway_order <- c("TP53", "Cell Cycle","WNT", "PI3K","RTK RAS", "MAPK", "SWI/SNF", "NOTCH", "MMR", "HIPPO", "TGF-Beta", "other")
  y <- factor(x, levels = c(pathway_order, pathways_x[!(pathways_x %in% pathway_order)]))
  return(y)
}

# take top ----------------------------------------------------------------
get_top_QT_by_id <- function(value_vector, id_vector, qt) {
  values2keep <- vector(mode = "logical", length = length(id_vector))
  for (id in unique(id_vector)) {
    values2extract <- value_vector[id_vector == id]
    qt_bottom <- quantile(x = values2extract, probs = qt, na.rm = T)
    values2keep[value_vector > qt_bottom & id_vector == id] <- T
  }
  return(values2keep)
}

get_bottom_QT_by_id <- function(value_vector, id_vector, qt) {
  values2keep <- vector(mode = "logical", length = length(id_vector))
  for (id in unique(id_vector)) {
    values2extract <- value_vector[id_vector == id]
    qt_bottom <- quantile(x = values2extract, probs = qt, na.rm = T)
    values2keep[value_vector < qt_bottom & id_vector == id] <- T
  }
  return(values2keep)
}

get_top_by_id <- function(value_vector, id_vector, num_top) {
  values2keep <- vector(mode = "logical", length = length(id_vector))
  for (id in unique(id_vector)) {
    values2extract <- value_vector[id_vector == id]
    values2extract_ordered <- values2extract[order(values2extract, decreasing = T)]
    values2keep[value_vector %in% values2extract_ordered[1:num_top] & id_vector == id] <- T
  }
  return(values2keep)
}

get_order_by_id <- function(value_vector, id_vector, num_top) {
  value_order <- vector(mode = "numeric", length = length(id_vector))
  for (id in unique(id_vector)) {
    values2extract <- value_vector[id_vector == id]
    value_order[id_vector == id] <- order(values2extract)
  }
  return(value_order)
}

get_SMG_by_cancer <- function(gene_vector, cancer_vector) {
  is.smg <- vector("logical", length = length(gene_vector))
  for (cancer in unique(cancer_vector)) {
    is.smg[cancer_vector ==  cancer] <- (gene_vector[cancer_vector ==  cancer] %in% SMGs[[cancer]])
  }
  return(is.smg)
}

get_CNA_genes_by_cancer <- function(gene_vector, cancer_vector) {
  is.smg <- vector("logical", length = length(gene_vector))
  for (cancer in unique(cancer_vector)) {
    is.smg[cancer_vector ==  cancer] <- (gene_vector[cancer_vector ==  cancer] %in% SMGs[[cancer]])
  }
  return(is.smg)
}

# number manipulation -----------------------------------------------------
remove_outliers_IQR <- function(x, out_thres, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- out_thres * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

remove_outliers_cap <- function(x, cap, na.rm = TRUE, ...) {
  y <- x
  y[x < (-cap)] <- NA
  y[x > (cap)] <- NA
  y
}

add_cap <- function(x, cap) {
  y <- x
  y[!is.na(x) & x < (-cap)] <- (-cap)
  y[!is.na(x) & x > (cap)] <- cap
  y
}

range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

FDR_by_id_columns <- function(p_vector, id_columns, df) {
  ## make a replicate the id columns and make it charactors
  df_id <- matrix(data = as.character(unlist(df[, id_columns])), nrow = nrow(df), ncol = length(id_columns), byrow = F)
  df_id <- data.frame(df_id); colnames(df_id) <- id_columns
  df_id_combo <- data.frame(table(df_id))
  df_id_combo <- df_id_combo[df_id_combo$Freq > 0,]
  ## give a number for each combo
  df_id_combo$combo_id <- 1:nrow(df_id_combo)
  df_id_combo
  df_id <- merge(df_id, df_id_combo[, c(id_columns, "combo_id")], all.x = T)
  if (any(is.na(df_id$combo_id))) {
    stop()
  }
  
  ## for every combo of values in id columns, adjust a set of pvalues
  fdr_vector <- vector(mode = "numeric", length = length(p_vector)) + NA
  for (i in 1:length(unique(df_id$combo_id))) {
    row2adjust <- (df_id$combo_id == i & !is.na(p_vector))
    fdr_vector[row2adjust] <- p.adjust(p = p_vector[row2adjust], method = "fdr")
  }
  return(fdr_vector)
}

