## shared objects by TP53 project


# gene list ---------------------------------------------------------------
TP53_GOF_interacting_TFs <- c("SP1", "ETS1", "ETS2", "YAP1", "NFYA", "CREBBP", "NFKB1", "VDR", "SREBP1", "SREBP2", "E2F1", "TP63", "TP73", "SMAD2", "SMAD3", "MRE11", "NFE2L2", "STAT3")

# load data ---------------------------------------------------------------
loadTP53Deletion <- function(cancer) {
  del_thres_TP53 <- c(-0.2, -0.05, -0.2, -0.1, -0.15, -0.4)
  names(del_thres_TP53) <- c("BRCA", "OV", "CO", "UCEC", "CCRCC", "LIHC")
  ## input CNA values
  if (cancer %in% c("BRCA", "OV", "CO")) {
    ### if it is breast cancer data, irregular columns names actually won't overlap with proteomics data
    cna <- fread(input = paste0(cptac2p_genomicD, "copy_number/gatk/v1.3.swap_contamination_fixed/prospective_somatic/gene_level/v1.3.CPTAC2_prospective.2018-03-19/", toupper(substr(cancer, start = 1, stop = 2)), "/gene_level_CNV.", substr(cancer, start = 1, stop = 2), ".v1.3.2018-03-19.tsv"), data.table = F)
  }
  if (cancer %in% c("UCEC", "CCRCC")) {
    cna <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/parse_", cancer, "_data_freeze/somatic_CNA.", cancer, ".partID.txt"), data.table = F)
  }
  if (cancer %in% c("LIHC")) {
    cna <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/parse_", "China_Liver", "/somatic_CNA.", cancer, ".partID.txt"), data.table = F)
  }
  cna <- cna[cna$gene == "TP53",]
  cna_head <- cna$gene
  cna_mat <- cna[, colnames(cna)[!(colnames(cna) %in% "gene")]]
  cna_status <- matrix(data = "neutral", nrow = nrow(cna_mat), ncol = ncol(cna_mat))
  cna_status[cna_mat < del_thres_TP53[cancer]] <- "deletion"
  cna_status <- data.frame(cbind(cna$gene, cna_status))
  colnames(cna_status) <- colnames(cna)
  return(cna_status)
}

loadTP53DeepDeletion <- function(cancer) {
  del_thres_TP53 <- c(-0.6, -0.75, -0.8, -0.95)
  names(del_thres_TP53) <- c("BRCA", "OV", "CO", "LIHC")

  ## input CNA values
  if (cancer %in% c("BRCA", "OV", "CO", "LIHC")) {
    ### if it is breast cancer data, irregular columns names actually won't overlap with proteomics data
    cna <- loadCNA(cancer)
  } else {
    return("No deep deletion in UCEC & CCRCC TCGA cohort")
  }
  cna <- cna[cna$gene == "TP53",]
  cna_head <- cna$gene
  cna_mat <- cna[, colnames(cna)[!(colnames(cna) %in% "gene")]]
  cna_status <- matrix(data = "neutral", nrow = nrow(cna_mat), ncol = ncol(cna_mat))
  cna_status[cna_mat < del_thres_TP53[cancer]] <- "deletion"
  cna_status <- data.frame(cbind(cna$gene, cna_status))
  colnames(cna_status) <- colnames(cna)
  return(cna_status)
}



# sorting -----------------------------------------------------------------
stat_comparisons <- list(c("Deletion_Truncation", "WT"),
                         c("Deletion_Missense", "WT"),
                         c("Truncation", "WT"), 
                         c("Missense", "WT"),
                         c("Deletion", "WT"))
names(stat_comparisons) <- c("Deletion_Truncation", "Deletion_Missense", "Truncation", "Missense", "Deletion")
sort_mutation_class_complex <- function(x) {
  y <- factor(x, levels = c("Deletion_Truncation", "Deletion_Missense", "Truncation", "Missense", "Missense_Truncation", "Deletion", "Other", "WT"))
  return(y)
}



# TP53 domain -------------------------------------------------------------
getTP53domain <- function(position) {
  if (is.na(position)) {
    return("other")
  } else if (position <= 43) {
    return("TAD1")
  } else if (position <= 63) {
    return("TAD2")
  } else if (position <= 92) {
    return("proline rich\ndomain")
  } else  if (position >= 102 & position <= 292) {
    return("DBD")
  } else if (position >= 316 & position <= 325) {
    return("NLS")
  } else if (position >= 325 & position <= 355) {
    return("tetramerization\ndomain")
  } else if (position >= 356 & position <= 393) {
    return("regulatory\ndomain")
  } else {
    return("other")
  }
}

getTP53domains <- function(positions) {
  domains <- sapply(positions, getTP53domain)
  return(domains)
}


# loading data ------------------------------------------------------------
loadTP53proteomics <- function(cancer, expression_type) {
  if (cancer %in% c("BRCA", "OV", "CO")) {
    exp_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = expression_type, sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
  } else if (cancer == "UCEC") {
    exp_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = expression_type, sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
  } else if (cancer == "CCRCC") {
    exp_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = expression_type, sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
  } else if (cancer == "LIHC") {
    exp_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = expression_type, sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
  }
  return(exp_tab)
}

# plotting variables ------------------------------------------------------
colors_mut <- c("#e41a1c", "#984ea3", "#fb9a99", "#cab2d6", "#377eb8", "#666666")
names(colors_mut) <- c("Deletion_Truncation", "Deletion_Missense", "Truncation", "Missense", "Deletion", "WT")
colors_mut_v2 <- c("#e41a1c", "#f46d43", "#fee090", "#cab2d6", "#377eb8", "#666666")
names(colors_mut_v2) <- c("Missense_protein_level_high", "Missense_protein_level_average", "Missense_protein_level_NA", "Truncation", "Deletion", "WT")
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

# plotting ----------------------------------------------------------------
library(RColorBrewer)


if(!("pheatmap" %in% installed.packages()[,"Package"])) {
  install.packages("pheatmap")
}
library(pheatmap)
save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf <- function(x, filename, width=6, height=6) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}