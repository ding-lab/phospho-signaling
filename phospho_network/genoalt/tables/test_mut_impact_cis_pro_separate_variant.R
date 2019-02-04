# Yige Wu @ WashU 2018 Aug
## for each kinase-affectedGene pairs, compare the affectedGene phosphorylation level of kinase-mutated tumors, kinase-amplified tumors, kinase-deleted tumors against other tumors
## using wilcox test
## reference: https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
# source("./cptac2p_analysis/phospho_network/genoalt/tables/generate_somatic_mutation_matrix.R")

# set variables -----------------------------------------------------------
num_control_thres <- 4
num_genoalt_thres <- 4
sample_type <- "tumor"
sig_thres <- 0.05
id_vars <- c("alteredGene", "affectedGene", "Phosphosite", "cancer")
prefix_vars <- c("p", "meddiff_bottom", "meddiff_upper", "meddiff", "num")

# inputs ------------------------------------------------------------------
## input alteredGene-affectedGene table
pair_tab <- data.frame(alteredGene = genes4mat, affectedGene = genes4mat)
print(nrow(pair_tab))

# not separate missense from truncations ----------------------------------
mut_cnv_cans <- NULL
for (cancer in c("UCEC", "CCRCC", "BRCA", "OV", "CO")) {
  exp_data <- loadParseProteomicsData(data_type = "PRO", cancer = cancer, sample_type = "tumor")
  partIDs <- colnames(exp_data)
  partIDs <- partIDs[!(partIDs %in% c("Gene", "Phosphosite", "Peptide_ID"))]
  
  ## input somatic mutation matrix
  mut_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_mutation_matrix/", cancer, "_somatic_mutation.txt"), data.table = F)
  
  ## get the overlap of the participant IDs
  partIDs_overlap <- intersect(partIDs, colnames(mut_mat))
  print(length(partIDs_overlap))
  
  ## the number of mutations
  mut_count <- rowSums(x = ((mut_mat[, -1] != "") & (mut_mat[, -1] != "Silent"))); names(mut_count) <- mut_mat$Hugo_Symbol
  
  ## get the genes to test
  genes_num_genoalt_thresholded <- names(mut_count)[mut_count >= num_genoalt_thres]
  print(length(genes_num_genoalt_thresholded))
  
  ## initiate identifier columns
  pair_tab_genoalt_thresholded <- pair_tab[pair_tab$alteredGene %in% genes_num_genoalt_thresholded,]
  print(nrow(pair_tab_genoalt_thresholded))
  
  ## initiate value columns
  num_control <- vector(mode = "numeric", length = nrow(pair_tab)); num_mut <- num_control; 
  meddiff_mut <- num_control; p_mut <- num_control; meddiff_bottom_mut <- num_control; meddiff_upper_mut <- num_control
  
  for (i in 1:nrow(pair_tab)) {
    alteredGene <- pair_tab[i, "alteredGene"]
    affectedGene <- pair_tab[i, "affectedGene"]
    mut_mat_en <- mut_mat[mut_mat$Hugo_Symbol == alteredGene,]
    if (nrow(mut_mat_en) > 0){
      partIDs_mut <- partIDs_overlap[(mut_mat_en[,partIDs_overlap] != "") & (mut_mat_en[, partIDs_overlap] != "Silent")]
    } else {
      partIDs_mut <- NULL 
    }
    
    mut_cnv_partIDs <- partIDs_mut
    control_partIDs <- partIDs_overlap[!(partIDs_overlap %in% mut_cnv_partIDs)]
    
    print(paste0("alteredGene:", alteredGene))
    print(paste0("affectedGene:", affectedGene))
    
    control_pro <- exp_data[exp_data$Gene == alteredGene, control_partIDs]
    control_pro <- control_pro[!is.na(control_pro)]
    num_control[i] <- length(control_pro)
    
    if (num_control[i] >= num_control_thres) {
      mut_pro <- exp_data[exp_data$Gene == alteredGene, partIDs_mut]
      mut_pro <- mut_pro[!is.na(mut_pro)]
      num_mut[i] <- length(mut_pro)
      if (num_mut[i] > 0) {
        meddiff_mut[i] <- median(mut_pro) - median(control_pro)
        if (num_mut[i] >= num_genoalt_thres){
          stat <- wilcox.test(x = mut_pro, y = control_pro, conf.int = T)
          p_mut[i] <- stat$p.value
          meddiff_bottom_mut[i] <- stat$conf.int[1]
          meddiff_upper_mut[i] <- stat$conf.int[2]
        }
      }
    }
    
  }

  mut_cnv_tab <- pair_tab
  mut_cnv_tab$Phosphosite  <- "protein"
  mut_cnv_tab <- cbind(mut_cnv_tab, data.frame(num_mut, num_control,
                                               p_mut, meddiff_bottom_mut, meddiff_upper_mut, meddiff_mut))
  mut_cnv_tab$cancer <- cancer
  mut_cnv_tab$pair_pro <- paste0(mut_cnv_tab$alteredGene, ":", mut_cnv_tab$affectedGene)
  mut_cnv_tab$pair <- paste0(mut_cnv_tab$alteredGene, ":", mut_cnv_tab$affectedGene, ":", mut_cnv_tab$Phosphosite)

  ## annotate cis and trans
  mut_cnv_tab$SELF <- ifelse(test = as.vector(mut_cnv_tab$alteredGene) == as.vector(mut_cnv_tab$affectedGene), yes = "cis", no = "trans")
  
  ## filter mutation-impacted kinase-affectedGene pairs
  mut_tab <- mut_cnv_tab
  mut_tab <- mut_tab[mut_tab$p_mut > 0,]
  
  mut_tab.f <- mut_tab[, c(id_vars, paste0(paste0(prefix_vars, "", sep = "_"), "mut"), "pair")]
  colnames(mut_tab.f) <- c(id_vars, prefix_vars, "pair")
  mut_tab.f$genoalt_type <- "mut"

  mut_cnv_can <- mut_tab.f
  mut_cnv_cans <- rbind(mut_cnv_cans, mut_cnv_can)
  mut_cnv_cans <- unique(mut_cnv_cans)
  
  write.table(x = mut_cnv_tab, file = paste0(makeOutDir(resultD = resultD), cancer, "_mut_cnv_tab.txt"), row.names = F, quote = F)
  write.table(x = mut_cnv_can, file = paste0(makeOutDir(resultD = resultD), cancer, "_mut_cnv_tab_formatted.txt"), row.names = F, quote = F, sep = "\t")
}
mut_cnv_cans$SELF <- ifelse(as.vector(mut_cnv_cans$alteredGene) == as.vector(mut_cnv_cans$affectedGene), "cis", "trans")
mut_cnv_cans$driver_gene_type.en <- ""
mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$alteredGene %in% oncogenes] <- "oncogene"
mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$alteredGene %in% tsgs] <- "tsg"
write.table(x = mut_cnv_cans, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_sig_cans.txt"), row.names = F, quote = F, sep = "\t")

# separate missense from truncations ----------------------------------
mut_cnv_cans <- NULL
for (cancer in c("UCEC", "CCRCC", "BRCA", "OV", "CO")) {
  exp_data <- loadParseProteomicsData(data_type = "PRO", cancer = cancer, sample_type = "tumor")
  partIDs <- colnames(exp_data)
  partIDs <- partIDs[!(partIDs %in% c("Gene", "Phosphosite", "Peptide_ID"))]
  
  ## input somatic mutation matrix
  mut_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_mutation_matrix/", cancer, "_somatic_mutation.txt"), data.table = F)
  
  ## get the overlap of the participant IDs
  partIDs_overlap <- intersect(partIDs, colnames(mut_mat))
  print(length(partIDs_overlap))
  
  ## the number of mutations
  mut_count <- rowSums(x = ((mut_mat[, -1] != "") & (mut_mat[, -1] != "Silent"))); names(mut_count) <- mut_mat$Hugo_Symbol
  
  ## get the genes to test
  genes_num_genoalt_thresholded <- names(mut_count)[mut_count >= num_genoalt_thres]
  print(length(genes_num_genoalt_thresholded))
  
  ## initiate identifier columns
  pair_tab_genoalt_thresholded <- pair_tab[pair_tab$alteredGene %in% genes_num_genoalt_thresholded,]
  print(nrow(pair_tab_genoalt_thresholded))
  
  ## initiate value columns
  num_control <- vector(mode = "numeric", length = nrow(pair_tab)); 
  num_mut <- num_control; meddiff_mut <- num_control; p_mut <- num_control; meddiff_bottom_mut <- num_control; meddiff_upper_mut <- num_control
  num_missense <- num_control; meddiff_missense <- num_control; p_missense <- num_control; meddiff_bottom_missense <- num_control; meddiff_upper_missense <- num_control
  num_truncation <- num_control; meddiff_truncation <- num_control; p_truncation <- num_control; meddiff_bottom_truncation <- num_control; meddiff_upper_truncation <- num_control
  
  for (i in 1:nrow(pair_tab)) {
    alteredGene <- pair_tab[i, "alteredGene"]
    affectedGene <- pair_tab[i, "affectedGene"]
    mut_mat_en <- mut_mat[mut_mat$Hugo_Symbol == alteredGene,]
    partIDs_mut <- partIDs_overlap[(mut_mat_en[,partIDs_overlap] != "") & (mut_mat_en[, partIDs_overlap] != "Silent")]
    partIDs_missense <- partIDs_overlap[grep(pattern = "Missense_Mutation", x = mut_mat_en[,partIDs_overlap])]
    partIDs_truncation <- partIDs_overlap[grep(pattern = "Frame_Shift_Del|Frame_Shift_Ins|Nonsense_Mutation|Nonstop_Mutation|Splice_Site", x = mut_mat_en[,partIDs_overlap])]
    
    mut_cnv_partIDs <- partIDs_mut
    control_partIDs <- partIDs_overlap[!(partIDs_overlap %in% mut_cnv_partIDs)]
    
    control_pro <- exp_data[exp_data$Gene == alteredGene, control_partIDs]
    control_pro <- control_pro[!is.na(control_pro)]
    num_control[i] <- length(control_pro)
    
    if (num_control[i] >= num_control_thres) {
      print(paste0("alteredGene:", alteredGene))
      print(paste0("affectedGene:", affectedGene))
      
      exp_missense <- exp_data[exp_data$Gene == alteredGene, partIDs_missense]
      exp_missense <- exp_missense[!is.na(exp_missense)]
      num_missense[i] <- length(exp_missense)
      if (num_missense[i] > 0) {
        meddiff_missense[i] <- median(exp_missense) - median(control_pro)
        if (num_missense[i] >= num_genoalt_thres){
          stat <- wilcox.test(x = exp_missense, y = control_pro, conf.int = T)
          p_missense[i] <- stat$p.value
          meddiff_bottom_missense[i] <- stat$conf.int[1]
          meddiff_upper_missense[i] <- stat$conf.int[2]
        }
      }
      
      exp_truncation <- exp_data[exp_data$Gene == alteredGene, partIDs_truncation]
      exp_truncation <- exp_truncation[!is.na(exp_truncation)]
      num_truncation[i] <- length(exp_truncation)
      if (num_truncation[i] > 0) {
        meddiff_truncation[i] <- median(exp_truncation) - median(control_pro)
        if (num_truncation[i] >= num_genoalt_thres){
          stat <- wilcox.test(x = exp_truncation, y = control_pro, conf.int = T)
          p_truncation[i] <- stat$p.value
          meddiff_bottom_truncation[i] <- stat$conf.int[1]
          meddiff_upper_truncation[i] <- stat$conf.int[2]
        }
      }
    }
  }
  
  mut_cnv_tab <- pair_tab
  mut_cnv_tab$Phosphosite  <- "protein"
  mut_cnv_tab <- cbind(mut_cnv_tab, data.frame(num_mut, num_control,
                                               p_mut, meddiff_bottom_mut, meddiff_upper_mut, meddiff_mut))
  mut_cnv_tab$cancer <- cancer
  mut_cnv_tab$pair_pro <- paste0(mut_cnv_tab$alteredGene, ":", mut_cnv_tab$affectedGene)
  mut_cnv_tab$pair <- paste0(mut_cnv_tab$alteredGene, ":", mut_cnv_tab$affectedGene, ":", mut_cnv_tab$Phosphosite)
  
  ## annotate cis and trans
  mut_cnv_tab$SELF <- ifelse(test = as.vector(mut_cnv_tab$alteredGene) == as.vector(mut_cnv_tab$affectedGene), yes = "cis", no = "trans")
  
  ## filter mutation-impacted kinase-affectedGene pairs
  mut_tab <- mut_cnv_tab
  mut_tab <- mut_tab[mut_tab$p_mut > 0,]
  
  mut_tab.f <- mut_tab[, c(id_vars, paste0(paste0(prefix_vars, "", sep = "_"), "mut"), "pair")]
  colnames(mut_tab.f) <- c(id_vars, prefix_vars, "pair")
  mut_tab.f$genoalt_type <- "mut"
  
  mut_cnv_can <- mut_tab.f
  mut_cnv_cans <- rbind(mut_cnv_cans, mut_cnv_can)
  mut_cnv_cans <- unique(mut_cnv_cans)
  
  write.table(x = mut_cnv_tab, file = paste0(makeOutDir(resultD = resultD), cancer, "_mut_cnv_tab.txt"), row.names = F, quote = F)
  write.table(x = mut_cnv_can, file = paste0(makeOutDir(resultD = resultD), cancer, "_mut_cnv_tab_formatted.txt"), row.names = F, quote = F, sep = "\t")
}
mut_cnv_cans$SELF <- ifelse(as.vector(mut_cnv_cans$alteredGene) == as.vector(mut_cnv_cans$affectedGene), "cis", "trans")
mut_cnv_cans$driver_gene_type.en <- ""
mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$alteredGene %in% oncogenes] <- "oncogene"
mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$alteredGene %in% tsgs] <- "tsg"
write.table(x = mut_cnv_cans, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_sig_cans.txt"), row.names = F, quote = F, sep = "\t")



stop("")
mut_cnv_cans <- NULL
for (cancer in c("UCEC")) {
  mut_cnv_can <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/genoalt/tables/test_mut_cna_impact_cis_pro_cptac3/", cancer, "_mut_cnv_tab_formatted.txt"), data.table = F)
  mut_cnv_cans <- rbind(mut_cnv_cans, mut_cnv_can)
}
mut_cnv_cans$SELF <- ifelse(as.vector(mut_cnv_cans$GENE) == as.vector(mut_cnv_cans$SUB_GENE), "cis", "trans")
mut_cnv_cans$driver_gene_type.en <- ""
mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$GENE %in% oncogenes] <- "oncogene"
mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$GENE %in% tsgs] <- "tsg"
write.table(x = mut_cnv_cans, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_cans.txt"), row.names = F, quote = F, sep = "\t")

