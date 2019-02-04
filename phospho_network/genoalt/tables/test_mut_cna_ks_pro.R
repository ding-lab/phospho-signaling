# Yige Wu @ WashU 2018 Aug
## for each complex containing protein A & B, compare the protein and phosphorylation of B between A-mutated and A-wildtype samples
## using wilcox test
## reference: https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/

# source ------------------------------------------------------------------
source("Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("Box Sync/cptac2p_analysis/phospho_network/genoalt/tables/generate_somatic_mutation_matrix.R")
source("Box Sync/cptac2p_analysis/phospho_network/genoalt/tables/generate_somatic_cna_matrix.R")

# set variables -----------------------------------------------------------
num_control_thres <- 4
num_genoalt_thres <- 4
sample_type <- "tumor"
sig_thres <- 0.05
prefix_vars <- c("p", "num", "meddiff_bottom", "meddiff_upper", "meddiff")
id_vars <- c("geneA", "geneB", "SUB_MOD_RSD", "num_control", "cancer")

# inputs ------------------------------------------------------------------
## input complex protein pair table
complex_pair_tab <- fread(input = "./cptac2p/analysis_results/phospho_network/compile_enzyme_substrate/tables/parse_corum_signor_reactome/sup_complex_pair_uniq.txt", data.table = F)
print(nrow(complex_pair_tab))
pair_tab <- complex_pair_tab


# Business  ---------------------------------------------------------------
testMutCnaComplexPro <- function(cancer) {
  affected_exp_type <- "protein"
  
  if (cancer %in% cancers_sort) {
    pro_data <- loadProteinNormalizedTumor(cancer = cancer)
    pho_data <- loadPhosphositeNormalizedTumor(cancer = cancer)
    sampIDs <- colnames(pho_data)[!(colnames(pho_data) %in% c("Gene", "Phosphosite"))]
    pho_data <- pho_data[rowSums(!is.na(pho_data[, sampIDs])) >= (num_control_thres + num_genoalt_thres),]
    site_head <- formatPhosphosite(pho_data$Phosphosite, pho_data$Gene)
    partIDs <- sampID2partID(sampleID_vector = sampIDs, sample_map = loadSampMap())
    names(pho_data) <- c("Gene", "Phosphosite", partIDs)
    names(pro_data) <- c("Gene", partIDs)
  }
  if (cancer == "UCEC") {
    ## input phosphorylation data
    pro_data <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1/UCEC_proteomics_V1.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
    pro_data <- data.frame(pro_data)
    colnames(pro_data)[1] <- "Gene"
    
    pho_data <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1/UCEC_phosphoproteomics_site_level_V1.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
    pho_data <- data.frame(pho_data)
    site_head <- data.frame(str_split_fixed(string = pho_data$idx, pattern = "-", n = 2))
    colnames(site_head) <- c("affected_gene", "SUB_MOD_RSD")
    
    samples <- colnames(pro_data)[!(colnames(pro_data) %in% c("idx")) & !grepl(pattern = "Mix", x = (colnames(pro_data)))  & !grepl(pattern = "rep", x = (colnames(pro_data)))]
    samples <- samples[grepl(pattern = paste0('\\.', toupper(substr(x = sample_type, start = 1, stop = 1))), x = samples)]
    sampIDs <- samples
    tmp <- str_split_fixed(string = samples, pattern = "\\.", n = 3)
    partIDs <- paste0(tmp[,1], "-", tmp[,2])
    
    pho_data <- pho_data[, sampIDs]
    colnames(pho_data) <- partIDs
    
    pro_data <- pro_data[, c("Gene", sampIDs)]
    colnames(pro_data) <- c("Gene", partIDs)
  }
  if (affected_exp_type == "protein") {
    affected_exp_data <- pro_data
    site_head <- data.frame(geneB = affected_exp_data$Gene, SUB_MOD_RSD = "protein")
  } else {
    stop("what kind of expression is affected?")
  }
  
  ## input somatic mutation matrix
  mut_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_mutation_matrix/", cancer, "_somatic_mutation.txt"), data.table = F)
  
  ## input CNA matrix
  shallow_del_thresholded_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_cna_matrix/", cancer, "_shallow_deletions_", del_thres_cans[cancer], ".txt"), data.table = F)
  shallow_amp_thresholded_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_cna_matrix/", cancer,"_shallow_amplifications_", amp_thres_cans[cancer], ".txt"), data.table = F)
  
  ## get the overlap of the participant IDs
  partIDs_overlap <- intersect(colnames(shallow_del_thresholded_mat), intersect(partIDs, colnames(mut_mat)))
  print(length(partIDs_overlap))
  
  ## the number of mutations
  mut_count <- rowSums(x = ((mut_mat[, -1] != "") & (mut_mat[, -1] != "Silent"))); names(mut_count) <- mut_mat$Hugo_Symbol
  
  ## the number of CNA
  shallow_del_count <- rowSums(shallow_del_thresholded_mat[, partIDs_overlap]); names(shallow_del_count) <- shallow_del_thresholded_mat$Hugo_Symbol
  shallow_amp_count <- rowSums(shallow_amp_thresholded_mat[, partIDs_overlap]); names(shallow_amp_count) <- shallow_amp_thresholded_mat$Hugo_Symbol
  
  ## get the genes to test
  genes_num_genoalt_thresholded <- union(names(mut_count)[mut_count >= num_genoalt_thres], names(shallow_del_count)[shallow_del_count >= num_genoalt_thres])
  genes_num_genoalt_thresholded <- union(genes_num_genoalt_thresholded, names(shallow_amp_count)[shallow_amp_count >= num_genoalt_thres])
  print(length(genes_num_genoalt_thresholded))
  
  ## initiate identifier columns
  pair_tab_genoalt_thresholded <- pair_tab[pair_tab$geneA %in% genes_num_genoalt_thresholded,]
  print(nrow(pair_tab_genoalt_thresholded))
  
  pair_tab_sites <- merge(pair_tab_genoalt_thresholded, site_head[, c("geneB", "SUB_MOD_RSD")], by.x = c("geneB"), by.y = c("geneB"), all.y = T)
  pair_tab_sites <- pair_tab_sites[!is.na(pair_tab_sites$geneA) & !is.na(pair_tab_sites$SUB_MOD_RSD),]
  pair_tab_sites <- rbind(pair_tab_sites[pair_tab_sites$geneA == "TP53",], pair_tab_sites[pair_tab_sites$geneA != "TP53",])
  print(nrow(pair_tab_sites))
  
  ## initiate value columns
  num_control <- vector(mode = "numeric", length = nrow(pair_tab_sites))
  num_mut <- num_control; meddiff_mut <- num_control; p_mut <- num_control; meddiff_bottom_mut <- num_control; meddiff_upper_mut <- num_control
  num_shallow_amp <- num_control; meddiff_shallow_amp <- num_control; p_shallow_amp <- num_control; meddiff_bottom_shallow_amp <- num_control; meddiff_upper_shallow_amp <- num_control
  num_shallow_del <- num_control; meddiff_shallow_del <- num_control; p_shallow_del <- num_control; meddiff_bottom_shallow_del <- num_control; meddiff_upper_shallow_del <- num_control
  
  for (altered_gene in unique(pair_tab_sites$geneA)) {
    mut_altered_gene <- mut_mat[mut_mat$Hugo_Symbol == altered_gene,]
    if (nrow(mut_altered_gene) > 0){
      mut_partIDs <- partIDs_overlap[(mut_altered_gene[,partIDs_overlap] != "") & (mut_altered_gene[, partIDs_overlap] != "Silent")]
    } else {
      mut_partIDs <- NULL 
    }
    shallow_amp_thresholded_mat_en <- shallow_amp_thresholded_mat[shallow_amp_thresholded_mat$Hugo_Symbol == altered_gene,]
    if (nrow(shallow_amp_thresholded_mat_en) > 0){
      shallow_amp_partIDs <- partIDs_overlap[(shallow_amp_thresholded_mat_en[,partIDs_overlap] == "TRUE")]
    } else {
      shallow_amp_partIDs <- NULL 
    }
    shallow_del_thresholded_mat_en <- shallow_del_thresholded_mat[shallow_del_thresholded_mat$Hugo_Symbol == altered_gene,]
    if (nrow(shallow_del_thresholded_mat_en) > 0){
      shallow_del_partIDs <- partIDs_overlap[(shallow_del_thresholded_mat_en[,partIDs_overlap] == "TRUE")]
    } else {
      shallow_del_partIDs <- NULL 
    }
    
    mut_cnv_partIDs <- c(mut_partIDs, shallow_amp_partIDs, shallow_amp_partIDs)
    control_partIDs <- partIDs_overlap[!(partIDs_overlap %in% mut_cnv_partIDs)]
    
    print(paste0("altered_gene:", altered_gene))
    for (affected_gene in unique(pair_tab_sites$geneB[pair_tab_sites$geneA == altered_gene])) {
      print(paste0("affected_gene:", affected_gene))
      
      for (site in unique(pair_tab_sites$SUB_MOD_RSD[pair_tab_sites$geneA == altered_gene & pair_tab_sites$geneB == affected_gene])) {
        print(paste0("site:", site))
        
        i <- (pair_tab_sites$geneA == altered_gene & pair_tab_sites$geneB == affected_gene & pair_tab_sites$SUB_MOD_RSD == site)
        exp_affected <- affected_exp_data[site_head$geneB == affected_gene & site_head$SUB_MOD_RSD == site,]
        exp_affected <- exp_affected[1,]
        exp_control <- exp_affected[, control_partIDs]
        exp_control <- exp_control[!is.na(exp_control)]
        num_control[i] <- length(exp_control)
        
        if (num_control[i] >= num_control_thres) {
          exp_mut <- exp_affected[, mut_partIDs]
          exp_mut <- exp_mut[!is.na(exp_mut)]
          num_mut[i] <- length(exp_mut)
          if (num_mut[i] > 0) {
            meddiff_mut[i] <- median(exp_mut) - median(exp_control)
            if (num_mut[i] >= num_genoalt_thres){
              stat <- wilcox.test(x = exp_mut, y = exp_control, conf.int = T)
              p_mut[i] <- stat$p.value
              meddiff_bottom_mut[i] <- stat$conf.int[1]
              meddiff_upper_mut[i] <- stat$conf.int[2]
            }
          }
          
          exp_shallow_amp <- exp_affected[, shallow_amp_partIDs]
          exp_shallow_amp <- exp_shallow_amp[!is.na(exp_shallow_amp)]
          num_shallow_amp[i] <- length(exp_shallow_amp)
          if (num_shallow_amp[i] > 0) {
            meddiff_shallow_amp[i] <- median(exp_shallow_amp) - median(exp_control)
            if (num_shallow_amp[i] >= num_genoalt_thres){
              stat <- wilcox.test(x = exp_shallow_amp, y = exp_control, conf.int = T)
              p_shallow_amp[i] <- stat$p.value
              meddiff_bottom_shallow_amp[i] <- stat$conf.int[1]
              meddiff_upper_shallow_amp[i] <- stat$conf.int[2]
            }
          }
          
          exp_shallow_del <- exp_affected[, shallow_del_partIDs]
          exp_shallow_del <- exp_shallow_del[!is.na(exp_shallow_del)]
          num_shallow_del[i] <- length(exp_shallow_del)
          if (num_shallow_del[i] > 0) {
            meddiff_shallow_del[i] <- median(exp_shallow_del) - median(exp_control)
            if (num_shallow_del[i] >= num_genoalt_thres){
              stat <- wilcox.test(x = exp_shallow_del, y = exp_control, conf.int = T)
              p_shallow_del[i] <- stat$p.value
              meddiff_bottom_shallow_del[i] <- stat$conf.int[1]
              meddiff_upper_shallow_del[i] <- stat$conf.int[2]
            }
          }
        }
      }
    }
    print(which(i))
  }
  mut_cnv_tab <- pair_tab_sites
  mut_cnv_tab <- cbind(mut_cnv_tab, data.frame(num_mut, num_shallow_amp, num_shallow_del, num_control,
                                               p_mut, meddiff_bottom_mut, meddiff_upper_mut, meddiff_mut,
                                               p_shallow_amp, meddiff_bottom_shallow_amp, meddiff_upper_shallow_amp, meddiff_shallow_amp,
                                               p_shallow_del, meddiff_bottom_shallow_del, meddiff_upper_shallow_del, meddiff_shallow_del))
  mut_cnv_tab$cancer <- cancer
  write.table(x = mut_cnv_tab, file = paste0(makeOutDir(resultD = resultD), cancer, "_mut_cnv_tab.txt"), col.names = T, row.names = F, quote = F)
}
lapply(c("BRCA", "OV", "CO", "UCEC"), FUN = testMutCnaComplexPro)

# parse all cancers -------------------------------------------------------
mut_cnv_cans <- NULL
for (cancer in c("BRCA", "OV", "CO", "UCEC")) {
  mut_cnv_tab <- fread(file = paste0(makeOutDir(resultD = resultD), cancer, "_mut_cnv_tab.txt"), data.table = F)
  ## filter mutation-impacted kinase-affected_gene pairs
  mut_tab <- mut_cnv_tab
  mut_tab <- mut_tab[mut_tab$p_mut > 0,]
  
  ## format the mutation table
  mut_tab.f <- mut_tab[, c(id_vars, paste0(paste0(prefix_vars, "", sep = "_"), "mut"))]
  colnames(mut_tab.f) <- c(id_vars, prefix_vars)
  mut_tab.f$genoalt_type <- "mut"
  
  mut_cnv_can <- mut_tab.f
  mut_cnv_cans <- rbind(mut_cnv_cans, mut_cnv_can)
  mut_cnv_cans <- unique(mut_cnv_cans)
}
mut_cnv_cans$driver_gene_type.en <- ""
mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$geneA %in% oncogenes] <- "oncogene"
mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$geneA %in% tsgs] <- "tsg"
write.table(x = mut_cnv_cans, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_sig_cans.txt"), row.names = F, quote = F, sep = "\t")

