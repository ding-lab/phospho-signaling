# Yige Wu @ WashU 2018 Aug
## for each complex containing protein A & B, compare the protein and phosphorylation of B between A-mutated and A-wildtype samples
## using wilcox test
## reference: https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/

# source ------------------------------------------------------------------
source("Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
num_control_thres <- 4
num_genoalt_thres <- 4
sample_type <- "tumor"
sig_thres <- 0.05
prefix_vars <- c("p", "num", "meddiff_bottom", "meddiff_upper", "meddiff")
id_vars <- c("geneA", "geneB", "SUB_MOD_RSD", "pair_pro")
# inputs ------------------------------------------------------------------
## input complex protein pair table
complex_pair_tab <- fread(input = "./cptac2p/analysis_results/phospho_network/compile_enzyme_substrate/tables/parse_corum_signor_reactome/sup_complex_pair_uniq.txt", data.table = F)
print(nrow(complex_pair_tab))

for (cancer in c("BRCA")) {
  if (cancer %in% cancers_sort) {
    pro_data <- loadProteinNormalizedTumor(cancer = cancer)
    pho_data <- loadPhosphositeNormalizedTumor(cancer = cancer)
    sampIDs <- colnames(pho_data)[!(colnames(pho_data) %in% c("Gene", "Phosphosite"))]
    pho_data <- pho_data[rowSums(!is.na(pho_data[, sampIDs])) >= (num_control_thres + num_genoalt_thres),]
    pho_head <- formatPhosphosite(pho_data$Phosphosite, pho_data$Gene)
    partIDs <- sampID2partID(sampleID_vector = sampIDs, sample_map = loadSampMap())
    names(pho_data) <- c("Gene", "Phosphosite", partIDs)
  }
  if (cancer == "UCEC") {
    ## input phosphorylation data
    pro_data <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1/UCEC_proteomics_V1.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
    pro_data <- data.frame(pro_data)
    colnames(pro_data)[1] <- "Gene"
    
    # pro_data <- loadProteinNormalizedTumor(cancer = cancer)
    pho_data <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1/UCEC_phosphoproteomics_site_level_V1.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
    pho_data <- data.frame(pho_data)
    # pho_data <- loadPhosphositeNormalizedTumor(cancer = cancer)
    pho_gdata <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1/UCEC_phosphoproteomics_gene_level_V1.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
    pho_gdata <- data.frame(pho_gdata)
    colnames(pho_gdata)[1] <- "Gene"
    # pho_gdata <- loadPhosphoproteinNormalizedTumor(cancer = cancer)
    pho_head <- data.frame(str_split_fixed(string = pho_data$idx, pattern = "-", n = 2))
    colnames(pho_head) <- c("SUBSTRATE", "SUB_MOD_RSD")
    
    samples <- colnames(pro_data)[!(colnames(pro_data) %in% c("idx")) & !grepl(pattern = "Mix", x = (colnames(pro_data)))  & !grepl(pattern = "rep", x = (colnames(pro_data)))]
    samples <- samples[grepl(pattern = paste0('\\.', toupper(substr(x = sample_type, start = 1, stop = 1))), x = samples)]
    sampIDs <- samples
    tmp <- str_split_fixed(string = samples, pattern = "\\.", n = 3)
    partIDs <- paste0(tmp[,1], "-", tmp[,2])
    print(length(partIDs))
    pho_data <- pho_data[, sampIDs]
    colnames(pho_data) <- partIDs
  }

  ## input somatic mutation matrix
  mut_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_mutation_matrix/", cancer, "_somatic_mutation_in_enzyme_substrate.txt"), data.table = F)
  
  ## get the overlap of the participant IDs
  partIDs_overlap <- intersect(partIDs, colnames(mut_mat))
  print(length(partIDs_overlap))
  
  ## the number of mutations
  mut_count <- rowSums(x = ((mut_mat[, -1] != "") & (mut_mat[, -1] != "Silent"))); names(mut_count) <- mut_mat$Hugo_Symbol
  
  ## get the genes to test
  genes_num_genoalt_thresholded <- names(mut_count)[mut_count >= num_genoalt_thres]
  print(length(genes_num_genoalt_thresholded))
  
  ## initiate identifier columns
  complex_pair_tab_genoalt_thresholded <- complex_pair_tab[complex_pair_tab$geneA %in% genes_num_genoalt_thresholded,]
  print(nrow(complex_pair_tab_genoalt_thresholded))
  
  complex_pair_tab_sites <- merge(complex_pair_tab_genoalt_thresholded, pho_head[, c("SUBSTRATE", "SUB_MOD_RSD")], by.x = c("geneB"), by.y = c("SUBSTRATE"), all.y = T)
  complex_pair_tab_sites <- complex_pair_tab_sites[!is.na(complex_pair_tab_sites$geneA) & !is.na(complex_pair_tab_sites$SUB_MOD_RSD),]
  complex_pair_tab_sites <- complex_pair_tab_sites[order(complex_pair_tab_sites$geneA, complex_pair_tab_sites$geneB),]
  print(nrow(complex_pair_tab_sites))
  
  ## initiate value columns
  num_control <- vector(mode = "numeric", length = nrow(complex_pair_tab_sites))
  num_mut <- num_control; meddiff_mut <- num_control; p_mut <- num_control; meddiff_bottom_mut <- num_control; meddiff_upper_mut <- num_control
  
  for (enzyme in unique(complex_pair_tab_sites$geneA)) {
    mut_mat_en <- mut_mat[mut_mat$Hugo_Symbol == enzyme,]
    if (nrow(mut_mat_en) > 0){
      mut_partIDs <- partIDs_overlap[(mut_mat_en[,partIDs_overlap] != "") & (mut_mat_en[, partIDs_overlap] != "Silent")]
    } else {
      mut_partIDs <- NULL 
    }
    
    mut_cnv_partIDs <- mut_partIDs
    control_partIDs <- partIDs_overlap[!(partIDs_overlap %in% mut_cnv_partIDs)]
    
    print(paste0("enzyme:", enzyme))
    for (substrate in unique(complex_pair_tab_sites$geneB[complex_pair_tab_sites$geneA == enzyme])) {
      print(paste0("substrate:", substrate))
      
      for (site in unique(complex_pair_tab_sites$SUB_MOD_RSD[complex_pair_tab_sites$geneA == enzyme & complex_pair_tab_sites$geneB == substrate])) {
        print(paste0("site:", site))
        i <- (complex_pair_tab_sites$geneA == enzyme & complex_pair_tab_sites$geneB == substrate & complex_pair_tab_sites$SUB_MOD_RSD == site)
        
        pho_sub <- pho_data[pho_head$SUBSTRATE == substrate & pho_head$SUB_MOD_RSD == site,]
        pho_sub <- pho_sub[1,]
        control_pho <- pho_sub[, control_partIDs]
        control_pho <- control_pho[!is.na(control_pho)]
        num_control[i] <- length(control_pho)
        
        if (num_control[i] >= num_control_thres) {
          mut_pho <- pho_sub[, mut_partIDs]
          mut_pho <- mut_pho[!is.na(mut_pho)]
          num_mut[i] <- length(mut_pho)
          if (num_mut[i] > 0) {
            meddiff_mut[i] <- median(mut_pho) - median(control_pho)
            if (num_mut[i] >= num_genoalt_thres){
              stat <- wilcox.test(x = mut_pho, y = control_pho, conf.int = T)
              p_mut[i] <- stat$p.value
              meddiff_bottom_mut[i] <- stat$conf.int[1]
              meddiff_upper_mut[i] <- stat$conf.int[2]
            }
          }
        }
      }
    }
    print(which(i))
  }
  mut_cnv_tab <- complex_pair_tab_sites
  mut_cnv_tab <- cbind(mut_cnv_tab, data.frame(num_mut, num_control,
                                               p_mut, meddiff_bottom_mut, meddiff_upper_mut, meddiff_mut))
  mut_cnv_tab$cancer <- cancer
  mut_cnv_tab$pair_pro <- paste0(mut_cnv_tab$geneA, ":", mut_cnv_tab$geneB)
  mut_cnv_tab$pair <- paste0(mut_cnv_tab$geneA, ":", mut_cnv_tab$geneB, ":", mut_cnv_tab$SUB_MOD_RSD)

  write.table(x = mut_cnv_tab, file = paste0(makeOutDir(resultD = resultD), cancer, "_mut_cnv_tab.txt"), col.names = T, row.names = F, quote = F)
}

for (cancer in c("OV", "CO", "UCEC")) {
  if (cancer %in% cancers_sort) {
    pro_data <- loadProteinNormalizedTumor(cancer = cancer)
    pho_data <- loadPhosphositeNormalizedTumor(cancer = cancer)
    sampIDs <- colnames(pho_data)[!(colnames(pho_data) %in% c("Gene", "Phosphosite"))]
    pho_data <- pho_data[rowSums(!is.na(pho_data[, sampIDs])) >= (num_control_thres + num_genoalt_thres),]
    pho_head <- formatPhosphosite(pho_data$Phosphosite, pho_data$Gene)
    partIDs <- sampID2partID(sampleID_vector = sampIDs, sample_map = loadSampMap())
    names(pho_data) <- c("Gene", "Phosphosite", partIDs)
  }
  if (cancer == "UCEC") {
    ## input phosphorylation data
    pro_data <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1/UCEC_proteomics_V1.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
    pro_data <- data.frame(pro_data)
    colnames(pro_data)[1] <- "Gene"
    
    # pro_data <- loadProteinNormalizedTumor(cancer = cancer)
    pho_data <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1/UCEC_phosphoproteomics_site_level_V1.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
    pho_data <- data.frame(pho_data)
    # pho_data <- loadPhosphositeNormalizedTumor(cancer = cancer)
    pho_gdata <- read_delim("Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1/UCEC_phosphoproteomics_gene_level_V1.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
    pho_gdata <- data.frame(pho_gdata)
    colnames(pho_gdata)[1] <- "Gene"
    # pho_gdata <- loadPhosphoproteinNormalizedTumor(cancer = cancer)
    pho_head <- data.frame(str_split_fixed(string = pho_data$idx, pattern = "-", n = 2))
    colnames(pho_head) <- c("SUBSTRATE", "SUB_MOD_RSD")
    
    samples <- colnames(pro_data)[!(colnames(pro_data) %in% c("idx")) & !grepl(pattern = "Mix", x = (colnames(pro_data)))  & !grepl(pattern = "rep", x = (colnames(pro_data)))]
    samples <- samples[grepl(pattern = paste0('\\.', toupper(substr(x = sample_type, start = 1, stop = 1))), x = samples)]
    sampIDs <- samples
    tmp <- str_split_fixed(string = samples, pattern = "\\.", n = 3)
    partIDs <- paste0(tmp[,1], "-", tmp[,2])
    print(length(partIDs))
    pho_data <- pho_data[, sampIDs]
    colnames(pho_data) <- partIDs
  }
  
  ## input somatic mutation matrix
  mut_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_mutation_matrix/", cancer, "_somatic_mutation_in_enzyme_substrate.txt"), data.table = F)
  
  ## get the overlap of the participant IDs
  partIDs_overlap <- intersect(partIDs, colnames(mut_mat))
  print(length(partIDs_overlap))
  
  ## the number of mutations
  mut_count <- rowSums(x = ((mut_mat[, -1] != "") & (mut_mat[, -1] != "Silent"))); names(mut_count) <- mut_mat$Hugo_Symbol
  
  ## get the genes to test
  genes_num_genoalt_thresholded <- names(mut_count)[mut_count >= num_genoalt_thres]
  print(length(genes_num_genoalt_thresholded))
  
  ## initiate identifier columns
  complex_pair_tab_genoalt_thresholded <- complex_pair_tab[complex_pair_tab$geneA %in% genes_num_genoalt_thresholded,]
  print(nrow(complex_pair_tab_genoalt_thresholded))
  
  complex_pair_tab_sites <- merge(complex_pair_tab_genoalt_thresholded, pho_head[, c("SUBSTRATE", "SUB_MOD_RSD")], by.x = c("geneB"), by.y = c("SUBSTRATE"), all.y = T)
  complex_pair_tab_sites <- complex_pair_tab_sites[!is.na(complex_pair_tab_sites$geneA) & !is.na(complex_pair_tab_sites$SUB_MOD_RSD),]
  complex_pair_tab_sites <- complex_pair_tab_sites[order(complex_pair_tab_sites$geneA, complex_pair_tab_sites$geneB),]
  print(nrow(complex_pair_tab_sites))
  
  ## initiate value columns
  num_control <- vector(mode = "numeric", length = nrow(complex_pair_tab_sites))
  num_mut <- num_control; meddiff_mut <- num_control; p_mut <- num_control; meddiff_bottom_mut <- num_control; meddiff_upper_mut <- num_control
  
  for (enzyme in unique(complex_pair_tab_sites$geneA)) {
    mut_mat_en <- mut_mat[mut_mat$Hugo_Symbol == enzyme,]
    if (nrow(mut_mat_en) > 0){
      mut_partIDs <- partIDs_overlap[(mut_mat_en[,partIDs_overlap] != "") & (mut_mat_en[, partIDs_overlap] != "Silent")]
    } else {
      mut_partIDs <- NULL 
    }
    
    mut_cnv_partIDs <- mut_partIDs
    control_partIDs <- partIDs_overlap[!(partIDs_overlap %in% mut_cnv_partIDs)]
    
    print(paste0("enzyme:", enzyme))
    for (substrate in unique(complex_pair_tab_sites$geneB[complex_pair_tab_sites$geneA == enzyme])) {
      print(paste0("substrate:", substrate))
      
      for (site in unique(complex_pair_tab_sites$SUB_MOD_RSD[complex_pair_tab_sites$geneA == enzyme & complex_pair_tab_sites$geneB == substrate])) {
        print(paste0("site:", site))
        i <- (complex_pair_tab_sites$geneA == enzyme & complex_pair_tab_sites$geneB == substrate & complex_pair_tab_sites$SUB_MOD_RSD == site)
        
        pho_sub <- pho_data[pho_head$SUBSTRATE == substrate & pho_head$SUB_MOD_RSD == site,]
        pho_sub <- pho_sub[1,]
        control_pho <- pho_sub[, control_partIDs]
        control_pho <- control_pho[!is.na(control_pho)]
        num_control[i] <- length(control_pho)
        
        if (num_control[i] >= num_control_thres) {
          mut_pho <- pho_sub[, mut_partIDs]
          mut_pho <- mut_pho[!is.na(mut_pho)]
          num_mut[i] <- length(mut_pho)
          if (num_mut[i] > 0) {
            meddiff_mut[i] <- median(mut_pho) - median(control_pho)
            if (num_mut[i] >= num_genoalt_thres){
              stat <- wilcox.test(x = mut_pho, y = control_pho, conf.int = T)
              p_mut[i] <- stat$p.value
              meddiff_bottom_mut[i] <- stat$conf.int[1]
              meddiff_upper_mut[i] <- stat$conf.int[2]
            }
          }
        }
      }
    }
    print(which(i))
  }
  mut_cnv_tab <- complex_pair_tab_sites
  mut_cnv_tab <- cbind(mut_cnv_tab, data.frame(num_mut, num_control,
                                               p_mut, meddiff_bottom_mut, meddiff_upper_mut, meddiff_mut))
  mut_cnv_tab$cancer <- cancer
  mut_cnv_tab$pair_pro <- paste0(mut_cnv_tab$geneA, ":", mut_cnv_tab$geneB)
  mut_cnv_tab$pair <- paste0(mut_cnv_tab$geneA, ":", mut_cnv_tab$geneB, ":", mut_cnv_tab$SUB_MOD_RSD)
  
  ## filter mutation-impacted kinase-substrate pairs
  mut_tab <- mut_cnv_tab
  mut_tab <- mut_tab[mut_tab$p_mut > 0,]
  
  ## format the mutation table
  mut_tab.f <- mut_tab[, c(id_vars, paste0(paste0(prefix_vars, "", sep = "_"), "mut"), "pair")]
  colnames(mut_tab.f) <- c(id_vars, prefix_vars, "pair")
  mut_tab.f$genoalt_type <- "mut"
  
  mut_cnv_can <- mut_tab.f
  mut_cnv_cans <- rbind(mut_cnv_cans, mut_cnv_can)
  mut_cnv_cans <- unique(mut_cnv_cans)
  
  write.table(x = mut_cnv_tab, file = paste0(makeOutDir(resultD = resultD), cancer, "_mut_cnv_tab.txt"), col.names = T, row.names = F, quote = F)
}


# parse all cancers -------------------------------------------------------
mut_cnv_cans <- NULL
for (cancer in c("BRCA", "OV", "CO", "UCEC")) {
  mut_cnv_tab <- fread(file = paste0(makeOutDir(resultD = resultD), cancer, "_mut_cnv_tab.txt"), data.table = F)
  ## filter mutation-impacted kinase-substrate pairs
  mut_tab <- mut_cnv_tab
  mut_tab <- mut_tab[mut_tab$p_mut > 0,]
  
  ## format the mutation table
  mut_tab.f <- mut_tab[, c(id_vars, paste0(paste0(prefix_vars, "", sep = "_"), "mut"), "pair", "cancer")]
  colnames(mut_tab.f) <- c(id_vars, prefix_vars, "pair", "cancer")
  mut_tab.f$genoalt_type <- "mut"
  
  mut_cnv_can <- mut_tab.f
  mut_cnv_cans <- rbind(mut_cnv_cans, mut_cnv_can)
  mut_cnv_cans <- unique(mut_cnv_cans)
}
mut_cnv_cans$driver_gene_type.en <- ""
mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$geneA %in% oncogenes] <- "oncogene"
mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$geneA %in% tsgs] <- "tsg"
write.table(x = mut_cnv_cans, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_sig_cans.txt"), row.names = F, quote = F, sep = "\t")


