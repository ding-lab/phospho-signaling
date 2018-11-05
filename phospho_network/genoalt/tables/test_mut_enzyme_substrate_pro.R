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
id_vars <- c("GENE", "SUB_GENE", "SUB_MOD_RSD", "num_control", "cancer")

# inputs ------------------------------------------------------------------
## input enzyme-substrate table
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source != "NetKIN" | (ptms_site_pairs_sup$Source == "NetKIN" & ptms_site_pairs_sup$networkin_score >= 5),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)
enzyme_sub <- unique(ptms_site_pairs_sup[, c("GENE", "SUB_GENE")])
enzyme_sub_rev <- data.frame(GENE = enzyme_sub$SUB_GENE, SUB_GENE = enzyme_sub$GENE)
enzyme_sub_cis <- data.frame(GENE = unique(c(as.vector(enzyme_sub$GENE), as.vector(enzyme_sub$SUB_GENE))), SUB_GENE = unique(c(as.vector(enzyme_sub$GENE), as.vector(enzyme_sub$SUB_GENE))))
enzyme_sub <- unique(rbind(enzyme_sub, enzyme_sub_rev, enzyme_sub_cis))

for (cancer in c("BRCA", "OV", "CO", "UCEC")) {
  if (cancer %in% cancers_sort) {
    pro_data <- loadProteinNormalizedTumor(cancer = cancer)
    pho_data <- loadPhosphositeNormalizedTumor(cancer = cancer)
    sampIDs <- colnames(pho_data)[!(colnames(pho_data) %in% c("Gene", "Phosphosite"))]
    pho_data <- pho_data[rowSums(!is.na(pho_data[, sampIDs])) >= (num_control_thres + num_genoalt_thres),]
    pho_head <- formatPhosphosite(pho_data$Phosphosite, pho_data$Gene)
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
    pho_head <- data.frame(str_split_fixed(string = pho_data$idx, pattern = "-", n = 2))
    colnames(pho_head) <- c("SUBSTRATE", "SUB_MOD_RSD")
    
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
  enzyme_sub_genoalt_thresholded <- enzyme_sub[enzyme_sub$GENE %in% genes_num_genoalt_thresholded,]
  print(nrow(enzyme_sub_genoalt_thresholded))
  
  pro_head <- data.frame(SUBSTRATE = pro_data$Gene, SUB_MOD_RSD = "protein")
  enzyme_sub_sites <- merge(enzyme_sub_genoalt_thresholded, pro_head, by.x = c("SUB_GENE"), by.y = c("SUBSTRATE"), all.y = T)
  enzyme_sub_sites <- enzyme_sub_sites[!is.na(enzyme_sub_sites$GENE) & !is.na(enzyme_sub_sites$SUB_MOD_RSD),]
  enzyme_sub_sites <- enzyme_sub_sites[order(enzyme_sub_sites$GENE, enzyme_sub_sites$SUB_GENE),]
  
  ## initiate value columns
  num_control <- vector(mode = "numeric", length = nrow(enzyme_sub_sites))
  num_mut <- num_control; meddiff_mut <- num_control; p_mut <- num_control; meddiff_bottom_mut <- num_control; meddiff_upper_mut <- num_control
  
  for (enzyme in unique(enzyme_sub_sites$GENE)) {
    mut_mat_en <- mut_mat[mut_mat$Hugo_Symbol == enzyme,]
    if (nrow(mut_mat_en) > 0){
      mut_partIDs <- partIDs_overlap[(mut_mat_en[,partIDs_overlap] != "") & (mut_mat_en[, partIDs_overlap] != "Silent")]
    } else {
      mut_partIDs <- NULL 
    }
    
    mut_cnv_partIDs <- mut_partIDs
    control_partIDs <- partIDs_overlap[!(partIDs_overlap %in% mut_cnv_partIDs)]
    
    print(paste0("enzyme:", enzyme))
    for (substrate in unique(enzyme_sub_sites$SUB_GENE[enzyme_sub_sites$GENE == enzyme])) {
      print(paste0("substrate:", substrate))
      i <- (enzyme_sub_sites$GENE == enzyme & enzyme_sub_sites$SUB_GENE == substrate)
      
      pro_sub <- pro_data[pro_data$Gene == substrate,]
      pro_sub <- pro_sub[1,]
      control_pro <- pro_sub[, control_partIDs]
      control_pro <- control_pro[!is.na(control_pro)]
      num_control[i] <- length(control_pro)
      
      if (num_control[i] >= num_control_thres) {
        mut_pro <- pro_sub[, mut_partIDs]
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
    print(which(i))
  }
  mut_cnv_tab <- enzyme_sub_sites
  mut_cnv_tab <- cbind(mut_cnv_tab, data.frame(num_mut, num_control,
                                               p_mut, meddiff_bottom_mut, meddiff_upper_mut, meddiff_mut))
  mut_cnv_tab$cancer <- cancer
  mut_cnv_tab$pair_pro <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE)
  mut_cnv_tab$pair <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE, ":", mut_cnv_tab$SUB_MOD_RSD)
  
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
  mut_tab.f <- mut_tab[, c(id_vars, paste0(paste0(prefix_vars, "", sep = "_"), "mut"))]
  colnames(mut_tab.f) <- c(id_vars, prefix_vars)
  mut_tab.f$genoalt_type <- "mut"
  
  mut_cnv_can <- mut_tab.f
  mut_cnv_cans <- rbind(mut_cnv_cans, mut_cnv_can)
  mut_cnv_cans <- unique(mut_cnv_cans)
}
mut_cnv_cans$driver_gene_type.en <- ""
mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$GENE %in% oncogenes] <- "oncogene"
mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$GENE %in% tsgs] <- "tsg"
write.table(x = mut_cnv_cans, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_sig_cans.txt"), row.names = F, quote = F, sep = "\t")

