# Yige Wu @ WashU 2018 Aug
## for each kinase-substrate pairs, compare the substrate phosphorylation level of kinase-mutated tumors, kinase-amplified tumors, kinase-deleted tumors against other tumors
## using wilcox test
## reference: https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/

# source ------------------------------------------------------------------
source("Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
num_control_thres <- 4
num_genoalt_thres <- 4
sample_type <- "tumor"
sig_thres <- 0.05
id_vars <- c("GENE", "SUB_GENE", "SUB_MOD_RSD", "cancer")
prefix_vars <- c("p", "meddiff_bottom", "meddiff_upper", "meddiff", "num")

# inputs ------------------------------------------------------------------
## input enzyme-substrate table
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source != "NetKIN" | (ptms_site_pairs_sup$Source == "NetKIN" & ptms_site_pairs_sup$networkin_score >= 5),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)
enzyme_sub <- unique(ptms_site_pairs_sup[, c("GENE", "SUB_GENE")])
enzyme_sub_cis <- data.frame(GENE = unique(c(as.vector(enzyme_sub$GENE), as.vector(enzyme_sub$SUB_GENE))), SUB_GENE = unique(c(as.vector(enzyme_sub$GENE), as.vector(enzyme_sub$SUB_GENE))))
enzyme_sub <- enzyme_sub_cis
print(nrow(enzyme_sub))

mut_cnv_cans <- NULL
count_enzyme_cans <- NULL
mut_cnv_sig_cans <- NULL

# UCEC --------------------------------------------------------------------
for (cancer in c("UCEC")) {
  pro_data <- read_delim("Ding_Lab/Projects_Current/CPTAC/cptac_shared/Endometrium_CPTAC3/01_Data_tables/UCEC_V1/UCEC_proteomics_V1.cct", "\t", escape_double = FALSE, trim_ws = TRUE)
  pro_data <- data.frame(pro_data)
  colnames(pro_data)[1] <- "Gene"
  # pro_data <- loadProteinNormalizedTumor(cancer = cancer)

  samples <- colnames(pro_data)[!(colnames(pro_data) %in% c("idx")) & !grepl(pattern = "Mix", x = (colnames(pro_data)))  & !grepl(pattern = "rep", x = (colnames(pro_data)))]
  samples <- samples[grepl(pattern = paste0('\\.', toupper(substr(x = sample_type, start = 1, stop = 1))), x = samples)]
  sampIDs <- samples
  tmp <- str_split_fixed(string = samples, pattern = "\\.", n = 3)
  partIDs <- paste0(tmp[,1], "-", tmp[,2])
  print(length(partIDs))
  pro_data <- pro_data[, c("Gene", sampIDs)]
  colnames(pro_data) <- c("Gene", partIDs)
  
  ## input somatic mutation matrix
  mut_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_mutation_matrix_UCEC/", cancer, "_somatic_mutation_in_enzyme_substrate.txt"), data.table = F)
  
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
  
  ## initiate value columns
  num_control <- vector(mode = "numeric", length = nrow(enzyme_sub))
  num_mut <- num_control; meddiff_mut <- num_control; p_mut <- num_control; meddiff_bottom_mut <- num_control; meddiff_upper_mut <- num_control
  enzymes <- vector(mode = "character", length = nrow(enzyme_sub))
  substrates <- vector(mode = "character", length = nrow(enzyme_sub))
  
  for (enzyme in unique(enzyme_sub$GENE)) {
    mut_mat_en <- mut_mat[mut_mat$Hugo_Symbol == enzyme,]
    if (nrow(mut_mat_en) > 0){
      mut_partIDs <- partIDs_overlap[(mut_mat_en[,partIDs_overlap] != "") & (mut_mat_en[, partIDs_overlap] != "Silent")]
    } else {
      mut_partIDs <- NULL 
    }
    
    mut_cnv_partIDs <- mut_partIDs
    control_partIDs <- partIDs_overlap[!(partIDs_overlap %in% mut_cnv_partIDs)]
    
    print(paste0("enzyme:", enzyme))
    
    for (substrate in unique(enzyme_sub$SUB_GENE[enzyme_sub$GENE == enzyme])) {
      print(paste0("substrate:", substrate))
      
      i <- (enzyme_sub$GENE == enzyme & enzyme_sub$SUB_GENE == substrate)
      enzymes[i] <- enzyme
      substrates[i] <- substrate
      
      control_pro <- pro_data[pro_data$Gene == enzyme, control_partIDs]
      control_pro <- control_pro[!is.na(control_pro)]
      num_control[i] <- length(control_pro)
      
      if (num_control[i] >= num_control_thres) {
        mut_pro <- pro_data[pro_data$Gene == enzyme, mut_partIDs]
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
  }
  mut_cnv_tab <- enzyme_sub
  mut_cnv_tab$SUB_MOD_RSD  <- "protein"
  mut_cnv_tab <- cbind(mut_cnv_tab, data.frame(num_mut, num_control,
                                               p_mut, meddiff_bottom_mut, meddiff_upper_mut, meddiff_mut))
  mut_cnv_tab$cancer <- cancer
  
  mut_cnv_tab$pair_pro <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE)
  mut_cnv_tab$pair <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE, ":", mut_cnv_tab$SUB_MOD_RSD)
  mut_cnv_tab <- merge(mut_cnv_tab, unique(ptms_site_pairs_sup[, c("GENE", "enzyme_type")]), all.x = T)
  
  ##filter out bad annotated enzyme type
  mut_cnv_tab <- mut_cnv_tab[!(mut_cnv_tab$GENE == "PTPRG" & mut_cnv_tab$enzyme_type == "kinase"),]
  
  ## annotate cis and trans
  mut_cnv_tab$SELF <- "trans"
  mut_cnv_tab$SELF[as.vector(mut_cnv_tab$GENE) == as.vector(mut_cnv_tab$SUB_GENE)] <- "cis"
  
  ## filter mutation-impacted kinase-substrate pairs
  mut_tab <- mut_cnv_tab
  mut_tab <- mut_tab[mut_tab$p_mut > 0,]
  
  sig_mut_tab <- mut_tab[mut_tab$p_mut < sig_thres,]
  sig_mut_tab <- sig_mut_tab[(sig_mut_tab$meddiff_bottom_mut > 0) == (sig_mut_tab$meddiff_upper_mut > 0),]
  sig_mut_tab <- sig_mut_tab[order(sig_mut_tab$p_mut),]
  
  ## count the occurence of each enzyme
  count_enzyme_mut <- data.frame(table(sig_mut_tab[, c("GENE", "SELF")]))
  if (nrow(count_enzyme_mut) > 0){
    count_enzyme_mut$genoalt_type <- "mutation"
    
    ## format the mutation table
    mut_tab.f <- mut_tab[, c(id_vars, paste0(paste0(prefix_vars, "", sep = "_"), "mut"), "pair")]
    colnames(mut_tab.f) <- c(id_vars, prefix_vars, "pair")
    mut_tab.f$genoalt_type <- "mutation"
    sig_mut_tab.f <- mut_tab.f[mut_tab.f$pair %in% sig_mut_tab$pair,]
  }
  
  count_enzyme_can <- count_enzyme_mut
  colnames(count_enzyme_can) <- c("GENE", "SELF", "Freq", "genoalt_type")
  count_enzyme_can <- merge(count_enzyme_can, unique(ptms_site_pairs_sup[, c("GENE", "enzyme_type")]), all.x = T)
  count_enzyme_can$cancer <- cancer
  count_enzyme_cans <- rbind(count_enzyme_cans, count_enzyme_can)
  
  mut_cnv_sig_can <- sig_mut_tab.f
  mut_cnv_sig_cans <- rbind(mut_cnv_sig_cans, mut_cnv_sig_can)
  mut_cnv_sig_cans <- unique(mut_cnv_sig_cans)
  
  mut_cnv_can <- mut_tab.f
  mut_cnv_cans <- rbind(mut_cnv_cans, mut_cnv_can)
  mut_cnv_cans <- unique(mut_cnv_cans)
  
  write.table(x = mut_cnv_tab, file = paste0(makeOutDir(resultD = resultD), cancer, "_mut_cnv_tab.txt"), row.names = F, quote = F)
  write.table(x = mut_cnv_can, file = paste0(makeOutDir(resultD = resultD), cancer, "_mut_cnv_tab_formatted.txt"), row.names = F, quote = F, sep = "\t")
}

count_enzyme_cans <- data.frame(count_enzyme_cans)
count_enzyme_cans$driver_gene_type.en <- ""
count_enzyme_cans$driver_gene_type.en[count_enzyme_cans$GENE %in% oncogenes] <- "oncogene"
count_enzyme_cans$driver_gene_type.en[count_enzyme_cans$GENE %in% tsgs] <- "tsg"

mut_cnv_sig_cans$SELF <- ifelse(as.vector(mut_cnv_sig_cans$GENE) == as.vector(mut_cnv_sig_cans$SUB_GENE), "cis", "trans")

mut_cnv_cans <- NULL
for (cancer in c("UCEC")) {
  mut_cnv_can <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/genoalt/tables/test_mut_cna_impact_cis_pro_cptac3/", cancer, "_mut_cnv_tab_formatted.txt"), data.table = F)
  mut_cnv_cans <- rbind(mut_cnv_cans, mut_cnv_can)
}
mut_cnv_cans$SELF <- ifelse(as.vector(mut_cnv_cans$GENE) == as.vector(mut_cnv_cans$SUB_GENE), "cis", "trans")
mut_cnv_cans$driver_gene_type.en <- ""
mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$GENE %in% oncogenes] <- "oncogene"
mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$GENE %in% tsgs] <- "tsg"

write.table(x = count_enzyme_cans, file = paste0(makeOutDir(resultD = resultD), "count_enzyme_cans_sig_thres", sig_thres, ".txt"), row.names = F, quote = F, sep = "\t")
write.table(x = mut_cnv_sig_cans, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_sig_cans_sig_thres", sig_thres, ".txt"), row.names = F, quote = F, sep = "\t")
write.table(x = mut_cnv_cans, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_sig_cans.txt"), row.names = F, quote = F, sep = "\t")


