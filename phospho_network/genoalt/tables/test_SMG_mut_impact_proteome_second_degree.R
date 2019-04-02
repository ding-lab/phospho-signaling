# Yige Wu @ WashU 2018 Aug
## test mutation impact on protein/phosphorylation within kinase-substrate pairs or protein complex pairs
## using wilcox test
## reference: https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/

# source ------------------------------------------------------------------
setwd(dir = "~/Box Sync/")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
sample_type <- "tumor"
id_vars <- c("GENE", "SUB_GENE", "SUB_MOD_RSD", "cancer")
prefix_vars <- c("p", "num", "meddiff_bottom", "meddiff_upper", "meddiff")
num_control_thres <- 5
num_genoalt_thres <- 5

# cancer types to process -------------------------------------------------
cancers2process <- c("BRCA", "CO", "OV", "UCEC", "CCRCC", "LIHC")

# expression data to process ----------------------------------------------
# affected_exp_type <- "PHO"
# affected_exp_type <- "PRO"
# affected_exp_type <- "RNA"
affected_exp_types2process <- c("RNA", "PRO", "PHO")


# variant types to test ---------------------------------------------------
## include all nonsynonymous mutations or just missense/truncation
# variant_class <- "truncation"
# variant_class <- "not_silent"
# variant_class <- "missense"
variant_classes2process <- c("truncation", "missense", "not_silent")
# variant_classes2process <- c("not_silent")


# input previous result ---------------------------------------------------
mut_impact_tab_previous <- fread(input = "./cptac2p/analysis_results/phospho_network/genoalt/tables/test_SMG_mut_impact_proteome/SMG_mut_impact_tab.txt", data.table = F)
mut_impact_tab_sig <- mut_impact_tab_previous %>%
  filter(fdr_by_gene < 0.1)

# input protein pair table ------------------------------------------------
pair_tab_annotated <- fread(input = "./cptac2p/analysis_results/dependencies/tables/compile_protein_pair_table/protein_pair_table.txt", data.table = F)
pair_tab <- mut_impact_tab_sig %>%
  select(GENE, SUB_GENE) %>%
  unique()
pair_tab2merge <- pair_tab_annotated[pair_tab_annotated$GENE %in% pair_tab$SUB_GENE, c("GENE", "SUB_GENE")]
pair_tab <- merge(pair_tab, pair_tab2merge, by.x = c("SUB_GENE"), by.y = c("GENE"), all.x = T)
pair_tab_previous <- data.frame(pair_pro_previous = paste0(pair_tab$GENE, ":", pair_tab$SUB_GENE), GENE = pair_tab$GENE, SUB_GENE = pair_tab$SUB_GENE.y)

pair_tab$SUB_GENE <- pair_tab$SUB_GENE.y
pair_tab$SUB_GENE.y <- NULL
pair_tab$pair_pro <- paste0(pair_tab$GENE, ":", pair_tab$SUB_GENE)
pair_tab <- pair_tab %>%
  filter(!(pair_pro %in% mut_impact_tab_previous$pair_pro))
pair_tab$pair_pro <- NULL
pair_tab <- unique(pair_tab)

## clean up
bad_strings <- c("DNA damage", "")


# Business  ---------------------------------------------------------------
for (affected_exp_type in affected_exp_types2process) {
  for (variant_class in variant_classes2process) {
    for (cancer in cancers2process) {
      fn <- paste0(makeOutDir(resultD = resultD), cancer, "_", variant_class,"_mut_impact_", affected_exp_type , "_tab.txt")
      if (!file.exists(fn)) {
        if (cancer %in% c("BRCA", "OV", "CO")) {
          pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
          pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
        } else if (cancer == "UCEC") {
          pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
          pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
        } else if (cancer == "CCRCC") {
          pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
          pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
        } else if (cancer == "LIHC"){
          pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
          pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
        }
        if (affected_exp_type == "PRO") {
          affected_exp_data <- pro_tab
          affected_exp_head <- data.frame(SUBSTRATE = affected_exp_data$Gene)
        }
        if (affected_exp_type == "RNA") {
          affected_exp_data <- loadRNA(cancer = cancer)
          affected_exp_data <- affected_exp_data[rowSums((affected_exp_data == 0 & !is.na(affected_exp_data)) | is.na(affected_exp_data)) < (as.numeric(ncol(affected_exp_data)) - num_genoalt_thres - num_control_thres -1),]
          affected_exp_head <- data.frame(SUBSTRATE = affected_exp_data$gene)
        }
        if (affected_exp_type == "PHO") {
          affected_exp_data <- pho_tab
          affected_exp_head <- data.frame(SUBSTRATE = affected_exp_data$Gene, SUB_MOD_RSD = affected_exp_data$Phosphosite)
        }
        
        ## input mutation matrix
        maf <- loadMaf(cancer = cancer, maf_files = maf_files)
        mut_mat <- generate_somatic_mutation_matrix(pair_tab = pair_tab$GENE, maf = maf)
        
        ## get the overlap of the participant IDs
        partIDs <- colnames(affected_exp_data)
        partIDs_overlap <- intersect(partIDs, colnames(mut_mat))
        print(length(partIDs_overlap))
        
        ## the number of mutations
        if (variant_class == "not_silent") {
          mut_count <- rowSums(x = ((mut_mat[, -1] != "") & (mut_mat[, -1] != "Silent"))); 
        }
        if (variant_class == "missense") {
          mut_count <- sapply(1:nrow(mut_mat), FUN = function(i, mut_mat) sum(grepl(x = mut_mat[i,], pattern = "Missense_Mutation")), mut_mat[,-1])
        }
        if (variant_class == "truncation") {
          mut_count <- sapply(1:nrow(mut_mat), FUN = function(i, mut_mat) sum(grepl(x = mut_mat[i,], pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del")), mut_mat[,-1])
        }
        names(mut_count) <- mut_mat$Hugo_Symbol
        
        ## get the genes to test
        genes_num_genoalt_thresholded <- names(mut_count)[mut_count >= num_genoalt_thres]
        print(length(genes_num_genoalt_thresholded))
        
        ## initiate identifier columns
        pair_tab_genoalt_thresholded <- pair_tab[pair_tab$GENE %in% genes_num_genoalt_thresholded,]
        print(nrow(pair_tab_genoalt_thresholded))
        
        if (affected_exp_type %in% c("PRO", "RNA")) {
          pair_tab_sites <- pair_tab_genoalt_thresholded
          pair_tab_sites$SUB_MOD_RSD <- affected_exp_type
        } else {
          pair_tab_sites <- merge(pair_tab_genoalt_thresholded, affected_exp_head[, c("SUBSTRATE", "SUB_MOD_RSD")], by.x = c("SUB_GENE"), by.y = c("SUBSTRATE"))
        }
        pair_tab_sites <- pair_tab_sites[!is.na(pair_tab_sites$GENE) & !is.na(pair_tab_sites$SUB_MOD_RSD),]
        pair_tab_sites <- pair_tab_sites[order(pair_tab_sites$GENE, pair_tab_sites$SUB_GENE),]
        print(nrow(pair_tab_sites))
        
        ## initiate value columns
        num_control <- vector(mode = "numeric", length = nrow(pair_tab_sites))
        num_mut <- num_control; meddiff_mut <- num_control; p_mut <- num_control; meddiff_bottom_mut <- num_control; meddiff_upper_mut <- num_control
        
        for (enzyme in unique(pair_tab_sites$GENE)) {
          mut_mat_en <- mut_mat[mut_mat$Hugo_Symbol == enzyme,]
          if (nrow(mut_mat_en) > 0){
            if (variant_class == "not_silent") {
              mut_partIDs <- partIDs_overlap[(mut_mat_en[,partIDs_overlap] != "") & (mut_mat_en[, partIDs_overlap] != "Silent")]
            }
            if (variant_class == "missense") {
              mut_partIDs <- partIDs_overlap[grepl(x = mut_mat_en[, partIDs_overlap], pattern = "Missense_Mutation")]
            }
            if (variant_class == "truncation") {
              mut_partIDs <- partIDs_overlap[grepl(x = mut_mat_en[, partIDs_overlap], pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del")]
            }
            
          } else {
            mut_partIDs <- NULL 
          }
          
          mut_cnv_partIDs <- mut_partIDs
          control_partIDs <- partIDs_overlap[(mut_mat_en[,partIDs_overlap] == "") | (mut_mat_en[, partIDs_overlap] == "Silent")]
          
          print(paste0("enzyme:", enzyme))
          for (substrate in unique(pair_tab_sites$SUB_GENE[pair_tab_sites$GENE == enzyme])) {
            print(paste0("substrate:", substrate))
            
            for (site in unique(pair_tab_sites$SUB_MOD_RSD[pair_tab_sites$GENE == enzyme & pair_tab_sites$SUB_GENE == substrate])) {
              print(paste0("site:", site))
              
              if (affected_exp_type  %in% c("PRO", "RNA")) {
                affected_exp_altered <- affected_exp_data[affected_exp_head$SUBSTRATE == substrate,]
                i <- (pair_tab_sites$GENE == enzyme & pair_tab_sites$SUB_GENE == substrate & pair_tab_sites$SUB_MOD_RSD == site)
              } else {
                affected_exp_altered <- affected_exp_data[affected_exp_head$SUBSTRATE == substrate & affected_exp_head$SUB_MOD_RSD == site,]
                i <- (pair_tab_sites$GENE == enzyme & pair_tab_sites$SUB_GENE == substrate & pair_tab_sites$SUB_MOD_RSD == site)
                
              }
              affected_exp_altered <- affected_exp_altered[1,]
              affected_exp_control <- affected_exp_altered[, control_partIDs]
              affected_exp_control <- affected_exp_control[!is.na(affected_exp_control)]
              num_control[i] <- length(affected_exp_control)
              
              if (num_control[i] >= num_control_thres) {
                mut_pho <- affected_exp_altered[, mut_partIDs]
                mut_pho <- mut_pho[!is.na(mut_pho)]
                num_mut[i] <- length(mut_pho)
                if (num_mut[i] > 0) {
                  meddiff_mut[i] <- median(mut_pho) - median(affected_exp_control)
                  if (num_mut[i] >= num_genoalt_thres){
                    stat <- wilcox.test(x = mut_pho, y = affected_exp_control, conf.int = T)
                    p_mut[i] <- stat$p.value
                    meddiff_bottom_mut[i] <- stat$conf.int[1]
                    meddiff_upper_mut[i] <- stat$conf.int[2]
                  }
                }
              }
              print(which(i))
            }
          }
          
        }
        mut_cnv_tab <- pair_tab_sites
        mut_cnv_tab <- cbind(mut_cnv_tab, 
                             data.frame(num = num_mut, 
                                        num_control = num_control,
                                        p = p_mut, 
                                        meddiff_bottom = meddiff_bottom_mut, 
                                        meddiff_upper = meddiff_upper_mut, 
                                        meddiff = meddiff_mut))
        mut_cnv_tab$cancer <- cancer
        mut_cnv_tab$pair <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE, ":", mut_cnv_tab$SUB_MOD_RSD)
        ## annotate cis and trans
        mut_cnv_tab$SELF <- "trans"
        mut_cnv_tab$SELF[as.vector(mut_cnv_tab$GENE) == as.vector(mut_cnv_tab$SUB_GENE)] <- "cis"
        mut_cnv_tab$genoalt_type <- "mut"
        
        ## clean up
        mut_cnv_tab <- mut_cnv_tab[!(mut_cnv_tab$GENE %in% bad_strings) & !(mut_cnv_tab$SUB_GENE %in% bad_strings),]
        write.table(x = mut_cnv_tab, file = fn, col.names = T, row.names = F, quote = F, sep = "\t")
        
      }
    }
  }
}

mut_cnv_cans <- NULL
for (affected_exp_type in c("PRO", "PHO", "RNA")) {
  for (variant_class in c("truncation", "missense", "not_silent")) {
    for (cancer in c("BRCA", "OV", "CO", "UCEC", "CCRCC", "LIHC")) {
      for (genoalt_type in c("mut")) {
        fn <- paste0(makeOutDir(resultD = resultD), cancer, "_", variant_class,"_mut_impact_", affected_exp_type , "_tab.txt")
        mut_cnv_can <- fread(input = fn, data.table = F, sep = "\t")
        mut_cnv_can$variant_class <- variant_class
        mut_cnv_can$affected_exp_type <- affected_exp_type
        mut_cnv_can$cancer <- cancer
        mut_cnv_cans <- rbind(mut_cnv_cans, mut_cnv_can)
      }
    }
  }
}
if (length(unique(mut_cnv_cans$cancer)) < length(cancers2process)) {
  stop("cancer type not enough!")
}
if (length(unique(mut_cnv_cans$variant_class)) < 3) {
  stop("variant_class not enough!")
}
if (length(unique(mut_cnv_cans$affected_exp_type)) < 3) {
  stop("affected_exp_type not enough!")
}
mut_cnv_cans <- mut_cnv_cans[mut_cnv_cans$num >= num_genoalt_thres,]
mut_cnv_cans <- unique(mut_cnv_cans)
mut_cnv_cans <- merge(mut_cnv_cans, pair_tab_previous, by = c("GENE", "SUB_GENE"), all.x = T)

mut_cnv_cans$fdr <- FDR_by_id_columns(p_vector = mut_cnv_cans$p, 
                                      id_columns = c("SELF", "cancer", "variant_class", "affected_exp_type"), 
                                      df = mut_cnv_cans)
mut_cnv_cans$fdr_by_gene <- FDR_by_id_columns(p_vector = mut_cnv_cans$p, 
                                              id_columns = c("SELF", "cancer", "variant_class", "affected_exp_type", "GENE"), 
                                              df = mut_cnv_cans)
write.table(x = mut_cnv_cans, file = paste0(makeOutDir(resultD = resultD), "SMG_mut_impact_tab.txt"), sep = "\t", row.names = F, quote = F)

