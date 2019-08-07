# Yige Wu @ WashU 2018 Aug
## test mutation impact on protein/phosphorylation within kinase-substrate pairs or protein complex pairs
## using wilcox test
## reference: https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/
## take out the samples deleted the altered gene

# source ------------------------------------------------------------------
code_top_dir <- "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/"
path2phospho_network_shared <- paste0(code_top_dir, "phospho_network/phospho_network_shared.R")
source(path2phospho_network_shared)

# set variables -----------------------------------------------------------
sample_type <- "tumor"
id_vars <- c("GENE", "SUB_GENE", "SUB_MOD_RSD", "cancer")
prefix_vars <- c("p", "num", "meddiff_bottom", "meddiff_upper", "meddiff")
num_control_thres <- 4
num_genoalt_thres <- 4
genes_altered <- c("PBRM1", "SETD2", "BAP1")
bad_strings <- c("DNA damage", "")

# cancer types to process -------------------------------------------------
# cancers2process <- c("BRCA", "OV", "CO")
# cancers2process <- c("CCRCC")
cancers2process <- c("CCRCC")
# cancers2process <- c("LIHC")
# cancers2process <- c("CO")

# expression data to process ----------------------------------------------
# affected_exp_type <- "PHO"
# affected_exp_type <- "PRO"
# affected_exp_type <- "RNA"
affected_exp_types2process <- c("RNA", "PRO", "PHO")
# affected_exp_types2process <- c("PRO", "PHO")
# affected_exp_types2process <- c("PRO", "PHO")

# variant types to test ---------------------------------------------------
## include all nonsynonymous mutations or just missense/truncation
variant_classes2process <- c("missense", "truncation", "not_silent")
variant_classes2process <- c("not_silent")

# Business  ---------------------------------------------------------------
for (gene_altered_tmp in genes_altered) {
  for (affected_exp_type in affected_exp_types2process) {
    for (variant_class in variant_classes2process) {
      for (cancer in cancers2process) {
        fn <- paste0(makeOutDir(), cancer, "_", variant_class, "_", gene_altered_tmp, "_mut_impact_", affected_exp_type , "_tab.txt")
        if (!file.exists(fn)) {
          # inputs ------------------------------------------------------------------
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
          ## construct the pairs to test
          pair_tab <- data.frame(GENE = gene_altered_tmp, SUB_GENE = unique(affected_exp_head$SUBSTRATE))
          
          ## input mutation matrix
          maf <- loadMaf(cancer = cancer, maf_files = maf_files)
          mut_mat <- generate_somatic_mutation_matrix(pair_tab = pair_tab, maf = maf)
          
          ## get the overlap of the participant IDs
          partIDs <- colnames(affected_exp_data)[!(colnames(affected_exp_data) %in% c("gene", "Gene", "Phosphosite", "Peptide_ID"))]
          partIDs_in_mut <- intersect(partIDs, colnames(mut_mat))
          print(length(partIDs_in_mut))
          
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
          
          ## initiate identifier columns
          pair_tab_genoalt_thresholded <- pair_tab
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
          
          # get the patient ids with mutations -------------------------------------
          
          for (enzyme in unique(pair_tab_sites$GENE)) {
            mut_mat_en <- mut_mat[mut_mat$Hugo_Symbol == enzyme,]
            if (nrow(mut_mat_en) > 0){
              if (variant_class == "not_silent") {
                mut_partIDs <- partIDs_in_mut[(mut_mat_en[,partIDs_in_mut] != "") & (mut_mat_en[, partIDs_in_mut] != "Silent")]
              }
              if (variant_class == "missense") {
                mut_partIDs <- partIDs_in_mut[grepl(x = mut_mat_en[, partIDs_in_mut], pattern = "Missense_Mutation")]
              }
              if (variant_class == "truncation") {
                mut_partIDs <- partIDs_in_mut[grepl(x = mut_mat_en[, partIDs_in_mut], pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del")]
              }
            } else {
              mut_partIDs <- NULL 
            }
            mut_silent_partIDs <- partIDs_in_mut[(mut_mat_en[, partIDs_in_mut] == "Silent")]
            
            # get the patient ids with TP53 deletion -------------------------------------
            cna_tab <- loadCNAstatus(cancer = cancer)
            # cna_tab <- loadTP53Deletion(cancer = cancer)
            
            cna_tab <- cna_tab[cna_tab$gene == gene_altered_tmp,]
            partIDs_in_cna <- colnames(cna_tab)[!(colnames(cna_tab) %in% c("gene"))]
            del_partIDs <- partIDs_in_cna[cna_tab[,partIDs_in_cna] == "deletion"]
            del_partIDs
            
            # get the patient ids of the controls -------------------------------------
            mut_cnv_partIDs <- unique(c(mut_partIDs, del_partIDs, mut_silent_partIDs))
            control_partIDs <- partIDs[!(partIDs %in% mut_cnv_partIDs)]
            
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
          write.table(x = mut_cnv_tab, file = fn, 
                      col.names = T, row.names = F, quote = F, sep = "\t")
          
        }
      }
    }
  }
  
}

stop("")
mut_cnv_cans <- NULL
for (gene_altered_tmp in genes_altered) {
  for (affected_exp_type in c("PRO", "PHO", "RNA")) {
    for (variant_class in c("truncation", "missense", "not_silent")) {
      for (cancer in c("CCRCC")) {
        for (genoalt_type in c("mut")) {
          fn <- paste0(resultD, "phospho_network/genoalt/tables/test_ccRCC_SMG_mut_impact_global/",  cancer, "_", variant_class, "_", gene_altered_tmp, "_mut_impact_", affected_exp_type , "_tab.txt")
          if (file.exists(fn)) {
            mut_cnv_can <- fread(input = fn, data.table = F, sep = "\t")
            # mut_cnv_can <- unique(mut_cnv_can)
            mut_cnv_can$variant_class <- variant_class
            mut_cnv_can$affected_exp_type <- affected_exp_type
            mut_cnv_cans <- rbind(mut_cnv_cans, mut_cnv_can)
          }
        }
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
mut_cnv_cans$fdr <- FDR_by_id_columns(p_vector = mut_cnv_cans$p, id_columns = c("cancer", "variant_class", "affected_exp_type", "SELF"), df = mut_cnv_cans)
mut_cnv_cans$fdr_by_gene <- FDR_by_id_columns(p_vector = mut_cnv_cans$p, id_columns = c("cancer", "variant_class", "affected_exp_type", "SELF", "GENE"), df = mut_cnv_cans)
 
write.table(x = mut_cnv_cans, file = paste0(makeOutDir(), paste0(genes_altered, collapse = "_"), "_mut_impact_proteome_RNA_cptac2p_cptac3_tab.txt"), sep = "\t", row.names = F, quote = F)
saveRDS(object = mut_cnv_cans, file = paste0(makeOutDir(), paste0(genes_altered, collapse = "_"), "_mut_impact_proteome_RNA_cptac2p_cptac3_tab.RDS"))
