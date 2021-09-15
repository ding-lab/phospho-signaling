# Yige Wu @ WashU 2021 Mar
## test mutation impact of SMGs on protein within kinase-gene_affected pairs or protein complex pairs 
## using wilcox test
## reference: https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/"
setwd(dir_base)
source("./phospho-signaling_analysis/load_pkgs.R")
source("./phospho-signaling_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# set variables -----------------------------------------------------------
sample_type <- "tumor"
id_vars <- c("GENE", "SUB_GENE", "SUB_MOD_RSD", "cancer")
prefix_vars <- c("p", "num", "meddiff_bottom", "meddiff_upper", "meddiff")
version_tmp <- 1
num_control_thres <- 5
num_genoalt_thres <- 5
cancer_process <- "ccRCC"
affected_exp_types2process <- c("PRO", "RNA")
# affected_exp_types2process <- c("PRO")
# affected_exp_types2process <- c("RNA")
variant_classes2process <- c("not_silent", "missense", "truncation")
# variant_classes2process <- c("not_silent")
## input SMGs to test
genes_mut_test <- c("VHL", "PBRM1", "BAP1", "TCEB1", "SETD2", "TP53", "PIK3CA", "MTOR", "KDM5C", "PTEN","TSC1",
                "CCNB2",
                "MSR1", "TXNIP", "BTNL3", "SLITRK6", "RHEB", "ARID1A", "NPNT", 
                "FPGT", "MUDENG", "TET2", "MUC4", "MLLT10", "MSGN1", "KRT32", "M6PR", "RPL14", "GRB7", "CSMD3", "DNHD1", "NLRP12", "VMO1", "OR4C13", "KCNMA1", "LMAN2L", "ZNF536", "YIPF3")

# input dependencies ------------------------------------------------------
## input the pair table
genepair_df <- fread(input = "./analysis_results/dependencies/tables/compile_protein_pair_table/protein_pair_table_v2.txt", data.table = F)
## input the mutation data
maf_df <- fread(data.table = F, input = "~/Box/Data_Freeze_1.0/CPTAC_ccRCC_Combined/Somatic_mutation/CPTAC_ccRCC_combined_somatic_mutation.maf_v1.0.tsv")
## input protein data
pro_df <- fread(data.table = F, input = "~/Box/Data_Freeze_1.0/CPTAC_ccRCC_Combined/Proteome/UMich_pipeline/ccRCC_DIA_Combined_norm_03042021.tsv", fill = T)
## input RNA data
rna_df <- fread(data.table = F, input = "~/Box/Data_Freeze_1.0/CPTAC_ccRCC_Combined/Gene_expression/CPTAC_ccRCC_combined_tumor_mRNA_FPKM_UQ_log2_v1.0.tsv")

# reprocess ---------------------------------------------------------------
## filter the pairs by SMGs
genepair_test_df <- genepair_df[genepair_df$GENE %in% genes_mut_test, c("GENE", "SUB_GENE")]
## geenrate mutation matrix to help devide mutant vs non-mutant
maf_df <- maf_df %>%
  mutate(SampleID = gsub(x = Tumor_Sample_Barcode, pattern = "_T", replacement = ""))
mut_mat <- get_mutation_class_matrix(pair_tab = genepair_test_df$GENE, maf = maf_df)
number_mutated <- rowSums(mut_mat[,-1] != "")
# genes_mut_test <- names(number_mutated[number_mutated >= num_genoalt_thres])
# View(mut_mat)

# Business  ---------------------------------------------------------------
result_df <- NULL
for (affected_exp_type in affected_exp_types2process) {
  if (affected_exp_type == "PRO") {
    affected_exp_head <- data.frame(gene_affected = pro_df$Genes)
    
    colnames_pro <- colnames(pro_df)
    affected_exp_data <- pro_df[, c(colnames_pro[grepl(pattern = "_T", x = colnames_pro)])]
    exp_sampleids <- gsub(x = colnames(affected_exp_data), pattern = "_T", replacement = "")
    exp_sampleids <- gsub(x = exp_sampleids, pattern = "\\.1", replacement = "")
    colnames(affected_exp_data) <- exp_sampleids
    affected_exp_data <- log2(affected_exp_data)
  }
  if (affected_exp_type == "RNA") {
    affected_exp_data <- rna_df
    affected_exp_head <- data.frame(gene_affected = rna_df$gene_name)
    exp_sampleids <- colnames(affected_exp_data)
    exp_sampleids <- exp_sampleids[!(exp_sampleids %in% c("gene_id", "gene_name"))]
    exp_sampleids <- gsub(pattern = "\\-T", replacement = "", x = exp_sampleids)
    colnames(affected_exp_data) <- c("gene_id", "gene_name", exp_sampleids)
  }

  for (variant_class in variant_classes2process) {
    fn <- paste0(dir_out, cancer_process, "_", variant_class,"_mut_impact_", affected_exp_type , ".tsv")
    if (!file.exists(fn)) {
      
      ## get the overlap of the participant IDs
      sampleids_overlap <- intersect(exp_sampleids, colnames(mut_mat))
      print(length(sampleids_overlap))
      
      ## the number of mutations
      if (variant_class == "not_silent") {
        mut_count <- rowSums(x = ((mut_mat[, sampleids_overlap] != "") & (mut_mat[, sampleids_overlap] != "Silent"))); 
      }
      if (variant_class == "missense") {
        mut_count <- sapply(1:nrow(mut_mat), FUN = function(i, mut_mat) sum(grepl(x = mut_mat[i,], pattern = "Missense", ignore.case = T)), mut_mat[,sampleids_overlap])
      }
      if (variant_class == "truncation") {
        mut_count <- sapply(1:nrow(mut_mat), FUN = function(i, mut_mat) sum(grepl(x = mut_mat[i,], pattern = "Nonsense|Frame_Shift_Ins|Frame_Shift_Del|Splice", ignore.case = T)), mut_mat[,sampleids_overlap])
      }
      names(mut_count) <- mut_mat$Hugo_Symbol;mut_count
      control_count <- length(exp_sampleids) - mut_count; control_count

      ## get the genes to test
      genes_num_genoalt_thresholded <- names(mut_count)[mut_count >= num_genoalt_thres & control_count >= num_control_thres]
      print(length(genes_num_genoalt_thresholded))
      
      ## initiate identifier columns
      genepair_test_df_genoalt_thresholded <- genepair_test_df[genepair_test_df$GENE %in% genes_num_genoalt_thresholded,]
      print(nrow(genepair_test_df_genoalt_thresholded))
      
      if (affected_exp_type %in% c("PRO", "RNA")) {
        genepair_test_df_sites <- genepair_test_df_genoalt_thresholded
        genepair_test_df_sites$SUB_MOD_RSD <- affected_exp_type
      } else {
        genepair_test_df_sites <- merge(genepair_test_df_genoalt_thresholded, affected_exp_head[, c("gene_affected", "SUB_MOD_RSD")], by.x = c("SUB_GENE"), by.y = c("gene_affected"))
      }
      genepair_test_df_sites <- genepair_test_df_sites[!is.na(genepair_test_df_sites$GENE) & !is.na(genepair_test_df_sites$SUB_MOD_RSD),]
      genepair_test_df_sites <- genepair_test_df_sites[order(genepair_test_df_sites$GENE, genepair_test_df_sites$SUB_GENE),]
      print(nrow(genepair_test_df_sites))
      
      ## initiate value columns
      num_control <- vector(mode = "numeric", length = nrow(genepair_test_df_sites))
      num_mut <- num_control; meddiff_mut <- num_control; p_mut <- num_control; meddiff_bottom_mut <- num_control; meddiff_upper_mut <- num_control
      
      for (gene_mut in unique(genepair_test_df_sites$GENE)) {
        mut_mat_en <- mut_mat[mut_mat$Hugo_Symbol == gene_mut,]
        if (nrow(mut_mat_en) > 0){
          if (variant_class == "not_silent") {
            mut_sampleids <- sampleids_overlap[(mut_mat_en[,sampleids_overlap] != "") & (mut_mat_en[, sampleids_overlap] != "Silent")]
          }
          if (variant_class == "missense") {
            mut_sampleids <- sampleids_overlap[grepl(x = mut_mat_en[, sampleids_overlap], pattern = "Missense_Mutation")]
          }
          if (variant_class == "truncation") {
            mut_sampleids <- sampleids_overlap[grepl(x = mut_mat_en[, sampleids_overlap], pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del|Splice_Site")]
          }
          
        } else {
          mut_sampleids <- NULL 
        }
        
        mut_cnv_sampleids <- mut_sampleids
        control_sampleids <- exp_sampleids[!(exp_sampleids %in% mut_sampleids)]
        
        print(paste0("gene_mut:", gene_mut))
        for (gene_affected in unique(genepair_test_df_sites$SUB_GENE[genepair_test_df_sites$GENE == gene_mut])) {
          print(paste0("gene_affected:", gene_affected))
          
          for (site in unique(genepair_test_df_sites$SUB_MOD_RSD[genepair_test_df_sites$GENE == gene_mut & genepair_test_df_sites$SUB_GENE == gene_affected])) {
            print(paste0("site:", site))
            
            if (affected_exp_type  %in% c("PRO", "RNA")) {
              affected_gene_exp_data <- affected_exp_data[affected_exp_head$gene_affected == gene_affected,]
              i <- (genepair_test_df_sites$GENE == gene_mut & genepair_test_df_sites$SUB_GENE == gene_affected & genepair_test_df_sites$SUB_MOD_RSD == site)
            } else {
              affected_gene_exp_data <- affected_exp_data[affected_exp_head$gene_affected == gene_affected & affected_exp_head$SUB_MOD_RSD == site,]
              i <- (genepair_test_df_sites$GENE == gene_mut & genepair_test_df_sites$SUB_GENE == gene_affected & genepair_test_df_sites$SUB_MOD_RSD == site)
              
            }
            affected_gene_exp_data <- affected_gene_exp_data[1,]
            affected_exp_control <- affected_gene_exp_data[, control_sampleids]
            affected_exp_control <- affected_exp_control[!is.na(affected_exp_control)]
            num_control[i] <- length(affected_exp_control)
            
            if (num_control[i] >= num_control_thres) {
              affected_exp_mut <- affected_gene_exp_data[, mut_sampleids]
              affected_exp_mut <- affected_exp_mut[!is.na(affected_exp_mut)]
              num_mut[i] <- length(affected_exp_mut)
              if (num_mut[i] > 0) {
                meddiff_mut[i] <- median(affected_exp_mut) - median(affected_exp_control)
                if (num_mut[i] >= num_genoalt_thres){
                  # stat <- wilcox.test(x = affected_exp_mut, y = affected_exp_control, conf.int = T)
                  stat <- tryCatch(expr = wilcox.test(x = affected_exp_mut, y = affected_exp_control, conf.int = T),
                                           error = function(e) {warning("wilcox.test failed.");return(NULL)})
                  if (!is.null(stat)) {
                    p_mut[i] <- stat$p.value
                    meddiff_bottom_mut[i] <- stat$conf.int[1]
                    meddiff_upper_mut[i] <- stat$conf.int[2]
                  } else {
                    p_mut[i] <- NA
                    meddiff_bottom_mut[i] <- NA
                    meddiff_upper_mut[i] <- NA
                  }
                }
              }
            }
            print(which(i))
          }
        }
        
      }
      mut_cnv_tab <- genepair_test_df_sites
      mut_cnv_tab <- cbind(mut_cnv_tab, 
                           data.frame(num = num_mut, 
                                      num_control = num_control,
                                      p = p_mut, 
                                      meddiff_bottom = meddiff_bottom_mut, 
                                      meddiff_upper = meddiff_upper_mut, 
                                      meddiff = meddiff_mut))
      mut_cnv_tab$cancer <- cancer_process
      mut_cnv_tab$pair_pro <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE)
      mut_cnv_tab$pair <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE, ":", mut_cnv_tab$SUB_MOD_RSD)
      ## annotate cis and trans
      mut_cnv_tab$SELF <- "trans"
      mut_cnv_tab$SELF[as.vector(mut_cnv_tab$GENE) == as.vector(mut_cnv_tab$SUB_GENE)] <- "cis"
      mut_cnv_tab$genoalt_type <- "mut"
      mut_cnv_tab$variant_class <- variant_class
      mut_cnv_tab$affected_exp_type <- affected_exp_type
      mut_cnv_tab <- mut_cnv_tab[mut_cnv_tab$num >= num_genoalt_thres,]
      mut_cnv_tab <- merge(mut_cnv_tab, genepair_df, by = c( "GENE", "SUB_GENE", "pair_pro"), all.x = T)
      
      mut_cnv_tab$fdr <- FDR_by_id_columns(p_vector = mut_cnv_tab$p, id_columns = c("SELF", "cancer", "variant_class", "affected_exp_type", colnames(genepair_df)[!(colnames(genepair_df) %in% c("GENE", "SUB_GENE", "pair_pro"))]), df = mut_cnv_tab)
      mut_cnv_tab$fdr_by_gene <- FDR_by_id_columns(p_vector = mut_cnv_tab$p, id_columns = c("GENE" , "SELF", "cancer", "variant_class", "affected_exp_type", colnames(genepair_df)[!(colnames(genepair_df) %in% c("GENE", "SUB_GENE", "pair_pro"))]), df = mut_cnv_tab)
      write.table(x = mut_cnv_tab, file = fn, col.names = T, row.names = F, quote = F, sep = "\t")
      
      ## store result
      result_df <- rbind(result_df, mut_cnv_tab)
    }
  }
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, cancer_process, "_mut_impact_RNA_PRO.tsv")
write.table(x = result_df, file = file2write, quote = F, sep = "\t", row.names = F)

