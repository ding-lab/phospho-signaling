# Yige Wu @ WashU 2019 Sep
## test copy number impact of copy number altered genes on protein/phosphorylation within kinase-substrate pairs or protein complex pairs for GBM
## using wilcox test
## reference: https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
sample_type <- "tumor"
id_vars <- c("GENE", "SUB_GENE", "SUB_MOD_RSD", "cancer")
prefix_vars <- c("p", "num", "meddiff_bottom", "meddiff_upper", "meddiff")
version_tmp <- 1
num_control_thres <- 4
num_genoalt_thres <- 4

# cancer types to process -------------------------------------------------
cancers2process <- c("GBM")

# expression data to process ----------------------------------------------
affected_exp_types2process <- c("PRO", "PHO")


# variant types to test ---------------------------------------------------
variant_classes2process <- c("gain", "loss")

# input protein pair table ------------------------------------------------
pair_tab_annotated <- fread(input = "./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/analysis_results/dependencies/tables/compile_protein_pair_table/protein_pair_table_v2.txt", data.table = F)

# Business  ---------------------------------------------------------------
for (affected_exp_type in affected_exp_types2process) {
  for (variant_class in variant_classes2process) {
    for (cancer in cancers2process) {
      pair_tab <- pair_tab_annotated[pair_tab_annotated$GENE %in% unlist(CNGs[[cancer]]), c("GENE", "SUB_GENE")]
      
      fn <- paste0(makeOutDir(), cancer, "_", variant_class,"_mut_impact_", affected_exp_type , ".", format(Sys.Date(), "%Y%m%d") , ".v", version_tmp , ".txt")
      
      if (!file.exists(fn)) {
        pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "d4")
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "d6")
        
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
        mut_mat <- generate_somatic_mutation_matrix(pair_tab = pair_tab, maf = maf)
        
        ## input cna matrix
        cna_tab <- loadCNAstatus(cancer = cancer)
        cna_tab <- cna_tab[!duplicated(cna_tab$gene),]
        rownames(cna_tab) <- as.vector(cna_tab$gene)
        
        ## get the overlap of the participant IDs
        partIDs <- colnames(affected_exp_data); partIDs <- partIDs[!(partIDs %in% c("Hugo_Symbol", "Gene", "gene", "Phosphosite", "Peptide_ID"))]
        partIDs_in_mut <- colnames(mut_mat); partIDs_in_mut <- partIDs_in_mut[!(partIDs_in_mut %in% c("Hugo_Symbol", "Gene", "gene", "Phosphosite", "Peptide_ID"))]
        partIDs_in_cna <- colnames(cna_tab); partIDs_in_cna <- partIDs_in_cna[!(partIDs_in_cna %in% c("Hugo_Symbol", "Gene", "gene", "Phosphosite", "Peptide_ID"))]
        
        ## get the count of smaples altered in each gene so as to narrow down to particular genes to test
        altered_count <- rowSums(x = cna_tab[, partIDs_in_cna] == variant_class); 
        
        ## get the genes to test
        genes_num_genoalt_thresholded <- names(altered_count)[altered_count >= num_genoalt_thres]
        print(length(genes_num_genoalt_thresholded))
        
        ## initiate identifier columns
        pair_tab_genoalt_thresholded <- pair_tab[pair_tab$GENE %in% genes_num_genoalt_thresholded,]
        print(nrow(pair_tab_genoalt_thresholded))
        
        if (affected_exp_type %in% c("PRO", "RNA")) {
          pair_tab_sites <- pair_tab_genoalt_thresholded
          pair_tab_sites$SUB_MOD_RSD <- affected_exp_type
        } else {
          pair_tab_sites <- merge(pair_tab_genoalt_thresholded, affected_exp_head[, c("SUBSTRATE", "SUB_MOD_RSD")], by.x = c("SUB_GENE"), by.y = c("SUBSTRATE"), all.y = T)
        }
        pair_tab_sites <- pair_tab_sites[!is.na(pair_tab_sites$GENE) & !is.na(pair_tab_sites$SUB_MOD_RSD),]
        pair_tab_sites <- pair_tab_sites[order(pair_tab_sites$GENE, pair_tab_sites$SUB_GENE),]
        print(nrow(pair_tab_sites))
        
        ## initiate value columns
        num_control <- vector(mode = "numeric", length = nrow(pair_tab_sites))
        num_mut <- num_control; meddiff_mut <- num_control; p_mut <- num_control; meddiff_bottom_mut <- num_control; meddiff_upper_mut <- num_control
        
        for (enzyme in unique(pair_tab_sites$GENE)) {
          # get the patient ids with somatic mutation -------------------------------
          mut_mat_en <- mut_mat[mut_mat$Hugo_Symbol == enzyme,]
          
          if (nrow(mut_mat_en) > 0){
            mut_partIDs <- partIDs_in_mut[(mut_mat_en[,partIDs_in_mut] != "") & (mut_mat_en[, partIDs_in_mut] != "Silent")]
          } else {
            mut_partIDs <- NULL 
          }
          
          # get the patient ids with cna -------------------------------
          cna_geneA <- cna_tab[cna_tab$gene == enzyme,]
          altered_partIDs <- partIDs_in_cna[cna_geneA[,partIDs_in_cna] == variant_class]; altered_partIDs <- intersect(partIDs, altered_partIDs)
          cna_partIDs <- partIDs_in_cna[cna_geneA[,partIDs_in_cna] != "neutral"]
          
          mut_cnv_partIDs <- unique(c(mut_partIDs, cna_partIDs))
          control_partIDs <- partIDs[!(partIDs %in% mut_cnv_partIDs)]
          
          print(paste0("enzyme:", enzyme))
          for (substrate in unique(pair_tab_sites$SUB_GENE[pair_tab_sites$GENE == enzyme])) {
            print(paste0("substrate:", substrate))
            
            for (site in unique(pair_tab_sites$SUB_MOD_RSD[pair_tab_sites$GENE == enzyme & pair_tab_sites$SUB_GENE == substrate])) {
              print(paste0("site:", site))
              
              if (affected_exp_type  %in% c("PRO", "RNA")) {
                affected_exp_geneB <- affected_exp_data[affected_exp_head$SUBSTRATE == substrate,]
                i <- (pair_tab_sites$GENE == enzyme & pair_tab_sites$SUB_GENE == substrate & pair_tab_sites$SUB_MOD_RSD == site)
              } else {
                affected_exp_geneB <- affected_exp_data[affected_exp_head$SUBSTRATE == substrate & affected_exp_head$SUB_MOD_RSD == site,]
                i <- (pair_tab_sites$GENE == enzyme & pair_tab_sites$SUB_GENE == substrate & pair_tab_sites$SUB_MOD_RSD == site)
                
              }
              affected_exp_geneB <- affected_exp_geneB[1,]
              affected_exp_control <- affected_exp_geneB[, control_partIDs]
              affected_exp_control <- affected_exp_control[!is.na(affected_exp_control)]
              num_control[i] <- length(affected_exp_control)
              
              if (num_control[i] >= num_control_thres) {
                affected_exp_altered <- affected_exp_geneB[, altered_partIDs]
                affected_exp_altered <- affected_exp_altered[!is.na(affected_exp_altered)]
                num_mut[i] <- length(affected_exp_altered)
                
                if (num_mut[i] > 0) {
                  meddiff_mut[i] <- median(affected_exp_altered) - median(affected_exp_control)
                  if (num_mut[i] >= num_genoalt_thres){
                    stat <- wilcox.test(x = affected_exp_altered, y = affected_exp_control, conf.int = T)
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
        mut_cnv_tab$pair_pro <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE)
        mut_cnv_tab$pair <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE, ":", mut_cnv_tab$SUB_MOD_RSD)
        ## annotate cis and trans
        mut_cnv_tab$SELF <- "trans"
        mut_cnv_tab$SELF[as.vector(mut_cnv_tab$GENE) == as.vector(mut_cnv_tab$SUB_GENE)] <- "cis"
        mut_cnv_tab$genoalt_type <- "cna"
        mut_cnv_tab$variant_class <- variant_class
        mut_cnv_tab$affected_exp_type <- affected_exp_type
        mut_cnv_tab <- mut_cnv_tab[mut_cnv_tab$num >= num_genoalt_thres,]
        mut_cnv_tab <- merge(mut_cnv_tab, pair_tab_annotated, by = c( "GENE", "SUB_GENE", "pair_pro"), all.x = T)
        
        mut_cnv_tab$fdr <- FDR_by_id_columns(p_vector = mut_cnv_tab$p, id_columns = c("SELF", "cancer", "variant_class", "affected_exp_type", colnames(pair_tab_annotated)[!(colnames(pair_tab_annotated) %in% c("GENE", "SUB_GENE", "pair_pro"))]), df = mut_cnv_tab)
        mut_cnv_tab$fdr_by_gene <- FDR_by_id_columns(p_vector = mut_cnv_tab$p, id_columns = c("GENE" , "SELF", "cancer", "variant_class", "affected_exp_type", colnames(pair_tab_annotated)[!(colnames(pair_tab_annotated) %in% c("GENE", "SUB_GENE", "pair_pro"))]), df = mut_cnv_tab)
        write.table(x = mut_cnv_tab, file = fn, col.names = T, row.names = F, quote = F, sep = "\t")
        
      }
    }
  }
}

# mut_cnv_cans <- NULL
# for (affected_exp_type in c("PRO", "PHO", "RNA")) {
#   for (variant_class in c("truncation", "missense", "not_silent")) {
#     for (cancer in c("BRCA", "OV", "CO", "UCEC", "CCRCC", "LIHC")) {
#       for (genoalt_type in c("mut")) {
#         mut_cnv_can <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_impact_proteome/", cancer, "_", variant_class, "_mut_impact_", affected_exp_type , "_tab.txt"), data.table = F, sep = "\t")
#         stop()
#         mut_cnv_can$variant_class <- variant_class
#         mut_cnv_cans <- rbind(mut_cnv_cans, mut_cnv_can)
#       }
#     }
#   }
# }
# mut_cnv_cans <- mut_cnv_cans[mut_cnv_cans$num >= num_genoalt_thres,]
# write.table(x = mut_cnv_cans, file = paste0(makeOutDir(), "mut_impact_proteome_RNA_cptac2p_cptac3_tab.txt"), sep = "\t", row.names = F, quote = F)
# 
