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
cancers2process <- c("Brain", "Kidney", "Lung", "Pancreas", "Uterus_Endometrium")
# cancers2process <- c("CCRCC")

# expression data to process ----------------------------------------------
# affected_exp_type <- "PHO"
# affected_exp_type <- "PRO"
# affected_exp_type <- "RNA"
# affected_exp_types2process <- c("RNA", "PRO")
# affected_exp_types2process <- c("PRO")
affected_exp_types2process <- c("RNA")


# variant types to test ---------------------------------------------------
variant_classes2process <- c("hyper_vs_normal")


# input methylation status matrix -----------------------------------------


# input protein pair table ------------------------------------------------
pair_tab_annotated <- fread(input = "./cptac2p/analysis_results/dependencies/tables/compile_protein_pair_table/protein_pair_table.txt", data.table = F)


# Business  ---------------------------------------------------------------
for (affected_exp_type in affected_exp_types2process) {
  for (variant_class in variant_classes2process) {
    for (cancer in cancers2process) {
      fn2w <- paste0(makeOutDir(resultD = resultD), cancer, "_", variant_class,"_mut_impact_", affected_exp_type , "_tab.txt")
      if (!file.exists(fn2w)) {
        methyl_status <- fread(input = paste0("./Ding_Lab/Projects_Current/Methylation/resources/methylation_status/", cancer, ".sil.protein.10.meth.tumor.status.matched.csv"), data.table = F)
        
        partIDs_methyl <- colnames(methyl_status)[-1]
        colnames(methyl_status)[1] <- c("probe")
        
        methyl_status <- methyl_status %>%
          mutate(Gene = str_split_fixed(string = probe, pattern = "\\@", n = 3)[,3])
        methyl_status$num_methyl_samples <- rowSums(methyl_status[, partIDs_methyl], na.rm = T)
        methyl_status <- methyl_status %>%
          filter(num_methyl_samples >= num_genoalt_thres)
        
        methyl_status$Gene %>%
          unique() %>%
          length()
        
        ### filter by the methylated gene to test
        pair_tab <- merge(methyl_status %>%
                            select(probe, Gene),
                          pair_tab_annotated,
                          by.x = c("Gene"),
                          by.y = c("GENE"))
        
        if (affected_exp_type == "PRO") {
          fn <- paste0("./Ding_Lab/Projects_Current/Methylation/resources/methylation_status/", cancer, ".sil.protein.07.transcriptome.csv")
          affected_exp_data <- fread(input = fn, data.table = F)
          affected_exp_head <- data.frame(probe = affected_exp_data$V1)
          affected_exp_head <- affected_exp_head %>%
            mutate(SUBSTRATE = str_split_fixed(string = probe, pattern = "\\@", n = 3)[,3])
          affected_exp_head$SUB_MOD_RSD <- affected_exp_type
          
        }
        if (affected_exp_type == "RNA") {
          fn <- paste0("./Ding_Lab/Projects_Current/Methylation/resources/methylation_status/", cancer, ".sil.rna.07.transcriptome.csv")
          affected_exp_data <- fread(input = fn, data.table = F)
          affected_exp_data <- affected_exp_data[rowSums((affected_exp_data == 0 & !is.na(affected_exp_data)) | is.na(affected_exp_data)) < (as.numeric(ncol(affected_exp_data)) - num_genoalt_thres - num_control_thres -1),]
          affected_exp_head <- data.frame(probe = affected_exp_data$V1)
          affected_exp_head <- affected_exp_head %>%
            mutate(SUBSTRATE = str_split_fixed(string = probe, pattern = "\\@", n = 3)[,3])
          affected_exp_head$SUB_MOD_RSD <- affected_exp_type
        }
        if (affected_exp_type == "PHO") {
          stop("meow")
        }
        
        ## get the overlap of the participant IDs
        partIDs <- colnames(affected_exp_data)
        partIDs_overlap <- intersect(partIDs, colnames(methyl_status))
        print(length(partIDs_overlap))

        if (affected_exp_type %in% c("PRO", "RNA")) {
          pair_tab_sites <- merge(pair_tab, affected_exp_head %>%
                                    select(SUBSTRATE, SUB_MOD_RSD) %>%
                                    unique(), 
                                  by.x = c("SUB_GENE"), by.y = c("SUBSTRATE"))
        } else {
          pair_tab_sites <- merge(pair_tab, affected_exp_head[, c("SUBSTRATE", "SUB_MOD_RSD")], by.x = c("SUB_GENE"), by.y = c("SUBSTRATE"))
        }
        # pair_tab_sites <- pair_tab_sites[!is.na(pair_tab_sites$GENE) & !is.na(pair_tab_sites$SUB_MOD_RSD),]
        # pair_tab_sites <- pair_tab_sites[order(pair_tab_sites$GENE, pair_tab_sites$SUB_GENE),]
        # print(nrow(pair_tab_sites))
        
        ## initiate value columns
        num_control <- vector(mode = "numeric", length = nrow(pair_tab_sites))
        num_mut <- num_control; meddiff_mut <- num_control; p_mut <- num_control; meddiff_bottom_mut <- num_control; meddiff_upper_mut <- num_control
        
        for (altered_probe in unique(pair_tab_sites$probe)) {
          methyl_status_altered_probe <- methyl_status %>%
            filter(probe == altered_probe)
            
          if (variant_class == "hyper_vs_normal") {
            altered_partIDs <- partIDs_overlap[(methyl_status_altered_probe[,partIDs_overlap] == 1) & !is.na(methyl_status_altered_probe[, partIDs_overlap])]
          }
          control_partIDs <- partIDs_overlap[(methyl_status_altered_probe[,partIDs_overlap] == 0) & !is.na(methyl_status_altered_probe[, partIDs_overlap])]
          
          print(paste0("altered_probe:", altered_probe))
          for (substrate in unique(pair_tab_sites$SUB_GENE[pair_tab_sites$probe == altered_probe])) {
            print(paste0("substrate:", substrate))
            
            for (site in unique(pair_tab_sites$SUB_MOD_RSD[pair_tab_sites$probe == altered_probe & pair_tab_sites$SUB_GENE == substrate])) {
              print(paste0("site:", site))
              
              if (affected_exp_type  %in% c("PRO", "RNA")) {
                affected_exp_altered <- head(affected_exp_data[affected_exp_head$SUBSTRATE == substrate,], n = 1)
                if (nrow(affected_exp_altered) == 0) {
                  next()
                }
                i <- (pair_tab_sites$probe == altered_probe & pair_tab_sites$SUB_GENE == substrate & pair_tab_sites$SUB_MOD_RSD == site)
              } else {
                affected_exp_altered <- affected_exp_data[affected_exp_head$SUBSTRATE == substrate & affected_exp_head$SUB_MOD_RSD == site,]
                i <- (pair_tab_sites$probe == altered_probe & pair_tab_sites$SUB_GENE == substrate & pair_tab_sites$SUB_MOD_RSD == site)
              }
              
              affected_exp_control <- affected_exp_altered[, control_partIDs]
              affected_exp_control <- affected_exp_control[!is.na(affected_exp_control)]
              num_control[i] <- length(affected_exp_control)
              # stop("meow")
              if (num_control[i] >= num_control_thres) {
                affected_exp_methyl <- affected_exp_altered[, altered_partIDs]
                affected_exp_methyl <- affected_exp_methyl[!is.na(affected_exp_methyl)]
                num_mut[i] <- length(affected_exp_methyl)
                
                if (num_mut[i] > 0) {
                  meddiff_mut[i] <- median(affected_exp_methyl) - median(affected_exp_control)
                  if (num_mut[i] >= num_genoalt_thres){
                    stat <- wilcox.test(x = affected_exp_methyl, y = affected_exp_control, conf.int = T)
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
        mut_cnv_tab$genoalt_type <- "methylation"
        
        ## clean up
        write.table(x = mut_cnv_tab, file = fn2w, col.names = T, row.names = F, quote = F, sep = "\t")
      }
    }
  }
}

# mut_cnv_cans <- NULL
# for (affected_exp_type in c("PRO", "RNA")) {
#   for (variant_class in c("hyper_vs_normal")) {
#     for (cancer in cancers2process) {
#       for (genoalt_type in c("methylation")) {
#         mut_cnv_can <- fread(input = paste0(makeOutDir(resultD = resultD), cancer, "_", variant_class, "_mut_impact_", affected_exp_type , "_tab.txt"), data.table = F, sep = "\t")
#         mut_cnv_can$variant_class <- variant_class
#         mut_cnv_can$affected_exp_type <- affected_exp_type
#         mut_cnv_can$cancer <- cancer
#         mut_cnv_cans <- rbind(mut_cnv_cans, mut_cnv_can)
#       }
#     }
#   }
# }
# if (length(unique(mut_cnv_cans$cancer)) < length(cancers2process)) {
#   stop("cancer type not enough!")
# }
# if (length(unique(mut_cnv_cans$variant_class)) < 3) {
#   stop("variant_class not enough!")
# }
# if (length(unique(mut_cnv_cans$affected_exp_type)) < 3) {
#   stop("affected_exp_type not enough!")
# }
# mut_cnv_cans <- mut_cnv_cans[mut_cnv_cans$num >= num_genoalt_thres,]
# mut_cnv_cans <- merge(mut_cnv_cans, pair_tab_annotated, by = c( "GENE", "SUB_GENE"), all.x = T)
# mut_cnv_cans$fdr <- FDR_by_id_columns(p_vector = mut_cnv_cans$p, id_columns = c("SELF", "cancer", "variant_class", "affected_exp_type", colnames(pair_tab_annotated)[!(colnames(pair_tab_annotated) %in% c("GENE", "SUB_GENE", "pair_pro"))]), df = mut_cnv_cans)
# mut_cnv_cans$fdr_by_gene <- FDR_by_id_columns(p_vector = mut_cnv_cans$p, id_columns = c("GENE" , "SELF", "cancer", "variant_class", "affected_exp_type", colnames(pair_tab_annotated)[!(colnames(pair_tab_annotated) %in% c("GENE", "SUB_GENE", "pair_pro"))]), df = mut_cnv_cans)
# write.table(x = mut_cnv_cans, file = paste0(makeOutDir(resultD = resultD), "SMG_mut_impact_tab.txt"), sep = "\t", row.names = F, quote = F)
# 
