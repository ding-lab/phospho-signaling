# Yige Wu @ WashU 2018 Aug
## test mutation impact on protein/phosphorylation within kinase-substrate pairs or protein complex pairs
## using wilcox test
## reference: https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/
## take out the samples deleted the altered gene

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
setwd(baseD)
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/p53/TP53_shared.R")

# set variables -----------------------------------------------------------
sample_type <- "tumor"
id_vars <- c("GENE", "SUB_GENE", "SUB_MOD_RSD", "cancer")
prefix_vars <- c("p", "num", "meddiff_bottom", "meddiff_upper", "meddiff")
num_control_thres <- 4
num_genoalt_thres <- 4
gene_altered <- "TP53"
# cancer types to process -------------------------------------------------
# cancers2process <- c("BRCA", "OV", "CO")
# cancers2process <- c("CCRCC")
cancers2process <- c("BRCA", "CO", "OV", "UCEC", "CCRCC", "LIHC")
# cancers2process <- c("BRCA")
# cancers2process <- c("CO")

# expression data to process ----------------------------------------------
# affected_exp_type <- "PHO"
# affected_exp_type <- "PRO"
# affected_exp_type <- "RNA"
affected_exp_types2process <- c("RNA", "PRO", "PHO")
# affected_exp_types2process <- c("PRO", "PHO")
# affected_exp_types2process <- c("PHO")

# variant types to test ---------------------------------------------------
## include all nonsynonymous mutations or just missense/truncation
# variant_class <- "truncation"
# variant_class <- "not_silent"
# variant_class <- "missense"
# variant_classes2process <- c("truncation", "missense")
# variant_classes2process <- c("not_silent")
# variant_classes2process <- c("missense", "truncation", "not_silent")
variant_classes2process <- c("Missense_protein_level_high_vs_Missense_protein_level_average")
variant_classes2process <- c("Missense_protein_level_high_vs_Truncation")

# input TF relations ------------------------------------------------------
tf_tab <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/PPI/TF_interactions.txt", data.table = F)
tf_tab %>%
  nrow()

tf_tab %>%
  select(tfregulons_level) %>%
  table()

tf_tab %>%
  select(sources) %>%
  table()

# Input TP53 interactor table ---------------------------------------------
TP53_pair_tab <- fread("./cptac2p/analysis_results/p53/tables/compile_TP53_interactor_pair_tab/TP53_interactor_pair_tab.txt", data.table = F)
TP53_pair_tab2add <- data.frame(GENE = "TP53", SUB_GENE = GOF_interacting_TFs)
pair_tab_annotated <- fread(input = "./cptac2p/analysis_results/dependencies/tables/compile_protein_pair_table/protein_pair_table.txt", data.table = F)

# get all the TF downstreams of TP53 interactors -------------------------------------------------
tf_tab2add <- tf_tab %>%
  filter(source_genesymbol %in% TP53_pair_tab$SUB_GENE)
tf_targets2add <- unique(tf_tab2add$target_genesymbol)

pair_tab <- data.frame(GENE = "TP53", SUB_GENE = tf_targets2add)

# Business  ---------------------------------------------------------------
for (affected_exp_type in affected_exp_types2process) {
  for (variant_class in variant_classes2process) {
    for (cancer in cancers2process) {
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
      
      ## input patient info table
      part_tab <- fread(input = paste0("./cptac2p/analysis_results/p53/tables/generate_patient_TP53_status/", cancer, "_tp53_mut_protein_classification_table.txt"), data.table = F)
      
      ## get the overlap of the participant IDs
      if (affected_exp_type %in% c("PRO", "RNA")) {
        pair_tab_sites <- pair_tab
        pair_tab_sites$SUB_MOD_RSD <- affected_exp_type
      } else {
        pair_tab_sites <- merge(pair_tab, affected_exp_head[, c("SUBSTRATE", "SUB_MOD_RSD")], by.x = c("SUB_GENE"), by.y = c("SUBSTRATE"), all.y = T)
      }
      pair_tab_sites <- pair_tab_sites[!is.na(pair_tab_sites$GENE) & !is.na(pair_tab_sites$SUB_MOD_RSD),]
      pair_tab_sites <- pair_tab_sites[order(pair_tab_sites$GENE, pair_tab_sites$SUB_GENE),]
      print(nrow(pair_tab_sites))
      
      ## initiate value columns
      num_control <- vector(mode = "numeric", length = nrow(pair_tab_sites))
      num_mut <- num_control; meddiff_mut <- num_control; p_mut <- num_control; meddiff_bottom_mut <- num_control; meddiff_upper_mut <- num_control
      
      # get the patient ids with mutations -------------------------------------
      
      for (enzyme in unique(pair_tab_sites$GENE)) {
        if (variant_class == "Missense_protein_level_high_vs_Missense_protein_level_average") {
          mut_partIDs <- part_tab$V1[part_tab$Classification_complex == "Missense_protein_level_high"]
          control_partIDs <- part_tab$V1[part_tab$Classification_complex == "Missense_protein_level_average"]
        }
        if (variant_class == "Missense_protein_level_high_vs_Truncation") {
          mut_partIDs <- part_tab$V1[part_tab$Classification_complex == "Missense_protein_level_high"]
          control_partIDs <- part_tab$V1[part_tab$Classification_complex == "Truncation"]
        }
        
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
            affected_exp_control <- affected_exp_altered[, intersect(colnames(affected_exp_altered), control_partIDs)]
            affected_exp_control <- affected_exp_control[!is.na(affected_exp_control)]
            if (affected_exp_type == "RNA") {
              affected_exp_control <- affected_exp_control[affected_exp_control != 0]
            }
            num_control[i] <- length(affected_exp_control)
            
            if (num_control[i] >= num_control_thres) {
              mut_pho <- affected_exp_altered[, mut_partIDs]
              mut_pho <- mut_pho[!is.na(mut_pho)]
              
              if (affected_exp_type == "RNA") {
                mut_pho <- mut_pho[mut_pho != 0]
              }
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
      
      write.table(x = mut_cnv_tab, file = paste0(makeOutDir(resultD = resultD), cancer, "_", variant_class, "_", gene_altered, "_mut_impact_", affected_exp_type , "_tab.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
      
    }
  }
}
# file.exists('./cptac2p/analysis_results/p53/tables/test_TP53_missense_impact_TF_downstream/BRCA_Missense_protein_level_high_vs_Missense_protein_level_average_mut_impact_PRO_tab.txt')
# file.exists('./cptac2p/analysis_results/p53/tables/test_TP53_missense_impact_TF_downstream/BRCA_Missense_protein_level_high_vs_Missense_protein_level_average_TP53_mut_impact_PRO_tab.txt')
# stop("meow")
mut_cnv_cans <- NULL
for (affected_exp_type in c("PRO", "PHO", "RNA")) {
  for (variant_class in c("Missense_protein_level_high_vs_Missense_protein_level_average",
                          "Missense_protein_level_high_vs_Truncation")) {
    for (cancer in c("BRCA", "OV", "CO", "UCEC", "CCRCC", "LIHC")) {
      for (genoalt_type in c("mut")) {
        fn <- paste0(makeOutDir(resultD = resultD), cancer, "_", variant_class, "_", gene_altered, "_mut_impact_", affected_exp_type , "_tab.txt")
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
# if (length(unique(mut_cnv_cans$variant_class)) < 3) {
#   stop("variant_class not enough!")
# }
if (length(unique(mut_cnv_cans$affected_exp_type)) < 3) {
  stop("affected_exp_type not enough!")
}
mut_cnv_cans <- mut_cnv_cans[mut_cnv_cans$num >= num_genoalt_thres,]
mut_cnv_cans <- unique(mut_cnv_cans)
mut_cnv_cans$fdr <- FDR_by_id_columns(p_vector = mut_cnv_cans$p,
                                      id_columns = c("SELF", "cancer", "variant_class", "affected_exp_type"),
                                      df = mut_cnv_cans)
mut_cnv_cans %>%
  filter(variant_class == "Missense_protein_level_high_vs_Truncation") %>%
  # filter(fdr < 0.1) %>%
  filter(SELF == "trans") %>%
  filter(affected_exp_type == "RNA") %>%
  filter(p < 0.05) %>%
  arrange(fdr)

write.table(x = mut_cnv_cans, file = paste0(makeOutDir(resultD = resultD), "TP53_missense_mut_impact_tab.txt"), sep = "\t", row.names = F, quote = F)

# anotate pair table to include the TFs of affected gene ------------------
pair_tab_previous <- mut_cnv_cans %>%
  select(GENE, SUB_GENE) %>%
  unique()
## add in the TFs of the affected gene
pair_tab_previous <- merge(pair_tab_previous, tf_tab2add %>%
                             select(target_genesymbol, source_genesymbol), by.x = c("SUB_GENE"), by.y = c("target_genesymbol"), all.x = T)
## filter for TFs that interact with TP53
pair_tab_previous <- pair_tab_previous %>%
  mutate(pair_pro_previous = paste0(GENE, ":", source_genesymbol)) %>%
  filter(source_genesymbol %in% TP53_pair_tab$SUB_GENE)

pair_tab_previous <- merge(pair_tab_previous, pair_tab_annotated %>%
                             filter(SUB_GENE.is_TF_downstream == T | SUB_GENE.is_complex_partner == T) %>%
                             select(pair_pro, SUB_GENE.is_TF_downstream, SUB_GENE.is_complex_partner),
                           by.x = c("pair_pro_previous"), by.y = c("pair_pro"), all.x = T)
pair_tab_previous$SUB_GENE.is_TF_downstream[is.na(pair_tab_previous$SUB_GENE.is_TF_downstream)] <- T
pair_tab_previous$SUB_GENE.is_complex_partner[is.na(pair_tab_previous$SUB_GENE.is_complex_partner)] <- F
pair_tab_previous <- pair_tab_previous %>%
  mutate(SUB_GENE_TF.is_TP53_downstream = SUB_GENE.is_TF_downstream) %>%
  mutate(SUB_GENE_TF.is_TP53_complex_partner = SUB_GENE.is_complex_partner) %>%
  mutate(pair_pro = paste0(GENE, ":", SUB_GENE)) %>%
  select(GENE, SUB_GENE, pair_pro, pair_pro_previous, SUB_GENE_TF.is_TP53_downstream, SUB_GENE_TF.is_TP53_complex_partner)
write.table(x = pair_tab_previous, file = paste0(makeOutDir(resultD = resultD), "pair_tab_w.TFs.txt"), sep = "\t", row.names = F, quote = F)

### check if SP1 and mutant TP53 interaction is within 

## check which cancer types have enough samples for the comparison