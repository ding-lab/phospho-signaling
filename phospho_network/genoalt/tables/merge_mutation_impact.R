# Yige Wu @ WashU 2018 July
## merging the mutation impact


# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# inputs ------------------------------------------------------------------
# driver_tab <- fread(input = "./PhosphoDrug/PhosphoDrug_shared_data/analysis_results/mutationimpact/07032018/Driver/correctedPvalueDriver", data.table = F)
# drug_tab <- fread(input = "./PhosphoDrug/PhosphoDrug_shared_data/analysis_results/mutationimpact/07032018/Druggable/correctedPvalueDruggable", data.table = F)

# driver_tab <- fread(input = "./PhosphoDrug/PhosphoDrug_shared_data/analysis_results/mutationimpact/07192018//Driver/correctedPvalueDriver", data.table = F)
drug_tab <- fread(input = "./PhosphoDrug/PhosphoDrug_shared_data/analysis_results/mutationimpact/07192018/Druggable/correctedPvalueDruggable", data.table = F)
# col_names <- colnames(driver_tab)[!(colnames(driver_tab) %in% c("adjusted_pvalue"))]
# mutimpact <-rbind(driver_tab[,col_names], drug_tab[,col_names])

# mutimpact <- fread(input = paste0("./PhosphoDrug/PhosphoDrug_shared_data/analysis_results/mutationimpact/07232018/correctedPvalueDriver"), data.table = F)
mutimpact <- fread(input = paste0("./PhosphoDrug/PhosphoDrug_shared_data/analysis_results/mutationimpact/07262018/correctedPvalueDriver"), data.table = F)
drivers <- read_delim("./TCGA_data/reference_files/Consensus.Genelist.full.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
oncogenes <- drivers$Gene[grepl(pattern = "oncogene", x = drivers$`Tumor suppressor or oncogene prediction (by 20/20+)`)]
tsgs <- drivers$Gene[grepl(pattern = "tsg", x = drivers$`Tumor suppressor or oncogene prediction (by 20/20+)`)]

## input enzyme-substrate table
ptms_sup_tab <- fread(input = "./cptac2p/analysis_results/phospho_network/compile_enzyme_substrate/tables/compile_omnipath_networkin_depod/omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned_addDEPOD.csv", data.table = F)
# set variables -----------------------------------------------------------
sig_thres <- 0.05

# annotate -----------------------------------------------------------------
mutimpact$id <- paste0(mutimpact$Cancer, mutimpact$Mutated_Gene, mutimpact$Substrate_Gene, mutimpact$Phosphosite)
mutimpact <- mutimpact[order(mutimpact$p_value),]
mutimpact <- mutimpact[!duplicated(mutimpact$id),]
mutimpact$id <- NULL
mutimpact$SELF <- ifelse(as.vector(mutimpact$Mutated_Gene) == as.vector(mutimpact$Substrate_Gene), "cis", "trans")
mutimpact$substrate_type <- ifelse(mutimpact$Phosphosite == "Protein", "protein", "phosphosite")
enzyme_driver_type <- vector(mode = "character", length = nrow(mutimpact))
enzyme_driver_type <- ""
enzyme_driver_type[mutimpact$Mutated_Gene %in% oncogenes]  <- "oncogene"
enzyme_driver_type[mutimpact$Mutated_Gene %in% tsgs]  <- "tsg"
mutimpact$enzyme_driver_type <- enzyme_driver_type
mutimpact$enzyme_is_druggable <- (mutimpact$Mutated_Gene %in% drug_tab$Mutated_Gene)
mutimpact$substrate_direction <- ifelse(mutimpact$Fold_Change > 0, "up", "down")
mutimpact$SUB_MOD_RSD <- toupper(mutimpact$Phosphosite)
mutimpact$pair <- paste0(mutimpact$Mutated_Gene, ":", mutimpact$Substrate_Gene, ":", mutimpact$SUB_MOR_RSD)
mutimpact$known_es_pair <- (mutimpact$pair %in% ptms_sup_tab$pair)
mutimpact_known <- mutimpact[mutimpact$known_es_pair,]
mutimpact_driver <- mutimpact[!is.na(mutimpact$enzyme_driver_type),]
mutimpact_driver_known <- mutimpact_known[!is.na(mutimpact_known$enzyme_driver_type),]
# adjust p-value ----------------------------------------------------------
## I think we should divide mutations affecting protein and phosphosites saperately and cis and trans separately, in concord with our adjusting method as regression
mutimpact_adjusted <- NULL
mutimpact_known_adjusted <- NULL
mutimpact_driver_adjusted <- NULL
mutimpact_driver_known_adjusted <- NULL

for (self in unique(mutimpact$SELF)) {
  for (substrate_type in unique(mutimpact$substrate_type)) {
    for (cancer in unique(mutimpact$Cancer)) {
      mutimpact2adjust <- mutimpact[mutimpact$SELF == self & mutimpact$substrate_type == substrate_type & mutimpact$Cancer == cancer,]
      mutimpact2adjust$adjusted_pvalue <- p.adjust(as.vector(mutimpact2adjust$p_value), method = "fdr")
      mutimpact_adjusted <- rbind(mutimpact_adjusted, mutimpact2adjust)
      
      mutimpact_known2adjust <- mutimpact_known[mutimpact_known$SELF == self & mutimpact_known$substrate_type == substrate_type & mutimpact_known$Cancer == cancer,]
      mutimpact_known2adjust$adjusted_pvalue <- p.adjust(as.vector(mutimpact_known2adjust$p_value), method = "fdr")
      mutimpact_known_adjusted <- rbind(mutimpact_known_adjusted, mutimpact_known2adjust)
      
      mutimpact_driver2adjust <- mutimpact_driver[mutimpact_driver$SELF == self & mutimpact_driver$substrate_type == substrate_type & mutimpact_driver$Cancer == cancer,]
      mutimpact_driver2adjust$adjusted_pvalue <- p.adjust(as.vector(mutimpact_driver2adjust$p_value), method = "fdr")
      mutimpact_driver_adjusted <- rbind(mutimpact_driver_adjusted, mutimpact_driver2adjust)
      
      mutimpact_driver_known2adjust <- mutimpact_driver_known[mutimpact_driver_known$SELF == self & mutimpact_driver_known$substrate_type == substrate_type & mutimpact_driver_known$Cancer == cancer,]
      mutimpact_driver_known2adjust$adjusted_pvalue <- p.adjust(as.vector(mutimpact_driver_known2adjust$p_value), method = "fdr")
      mutimpact_driver_known_adjusted <- rbind(mutimpact_driver_known_adjusted, mutimpact_driver_known2adjust)
    }
  }
}
mutimpact_adjusted <- mutimpact_adjusted[order(mutimpact_adjusted$adjusted_pvalue),]
write.table(x = mutimpact, file = paste0(makeOutDir(resultD = resultD), "mutation_impact.txt"), quote = F, row.names = F)
mutimpact_known_adjusted <- mutimpact_known_adjusted[order(mutimpact_known_adjusted$adjusted_pvalue),]
write.table(x = mutimpact_known, file = paste0(makeOutDir(resultD = resultD), "mutation_impact_known.txt"), quote = F, row.names = F)
mutimpact_driver_adjusted <- mutimpact_driver_adjusted[order(mutimpact_driver_adjusted$adjusted_pvalue),]
mutimpact_driver_known_adjusted <- mutimpact_driver_known_adjusted[order(mutimpact_driver_known_adjusted$adjusted_pvalue),]
