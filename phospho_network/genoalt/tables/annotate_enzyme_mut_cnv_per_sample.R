# Yige Wu @ WashU 2018 Aug
## for each kinase/phosphatase, annotate samples into kinase-mutated tumors, kinase-amplified tumors, kinase-deleted tumors against other tumors

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
amp_log2cn <- 0.3
del_log2cn <- -0.3
# inputs ------------------------------------------------------------------
## input enzyme-substrate table
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod/omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned_addDEPOD_extended.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source == "NetKIN" | (ptms_site_pairs_sup$Source != "NetKIN" & ptms_site_pairs_sup$networkin_score >= 5),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)
enzymes <- unique(ptms_site_pairs_sup$GENE)


cnv <- fread(input = "./cptac2p/copy_number/gatk4wxscnv/v1.3.swap_contamination_fixed/prospective_somatic/deliverables/ov/gene_level_CNV.ov.v1.3.2018-03-19.tsv", data.table = F)
int_enzyme_cnv <- intersect(unique(cnv$gene), enzymes)


# annotate samples --------------------------------------------
mut_cnv_num_cans <- NULL
mut_cnv_partIDs_cans <- list()
for (cancer in cancers_sort) {
  ## initiate columns
  num_mut <- vector(mode = "numeric", length = length(int_enzyme_cnv))
  num_amp <- vector(mode = "numeric", length = length(int_enzyme_cnv))
  num_del <- vector(mode = "numeric", length = length(int_enzyme_cnv))
  mut_cnv_partIDs_cans[[cancer]] <- list()
  
  ## input gene-level CNV
  if (cancer %in% c("OV", "CO")) {
    cnv <- fread(input = paste0("./cptac2p/copy_number/gatk4wxscnv/v1.3.swap_contamination_fixed/prospective_somatic/deliverables/", tolower(cancer),"/gene_level_CNV.", tolower(cancer),".v1.3.2018-03-19.tsv"), data.table = F)
  }
  if (cancer %in% c("BRCA")) {
    cnv <- fread(input = "./cptac2p/copy_number/gatk4wxscnv/v1.3.swap_contamination_fixed/prospective_somatic/deliverables/br/gene_level_CNV.br.v1.3.2018-03-19.tsv", data.table = F)
  }
  maf <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/", maf_files[grepl(x = maf_files, pattern = cancer)]), data.table = F)
  maf$partID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", n = 2)[,1]
  
  for (i in 1:length(int_enzyme_cnv)) {
    gene <- int_enzyme_cnv[i]
    mut_cnv_partIDs_cans[[cancer]][[gene]] <- list()
    cnv_gene <- cnv[cnv$gene == gene, -1]
    
    amp_vals <- cnv_gene[1, cnv_gene >= amp_log2cn]
    amp_partIDs <- names(amp_vals)
    mut_cnv_partIDs_cans[[cancer]][[gene]][["amp"]] <- amp_partIDs
    num_amp[i] <- length(amp_partIDs)
    
    del_vals <- cnv_gene[1, cnv_gene <= del_log2cn]
    del_partIDs <- names(del_vals)
    mut_cnv_partIDs_cans[[cancer]][[gene]][["del"]] <- del_partIDs
    num_del[i] <- length(del_partIDs)
    
    mut_gene <- maf[maf$Hugo_Symbol == gene & maf$Variant_Classification != "Silent",]
    mut_vals <- as.vector(mut_gene$Variant_Classification)
    mut_partIDs <- as.vector(mut_gene$partID)
    mut_cnv_partIDs_cans[[cancer]][[gene]][["mut"]] <- mut_partIDs
    mut_partIDs <- mut_partIDs[!(mut_partIDs %in% del_partIDs) | mut_partIDs %in% amp_partIDs]
    num_mut[i] <- length(mut_partIDs)
  }
  mut_cnv_num_can <- data.frame(gene = int_enzyme_cnv, num_mut = num_mut, num_amp = num_amp, num_del = num_del, cancer = cancer)
  mut_cnv_num_cans <- rbind(mut_cnv_num_cans, mut_cnv_num_can)
}
saveRDS(object = mut_cnv_partIDs_cans, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_partIDs_cans.RDS"))
write.table(x = mut_cnv_num_cans, file = paste0(makeOutDir(resultD = resultD),"mut_cnv_num_cans.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
