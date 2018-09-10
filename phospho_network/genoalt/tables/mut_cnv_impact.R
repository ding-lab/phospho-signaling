# Yige Wu @ WashU 2018 Aug
## for each kinase-substrate pairs, compare the substrate phosphorylation level of kinase-mutated tumors, kinase-amplified tumors, kinase-deleted tumors against other tumors
## using wilcox test
## reference: https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
amp_log2cn <- 0.3
del_log2cn <- -0.3
num_control_thres <- 4
num_genoalt_thres <- 4

# inputs ------------------------------------------------------------------
## input enzyme-substrate table
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod/omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned_addDEPOD_extended.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source == "NetKIN" | (ptms_site_pairs_sup$Source != "NetKIN" & ptms_site_pairs_sup$networkin_score >= 5),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)

mut_cnv_partIDs_cans <- readRDS(file = "./cptac2p/analysis_results/phospho_network/genoalt/tables/annotate_enzyme_mut_cnv_per_sample/mut_cnv_partIDs_cans.RDS")
enzyme_sub <- unique(ptms_site_pairs_sup[ptms_site_pairs_sup$GENE %in% names(mut_cnv_partIDs_cans[["BRCA"]]), c("GENE", "SUB_GENE")])

for (cancer in c("BRCA")) {
  pho_data <- fread(input = paste0(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt"), data.table = F)
  sampIDs <- colnames(pho_data)[!(colnames(pho_data) %in% c("Gene", "Phosphosite"))]
  pho_data <- pho_data[rowSums(!is.na(pho_data[, sampIDs])) >= (num_control_thres + num_genoalt_thres),]
  pho_head <- formatPhosphosite(pho_data$Phosphosite, pho_data$Gene)
  partIDs <- sampID2partID(sampleID_vector = sampIDs, sample_map = loadSampMap())
  names(pho_data) <- c("Gene", "Phosphosite", partIDs)
  
  ## initiate columns
  enzyme_sub_sites <- merge(enzyme_sub, pho_head[, c("SUBSTRATE", "SUB_MOD_RSD")], by.x = c("SUB_GENE"), by.y = c("SUBSTRATE"), all.y = T)
  enzyme_sub_sites <- enzyme_sub_sites[!is.na(enzyme_sub_sites$GENE) & !is.na(enzyme_sub_sites$SUB_MOD_RSD),]
  enzyme_sub_sites <- enzyme_sub_sites[order(enzyme_sub_sites$GENE, enzyme_sub_sites$SUB_GENE),]
  num_control <- vector(mode = "numeric", length = nrow(enzyme_sub_sites))
  
  num_mut <- num_control; meddiff_mut <- num_control; p_mut <- num_control; meddiff_bottom_mut <- num_control; meddiff_upper_mut <- num_control
  num_amp <- num_control; meddiff_amp <- num_control; p_amp <- num_control; meddiff_bottom_amp <- num_control; meddiff_upper_amp <- num_control
  num_del <- num_control; meddiff_del <- num_control; p_del <- num_control; meddiff_bottom_del <- num_control; meddiff_upper_del <- num_control
  
  for (enzyme in unique(enzyme_sub_sites$GENE)) {
    mut_cnv_partIDs <- c(mut_cnv_partIDs_cans[[cancer]][[enzyme]][["mut"]], mut_cnv_partIDs_cans[[cancer]][[enzyme]][["amp"]], mut_cnv_partIDs_cans[[cancer]][[enzyme]][["del"]])
    control_partIDs <- partIDs[!(partIDs %in% mut_cnv_partIDs)]
    
    mut_partIDs1 <- intersect(partIDs, mut_cnv_partIDs_cans[[cancer]][[enzyme]][["mut"]])
    amp_partIDs1 <- intersect(partIDs, mut_cnv_partIDs_cans[[cancer]][[enzyme]][["amp"]])
    del_partIDs1 <- intersect(partIDs, mut_cnv_partIDs_cans[[cancer]][[enzyme]][["del"]])
    
    print(paste0("enzyme:", enzyme))
    for (substrate in unique(enzyme_sub_sites$SUB_GENE[enzyme_sub_sites$GENE == enzyme])) {
      # print(paste0("substrate:", substrate))
      
      for (site in unique(enzyme_sub_sites$SUB_MOD_RSD[enzyme_sub_sites$GENE == enzyme & enzyme_sub_sites$SUB_GENE == substrate])) {
        # print(paste0("site:", site))
        
        i <- (enzyme_sub_sites$GENE == enzyme & enzyme_sub_sites$SUB_GENE == substrate & enzyme_sub_sites$SUB_MOD_RSD == site)
        
        pho_sub <- pho_data[pho_head$SUBSTRATE == substrate & pho_head$SUB_MOD_RSD == site, partIDs]
        pho_sub <- pho_sub[1,]
        control_pho <- pho_sub[, control_partIDs]
        control_pho <- control_pho[!is.na(control_pho)]
        num_control[i] <- length(control_pho)
        
        if (num_control[i] >= num_control_thres) {
          mut_pho <- pho_sub[, mut_partIDs1]
          mut_pho <- mut_pho[!is.na(mut_pho)]
          # col_tmp <- which(!is.na(mut_pho)); 
          num_mut[i] <- length(mut_pho)
          if (num_mut[i] > 0) {
            meddiff_mut[i] <- median(mut_pho) - median(control_pho)
            if (num_mut[i] >= num_genoalt_thres){
              stat <- wilcox.test(x = mut_pho, y = control_pho, conf.int = T)
              p_mut[i] <- stat$p.value
              meddiff_bottom_mut[i] <- stat$conf.int[1]
              meddiff_upper_mut[i] <- stat$conf.int[2]
            }
          }
          amp_pho <- pho_sub[, amp_partIDs1]
          # col_tmp <- which(!is.na(amp_pho))
          amp_pho <- amp_pho[!is.na(amp_pho)]
          num_amp[i] <- length(amp_pho)
          if (num_amp[i] > 0) {
            meddiff_amp[i] <- median(amp_pho) - median(control_pho)
            if (num_amp[i] >= num_genoalt_thres){
              stat <- wilcox.test(x = amp_pho, y = control_pho, conf.int = T)
              p_amp[i] <- stat$p.value
              meddiff_bottom_amp[i] <- stat$conf.int[1]
              meddiff_upper_amp[i] <- stat$conf.int[2]
            }
          }
          del_pho <- pho_sub[, del_partIDs1]
          # col_tmp <- which(!is.na(del_pho))
          del_pho <- del_pho[!is.na(del_pho)]
          num_del[i] <- length(del_pho)
          if (num_del[i] > 0) {
            meddiff_del[i] <- median(del_pho) - median(control_pho)
            if (num_del[i] >= num_genoalt_thres){
              stat <- wilcox.test(x = del_pho, y = control_pho, conf.int = T)
              p_del[i] <- stat$p.value
              meddiff_bottom_del[i] <- stat$conf.int[1]
              meddiff_upper_del[i] <- stat$conf.int[2]
            }
          }
        }
        
      }
    }
    print(which(i))
  }
  mut_cnv_tab <- enzyme_sub_sites
  mut_cnv_tab <- cbind(mut_cnv_tab, data.frame(num_mut, num_amp, num_del, num_control,
                                               p_mut, meddiff_bottom_mut, meddiff_upper_mut, meddiff_mut,
                                               p_amp, meddiff_bottom_amp, meddiff_upper_amp, meddiff_amp,
                                               p_del, meddiff_bottom_del, meddiff_upper_del, meddiff_del))
  mut_cnv_tab$cancer <- cancer
  write.table(x = mut_cnv_tab, file = paste0(makeOutDir(resultD = resultD), cancer, "_mut_cnv_tab.txt"), col.names = T, row.names = F, quote = F)
}

for (cancer in c("OV")) {
  pho_data <- fread(input = paste0(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt"), data.table = F)
  sampIDs <- colnames(pho_data)[!(colnames(pho_data) %in% c("Gene", "Phosphosite"))]
  pho_data <- pho_data[rowSums(!is.na(pho_data[, sampIDs])) >= (num_control_thres + num_genoalt_thres),]
  pho_head <- formatPhosphosite(pho_data$Phosphosite, pho_data$Gene)
  partIDs <- sampID2partID(sampleID_vector = sampIDs, sample_map = loadSampMap())
  names(pho_data) <- c("Gene", "Phosphosite", partIDs)
  
  ## initiate columns
  enzyme_sub_sites <- merge(enzyme_sub, pho_head[, c("SUBSTRATE", "SUB_MOD_RSD")], by.x = c("SUB_GENE"), by.y = c("SUBSTRATE"), all.y = T)
  enzyme_sub_sites <- enzyme_sub_sites[!is.na(enzyme_sub_sites$GENE) & !is.na(enzyme_sub_sites$SUB_MOD_RSD),]
  enzyme_sub_sites <- enzyme_sub_sites[order(enzyme_sub_sites$GENE, enzyme_sub_sites$SUB_GENE),]
  num_control <- vector(mode = "numeric", length = nrow(enzyme_sub_sites))
  
  num_mut <- num_control; meddiff_mut <- num_control; p_mut <- num_control; meddiff_bottom_mut <- num_control; meddiff_upper_mut <- num_control
  num_amp <- num_control; meddiff_amp <- num_control; p_amp <- num_control; meddiff_bottom_amp <- num_control; meddiff_upper_amp <- num_control
  num_del <- num_control; meddiff_del <- num_control; p_del <- num_control; meddiff_bottom_del <- num_control; meddiff_upper_del <- num_control
  
  for (enzyme in unique(enzyme_sub_sites$GENE)) {
    mut_cnv_partIDs <- c(mut_cnv_partIDs_cans[[cancer]][[enzyme]][["mut"]], mut_cnv_partIDs_cans[[cancer]][[enzyme]][["amp"]], mut_cnv_partIDs_cans[[cancer]][[enzyme]][["del"]])
    control_partIDs <- partIDs[!(partIDs %in% mut_cnv_partIDs)]
    
    mut_partIDs1 <- intersect(partIDs, mut_cnv_partIDs_cans[[cancer]][[enzyme]][["mut"]])
    amp_partIDs1 <- intersect(partIDs, mut_cnv_partIDs_cans[[cancer]][[enzyme]][["amp"]])
    del_partIDs1 <- intersect(partIDs, mut_cnv_partIDs_cans[[cancer]][[enzyme]][["del"]])
    
    print(paste0("enzyme:", enzyme))
    for (substrate in unique(enzyme_sub_sites$SUB_GENE[enzyme_sub_sites$GENE == enzyme])) {
      # print(paste0("substrate:", substrate))
      
      for (site in unique(enzyme_sub_sites$SUB_MOD_RSD[enzyme_sub_sites$GENE == enzyme & enzyme_sub_sites$SUB_GENE == substrate])) {
        # print(paste0("site:", site))
        
        i <- (enzyme_sub_sites$GENE == enzyme & enzyme_sub_sites$SUB_GENE == substrate & enzyme_sub_sites$SUB_MOD_RSD == site)
        
        pho_sub <- pho_data[pho_head$SUBSTRATE == substrate & pho_head$SUB_MOD_RSD == site, partIDs]
        pho_sub <- pho_sub[1,]
        control_pho <- pho_sub[, control_partIDs]
        control_pho <- control_pho[!is.na(control_pho)]
        num_control[i] <- length(control_pho)
        
        if (num_control[i] >= num_control_thres) {
          mut_pho <- pho_sub[, mut_partIDs1]
          mut_pho <- mut_pho[!is.na(mut_pho)]
          # col_tmp <- which(!is.na(mut_pho)); 
          num_mut[i] <- length(mut_pho)
          if (num_mut[i] > 0) {
            meddiff_mut[i] <- median(mut_pho) - median(control_pho)
            if (num_mut[i] >= num_genoalt_thres){
              stat <- wilcox.test(x = mut_pho, y = control_pho, conf.int = T)
              p_mut[i] <- stat$p.value
              meddiff_bottom_mut[i] <- stat$conf.int[1]
              meddiff_upper_mut[i] <- stat$conf.int[2]
            }
          }
          amp_pho <- pho_sub[, amp_partIDs1]
          # col_tmp <- which(!is.na(amp_pho))
          amp_pho <- amp_pho[!is.na(amp_pho)]
          num_amp[i] <- length(amp_pho)
          if (num_amp[i] > 0) {
            meddiff_amp[i] <- median(amp_pho) - median(control_pho)
            if (num_amp[i] >= num_genoalt_thres){
              stat <- wilcox.test(x = amp_pho, y = control_pho, conf.int = T)
              p_amp[i] <- stat$p.value
              meddiff_bottom_amp[i] <- stat$conf.int[1]
              meddiff_upper_amp[i] <- stat$conf.int[2]
            }
          }
          del_pho <- pho_sub[, del_partIDs1]
          # col_tmp <- which(!is.na(del_pho))
          del_pho <- del_pho[!is.na(del_pho)]
          num_del[i] <- length(del_pho)
          if (num_del[i] > 0) {
            meddiff_del[i] <- median(del_pho) - median(control_pho)
            if (num_del[i] >= num_genoalt_thres){
              stat <- wilcox.test(x = del_pho, y = control_pho, conf.int = T)
              p_del[i] <- stat$p.value
              meddiff_bottom_del[i] <- stat$conf.int[1]
              meddiff_upper_del[i] <- stat$conf.int[2]
            }
          }
        }
        
      }
    }
    print(which(i))
  }
  mut_cnv_tab <- enzyme_sub_sites
  mut_cnv_tab <- cbind(mut_cnv_tab, data.frame(num_mut, num_amp, num_del, num_control, meddiff_mut, p_mut, meddiff_amp, p_amp, meddiff_del, p_del))
  mut_cnv_tab$cancer <- cancer
  write.table(x = mut_cnv_tab, file = paste0(makeOutDir(resultD = resultD), cancer, "_mut_cnv_tab.txt"), col.names = T, row.names = F, quote = F)
}

for (cancer in c("CO")) {
  pho_data <- fread(input = paste0(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt"), data.table = F)
  sampIDs <- colnames(pho_data)[!(colnames(pho_data) %in% c("Gene", "Phosphosite"))]
  pho_data <- pho_data[rowSums(!is.na(pho_data[, sampIDs])) >= (num_control_thres + num_genoalt_thres),]
  pho_head <- formatPhosphosite(pho_data$Phosphosite, pho_data$Gene)
  partIDs <- sampID2partID(sampleID_vector = sampIDs, sample_map = loadSampMap())
  names(pho_data) <- c("Gene", "Phosphosite", partIDs)
  
  ## initiate columns
  enzyme_sub_sites <- merge(enzyme_sub, pho_head[, c("SUBSTRATE", "SUB_MOD_RSD")], by.x = c("SUB_GENE"), by.y = c("SUBSTRATE"), all.y = T)
  enzyme_sub_sites <- enzyme_sub_sites[!is.na(enzyme_sub_sites$GENE) & !is.na(enzyme_sub_sites$SUB_MOD_RSD),]
  enzyme_sub_sites <- enzyme_sub_sites[order(enzyme_sub_sites$GENE, enzyme_sub_sites$SUB_GENE),]
  num_control <- vector(mode = "numeric", length = nrow(enzyme_sub_sites))
  
  num_mut <- num_control; meddiff_mut <- num_control; p_mut <- num_control; meddiff_bottom_mut <- num_control; meddiff_upper_mut <- num_control
  num_amp <- num_control; meddiff_amp <- num_control; p_amp <- num_control; meddiff_bottom_amp <- num_control; meddiff_upper_amp <- num_control
  num_del <- num_control; meddiff_del <- num_control; p_del <- num_control; meddiff_bottom_del <- num_control; meddiff_upper_del <- num_control
  
  for (enzyme in unique(enzyme_sub_sites$GENE)) {
    mut_cnv_partIDs <- c(mut_cnv_partIDs_cans[[cancer]][[enzyme]][["mut"]], mut_cnv_partIDs_cans[[cancer]][[enzyme]][["amp"]], mut_cnv_partIDs_cans[[cancer]][[enzyme]][["del"]])
    control_partIDs <- partIDs[!(partIDs %in% mut_cnv_partIDs)]
    
    mut_partIDs1 <- intersect(partIDs, mut_cnv_partIDs_cans[[cancer]][[enzyme]][["mut"]])
    amp_partIDs1 <- intersect(partIDs, mut_cnv_partIDs_cans[[cancer]][[enzyme]][["amp"]])
    del_partIDs1 <- intersect(partIDs, mut_cnv_partIDs_cans[[cancer]][[enzyme]][["del"]])
    
    print(paste0("enzyme:", enzyme))
    for (substrate in unique(enzyme_sub_sites$SUB_GENE[enzyme_sub_sites$GENE == enzyme])) {
      # print(paste0("substrate:", substrate))
      
      for (site in unique(enzyme_sub_sites$SUB_MOD_RSD[enzyme_sub_sites$GENE == enzyme & enzyme_sub_sites$SUB_GENE == substrate])) {
        # print(paste0("site:", site))
        
        i <- (enzyme_sub_sites$GENE == enzyme & enzyme_sub_sites$SUB_GENE == substrate & enzyme_sub_sites$SUB_MOD_RSD == site)
        
        pho_sub <- pho_data[pho_head$SUBSTRATE == substrate & pho_head$SUB_MOD_RSD == site, partIDs]
        pho_sub <- pho_sub[1,]
        control_pho <- pho_sub[, control_partIDs]
        control_pho <- control_pho[!is.na(control_pho)]
        num_control[i] <- length(control_pho)
        
        if (num_control[i] >= num_control_thres) {
          mut_pho <- pho_sub[, mut_partIDs1]
          mut_pho <- mut_pho[!is.na(mut_pho)]
          # col_tmp <- which(!is.na(mut_pho)); 
          num_mut[i] <- length(mut_pho)
          if (num_mut[i] > 0) {
            meddiff_mut[i] <- median(mut_pho) - median(control_pho)
            if (num_mut[i] >= num_genoalt_thres){
              stat <- wilcox.test(x = mut_pho, y = control_pho, conf.int = T)
              p_mut[i] <- stat$p.value
              meddiff_bottom_mut[i] <- stat$conf.int[1]
              meddiff_upper_mut[i] <- stat$conf.int[2]
            }
          }
          amp_pho <- pho_sub[, amp_partIDs1]
          # col_tmp <- which(!is.na(amp_pho))
          amp_pho <- amp_pho[!is.na(amp_pho)]
          num_amp[i] <- length(amp_pho)
          if (num_amp[i] > 0) {
            meddiff_amp[i] <- median(amp_pho) - median(control_pho)
            if (num_amp[i] >= num_genoalt_thres){
              stat <- wilcox.test(x = amp_pho, y = control_pho, conf.int = T)
              p_amp[i] <- stat$p.value
              meddiff_bottom_amp[i] <- stat$conf.int[1]
              meddiff_upper_amp[i] <- stat$conf.int[2]
            }
          }
          del_pho <- pho_sub[, del_partIDs1]
          # col_tmp <- which(!is.na(del_pho))
          del_pho <- del_pho[!is.na(del_pho)]
          num_del[i] <- length(del_pho)
          if (num_del[i] > 0) {
            meddiff_del[i] <- median(del_pho) - median(control_pho)
            if (num_del[i] >= num_genoalt_thres){
              stat <- wilcox.test(x = del_pho, y = control_pho, conf.int = T)
              p_del[i] <- stat$p.value
              meddiff_bottom_del[i] <- stat$conf.int[1]
              meddiff_upper_del[i] <- stat$conf.int[2]
            }
          }
        }
        
      }
    }
    print(which(i))
  }
  mut_cnv_tab <- enzyme_sub_sites
  mut_cnv_tab <- cbind(mut_cnv_tab, data.frame(num_mut, num_amp, num_del, num_control,
                                               p_mut, meddiff_bottom_mut, meddiff_upper_mut, meddiff_mut,
                                               p_amp, meddiff_bottom_amp, meddiff_upper_amp, meddiff_amp,
                                               p_del, meddiff_bottom_del, meddiff_upper_del, meddiff_del))
  mut_cnv_tab$cancer <- cancer
  write.table(x = mut_cnv_tab, file = paste0(makeOutDir(resultD = resultD), cancer, "_mut_cnv_tab.txt"), 
              col.names = T, row.names = F, quote = F, sep = "\t")
}


