# Yige Wu @ WashU 2018 July
## visualize the mutation-impacted pairs, whether different kinase mutations are in the same samples, especially PRKDC, MAP2K4 and AKT1
## not the same samples

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
library(ggrepel)
library(dplyr)
library(readr)
library(gtable)

# set variables -----------------------------------------------------------
sig_thres <- 0.1
mut_enzymes <- c("AKT1", "PRKDC")
mut_enzymes <- unique(mutimpact_sig_driver_pho$Mutated_Gene)
mut_enzymes <- c("AKT1", mut_enzymes[mut_enzymes != "AKT1"])
substrates2p <- c("BAD", "GSK3B")
substrates2p <- c("BAD", "GSK3B", "GAB1")
genes2p <- c(mut_enzymes, substrates2p)
maf_files <- list.files(path = "./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/")

# inputs ------------------------------------------------------------------
mutimpact <- fread(input = paste0(ppnD, "genoalt/tables/merge_mutation_impact/mutation_impact.txt"), data.table = F)
mutimpact_sig <- mutimpact[mutimpact$p_value < sig_thres,]
mutimpact_sig_driver <- mutimpact_sig[(mutimpact_sig$enzyme_driver_type == "oncogene" & mutimpact_sig$Fold_Change > 0) | (mutimpact_sig$enzyme_driver_type == "tsg" & mutimpact_sig$Fold_Change < 0),]
mutimpact_sig_driver_pho <- mutimpact_sig_driver[mutimpact_sig_driver$substrate_type == "phosphosite",]
mutimpact_sig_driver_pho <- mutimpact_sig_driver_pho[order(mutimpact_sig_driver_pho$p_value, decreasing = F),]
clinical <- read.delim(paste0(cptac_sharedD, "Specimen_Data_20161005_Yige_20180307.txt"), sep = "\t")
clinical <- data.frame(clinical)
driver_genes <- read.delim(file = "./TCGA_data/reference_files/Consensus.Genelist.full.txt")

for (cancer in c("BRCA")) {
  maf <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/", maf_files[grepl(x = maf_files, pattern = cancer)]), data.table = F)
  maf_enzyme_sub <- maf[(maf$Hugo_Symbol %in% genes2p) & maf$Variant_Classification != "Silent",]
  
  ## order the samples by mutations
  sample_tab <- data.frame(table(maf_enzyme_sub[,c("Tumor_Sample_Barcode", "Hugo_Symbol")]))
  mut_score <- vector(mode = "numeric", length = nrow(sample_tab))
  for (i in 1:length(genes2p)) {
    mut_score[sample_tab$Hugo_Symbol == genes2p[i] & sample_tab$Freq > 0] <- length(genes2p) + 1 - i
  }
  sample_tab$mut_score <- mut_score
  sample_tab_tmp <- group_by(sample_tab, Tumor_Sample_Barcode)
  sample_scores <- summarise(sample_tab_tmp, sample_score = sum(mut_score))
  samples_ordered <- sample_scores$Tumor_Sample_Barcode[order(sample_scores$sample_score, decreasing = T)]
  maf_enzyme_sub$partID <- str_split_fixed(string = maf_enzyme_sub$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  
  mutimpact_sub <- mutimpact_sig_driver_pho[mutimpact_sig_driver_pho$Mutated_Gene %in% mut_enzymes & mutimpact_sig_driver_pho$Substrate_Gene %in% substrates2p,]
  mutimpact_sub$SUB_MOD_RSD <- toupper(mutimpact_sub$Phosphosite)
  pho_sub <- NULL
  ## get phosphorylation level of these phosphosites
  for (i in 1:nrow(mutimpact_sub)) {
    ## input phosphorylation data
    pho <- fread(paste0(cptac_sharedD, cancer, "/", prefix[cancer], "_", "PHO", "_formatted_normalized_noControl.txt"),
                 data.table = F)
    pho_head <- formatPhosphosite(phosphosite_vector = pho$Phosphosite, gene_vector = pho$Gene)
    
    sub <- mutimpact_sub[i, "Substrate_Gene"]
    rsd <- mutimpact_sub[i, "SUB_MOD_RSD"]
    
    pho_row <- pho[pho_head$SUBSTRATE == sub & pho_head$SUB_MOD_RSD == rsd,]
    pho_tmp <- melt(pho_row)
    pho_tmp$SUB_MOD_RSD <- rsd
    ## map sample IDs to patient IDs
    pho_tmp$partID <- sampID2partID(sampleID_vector = as.vector(pho_tmp$variable), sample_map = clinical)
    pho_sub <- rbind(pho_sub, pho_tmp)
  }
  pho_sub$site <- paste0(pho_sub$Gene, ":", pho_sub$SUB_MOD_RSD)
  lim = max(abs(max(pho_sub$value, na.rm = T)), abs(min(pho_sub$value, na.rm = T)))
  cap <- min(2, lim)
  pho_sub$pho_capped <- as.numeric(pho_sub$value)
  pho_sub$pho_capped[pho_sub$pho_capped > cap] <- cap
  pho_sub$pho_capped[pho_sub$pho_capped < (-cap)] <- (-cap)
  all_partIDs <- unique(pho_sub$partID)
  
  ## add patients without mutations
  maf_tmp <- maf_enzyme_sub[!duplicated(maf_enzyme_sub$Hugo_Symbol),]
  maf_tmp  <- maf_tmp[rep(1:length(unique(maf_enzyme_sub$Hugo_Symbol)), length(all_partIDs[!(all_partIDs %in% maf_enzyme_sub$partID)])),]
  maf_tmp$partID <- rep(all_partIDs[!(all_partIDs %in% maf_enzyme_sub$partID)], length(unique(maf_enzyme_sub$Hugo_Symbol)))
  maf_tmp$Variant_Classification <- NA
  maf_enzyme_sub$is_mut <- TRUE
  maf_tmp$is_mut <- FALSE
  maf2p <- rbind(maf_enzyme_sub, maf_tmp)
  partIDs_ordered <- str_split_fixed(string = samples_ordered, pattern = "_", 2)[,1]
  partIDs_ordered <- c(partIDs_ordered, all_partIDs[!(all_partIDs %in% partIDs_ordered)])
  
  tab2p <- maf2p
  tab2p$partID <- factor(tab2p$partID, levels = partIDs_ordered)
  tab2p$Hugo_Symbol <- factor(tab2p$Hugo_Symbol, levels = rev(genes2p))
  p <- ggplot()
  p <- p + geom_tile(data = tab2p, mapping = aes(x = partID, y = Hugo_Symbol, fill = Variant_Classification, color = is_mut))
  p <- p + scale_fill_manual(na.value = NA, values = set1)
  p <- p + scale_color_manual(na.value = NA, values = c("TRUE" = "black", "FALSE" = NA))
  p = p + theme(axis.text.x = element_text(angle = 90, size = 5))
  p
  p_mut <- p
  gp_mut <- ggplot_build(p_mut)
  gtabp_mut <- ggplot_gtable(gp_mut)
  np_mut <- length(gtabp_mut$layout$panel_ranges[[1]]$y.labels)
  
  tab2p <- pho_sub
  tab2p$partID <- factor(tab2p$partID, levels = partIDs_ordered)
  tab2p$Gene <- factor(tab2p$Gene, levels = rev(substrates2p))
  p <- ggplot()
  p <- p + geom_tile(data = tab2p, mapping = aes(x = partID, y = site, fill = pho_capped))
  p = p + scale_fill_gradientn(name= "phospho_level", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
  p = p + theme(axis.text.x = element_text(angle = 90, size = 5))
  p
  p_pho_sub <- p
  gp_pho_sub <- ggplot_build(p_pho_sub)
  gtabp_pho_sub <- ggplot_gtable(gp_pho_sub)
  np_pho_sub <- length(gtabp_pho_sub$layout$panel_ranges[[1]]$y.labels)
  
  g <- gtable:::rbind_gtable(gtabp_mut, gtabp_pho_sub, "last")

  panels <- g$layout$t[grep("panel", g$layout$name)]
  g$layout$name[grep("panel", g$layout$name)]
  
  ## adjust height of each panel
  g$heights[panels] <- unit(x = c(3, 3), units = "null")
  fn = paste0(makeOutDir(resultD = resultD), "test.pdf")
  grid.newpage()
  pdf(fn, height = 3, width = 10)
  grid.draw(g)
  dev.off()
}


for (cancer in c("BRCA")) {
  maf <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/", maf_files[grepl(x = maf_files, pattern = cancer)]), data.table = F)
  maf_enzyme <- maf[maf$Hugo_Symbol %in% mutimpact$Mutated_Gene & maf$Variant_Classification != "Silent",]

  tab2p <- maf_enzyme
  p <- ggplot()
  p <- p + geom_tile(data = tab2p, mapping = aes(x = Tumor_Sample_Barcode, y = Hugo_Symbol, fill = Variant_Classification))
  p 
  ggsave(filename = paste0(makeOutDir(resultD = resultD), cancer, "_drivers_somatic_mutations.pdf"), width = 10, height = 5)
  
  for (enzyme in unique(mutimpact_sig$Mutated_Gene)) {
    maf_sub <- maf[maf$Hugo_Symbol %in% mutimpact_sig$Substrate_Gene[mutimpact_sig$Mutated_Gene == enzyme] & maf$Variant_Classification != "Silent",]
    maf_enzyme_sub <- rbind(maf_enzyme[maf_enzyme$Hugo_Symbol == enzyme,], maf_sub)
    tab2p <- maf_enzyme_sub
    p <- ggplot()
    p <- p + geom_tile(data = tab2p, mapping = aes(x = Tumor_Sample_Barcode, y = Hugo_Symbol, fill = Variant_Classification))
    p 
    ggsave(filename = paste0(makeOutDir(resultD = resultD), cancer, "_", enzyme, "_and_substrates_somatic_mutations.pdf"), width = 10, height = 5)
  }
}

for (cancer in c("OV", "BRCA")) {
  maf <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/", maf_files[grepl(x = maf_files, pattern = cancer)]), data.table = F)
  maf_enzyme <- maf[maf$Hugo_Symbol %in% driver_genes$Gene[driver_genes$Cancer.type == cancer] & maf$Variant_Classification != "Silent",]
  
  tab2p <- maf_enzyme
  p <- ggplot()
  p <- p + geom_tile(data = tab2p, mapping = aes(x = Tumor_Sample_Barcode, y = Hugo_Symbol, fill = Variant_Classification))
  p 
  ggsave(filename = paste0(makeOutDir(resultD = resultD), cancer, "_drivers_somatic_mutations.pdf"), width = 10, height = 5)
}


for (cancer in c("CO")) {
  maf <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/", maf_files[grepl(x = maf_files, pattern = cancer)]), data.table = F)
  maf_enzyme <- maf[maf$Hugo_Symbol %in% driver_genes$Gene[driver_genes$Cancer.type == "COADREAD"] & maf$Variant_Classification != "Silent",]
  
  tab2p <- maf_enzyme
  p <- ggplot()
  p <- p + geom_tile(data = tab2p, mapping = aes(x = Tumor_Sample_Barcode, y = Hugo_Symbol, fill = Variant_Classification))
  p 
  ggsave(filename = paste0(makeOutDir(resultD = resultD), cancer, "_drivers_somatic_mutations.pdf"), width = 10, height = 5)
}