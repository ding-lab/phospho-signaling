# Yige Wu @ WashU 2018 Jul
## check the phosphosite abundance of samples with/without hotpho mutations in the same hybrid cluster

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# remove.packages(c("ggplot2", "ggrepel"))

library(ggrepel)
library(ggplot2)
library(readxl)
# inputs ------------------------------------------------------------------
color_direction <- c(set1[1], set1[2], "black")
names(color_direction) <- c("up", "down", "NA")
clinical <- fread(input = paste0(cptac_sharedD, "Specimen_Data_20161005_Yige_20180307.txt"), data.table = F)
maf_files <- list.files(path = "./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/")
subtype_map <- read_excel(paste0(cptac_sharedD, "5_CPTAC2_Breast_Prospective_Collection_BI/proteome-data-v1.01011/data/20171106_CPTAC2_ProspectiveBreastCancer_Sample_Annotations_Table_v79.xlsx"))


# input hotpho hybrid cluster annotation file -----------------------------
hotpho <- read_delim(file = paste0(ppnD, 'genoalt/tables/check_hotpho_cptac2p/', 'Data_201807_cc.p0.05.cluster_transcriptSynced_CPTAC2p_MAF_availability_annotated.tsv'),
                     "\t", escape_double = FALSE, col_types = cols(Cluster = col_character()), trim_ws = TRUE)
cptac_samples <- readRDS(file = paste0(ppnD, 'genoalt/tables/check_hotpho_cptac2p/', 'hotpho_hybrid_clusters_in_cptac_samples.RDS'))

# input phosphosite genomic coordinates (hg38) -----------------------------------
pho_head_geno_pos <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid.txt"), data.table = F)
# pho_head_geno_pos <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid_hg38_hg19.txt"), data.table = F)
pho_head_geno_pos$id <- paste0(pho_head_geno_pos$SUBSTRATE, ":", pho_head_geno_pos$SUB_MOD_RSD)
pho_head_geno_pos <- pho_head_geno_pos[!duplicated(pho_head_geno_pos$id),]
rownames(pho_head_geno_pos) <- pho_head_geno_pos$id
pho_head_geno_pos$geno_pos <- paste0(pho_head_geno_pos$seq_region_name, ":", pho_head_geno_pos$start, "-", pho_head_geno_pos$end)

# trans-regulated pairs overlap with driver mutations --------------
for (cancer in c("BRCA", "CO")) {
  hotpho_mut_ptm <- hotpho[hotpho[, paste0("mutation_and_ptm_in_", cancer)] & !is.na(hotpho[, paste0("mutation_and_ptm_in_", cancer)]),]
  clusters_mut_ptm <- unique(hotpho_mut_ptm$Cluster)
  
  ## input phospho level
  pho_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
  pho_head <- formatPhosphosite(phosphosite_vector = as.vector(pho_data$Phosphosite), gene_vector = as.vector(pho_data$Gene))
  pho_head$id <- paste0(pho_head$SUBSTRATE, ":", pho_head$SUB_MOD_RSD)
  pho_head$geno_pos <- pho_head_geno_pos[as.vector(pho_head$id), "geno_pos"]
  
  ## input mutation data
  maf <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/", maf_files[grepl(x = maf_files, pattern = cancer)]), data.table = F)
  maf$geno_pos <- paste0(str_split_fixed(string = maf$Chromosome, pattern = "chr", 2)[,2], ":", maf$Start_Position, "-", maf$End_Position)

  for (cluster in clusters_mut_ptm) {
    hotpho_cluster <- hotpho_mut_ptm[hotpho_mut_ptm$Cluster == cluster,]
    hotpho_cluster_ptm <- hotpho_cluster[hotpho_cluster$Alternate == "ptm",]
    hotpho_cluster_mut <- hotpho_cluster[hotpho_cluster$Alternate != "ptm",]
    for (i in 1:nrow(hotpho_cluster_ptm)) {
      gene_ptm <- as.character(hotpho_cluster_ptm[i, "Gene_Drug"])
      rsd_ptm <- as.character(hotpho_cluster_ptm[i, "SUB_MOD_RSD"])
      geno_ptm <- as.character(hotpho_cluster_ptm[i, "geno_pos"])
      for (j in 1:nrow(hotpho_cluster_mut)) {
        gene_mut <- as.character(hotpho_cluster_mut[j, "Gene_Drug"])
        aa_mut <- as.character(hotpho_cluster_mut[j, "Mutation_Gene"])
        geno_mut <- as.character(hotpho_cluster_mut[j, "geno_pos"])
        
        pho_ptm <- pho_data[pho_head$SUBSTRATE == gene_ptm & !is.na(pho_head$geno_pos) & pho_head$geno_pos == geno_ptm,]
        pho_ptm.m <- melt(pho_ptm)
        pho_ptm.m$partID <- sampID2partID(sampleID_vector = as.vector(pho_ptm.m$variable), sample_map =  clinical)
        
        maf_mut <- maf[maf$Hugo_Symbol == gene_mut & maf$geno_pos == geno_mut,]
        maf_mut$partID <- str_split_fixed(maf_mut$Tumor_Sample_Barcode, "_", 2)[,1]
        pho_ptm.m$mutation <- ifelse(test = (pho_ptm.m$partID %in% maf_mut$partID), yes = aa_mut, no = "")
        pho_ptm.m$cancer <- cancer
        
        tab2p <- pho_ptm.m
        p = ggplot(tab2p, aes(x=cancer, y=value, color = mutation, label= as.character(mutation)))
        p = p + geom_point(aes(shape = ifelse(mutation != "", "b", "a")), 
                           position = position_jitterdodge(), stroke = 0, alpha = 0.8)
        p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
        p = p + geom_text_repel(aes(segment.color = mutation),
                                force = 1,
                                segment.size = 0.5, segment.alpha = 1,
                                size=2,alpha=0.8)
        p = p + labs(y=paste0(gene_ptm, " ", rsd_ptm, " phosphorylation abundance(log2 ratio"))
        p = p + theme_nogrid()
        p = p + theme(axis.title = element_text(size=6), legend.position = 'none',
                      axis.text.x = element_text(colour="black", size=8, vjust=0.5),
                      axis.text.y = element_text(colour="black", size=8))
        p = p + theme(title = element_text(size = 8))
        p = p + theme(axis.title.x = element_blank())
        p = p + theme(axis.title.y = element_text(size = 10))
        fn = paste0(makeOutDir(resultD = resultD), "cluster_", cluster, "_", gene_mut, "_", aa_mut, "-", gene_ptm, "_", rsd_ptm, ".pdf")
        pdf(file = fn, height=4, width=4)
        print(p)
        dev.off()
      }
    }
  }
}



# plot downstream of AKT1 -------------------------------------------------
for (cancer in c("BRCA")) {
  hotpho_mut_ptm <- hotpho[hotpho[, paste0("mutation_and_ptm_in_", cancer)] & !is.na(hotpho[, paste0("mutation_and_ptm_in_", cancer)]),]
  clusters_mut_ptm <- unique(hotpho_mut_ptm$Cluster)
  
  ## input phospho level
  pho_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
  pho_head <- formatPhosphosite(phosphosite_vector = as.vector(pho_data$Phosphosite), gene_vector = as.vector(pho_data$Gene))
  pho_head$id <- paste0(pho_head$SUBSTRATE, ":", pho_head$SUB_MOD_RSD)
  pho_head$geno_pos <- pho_head_geno_pos[as.vector(pho_head$id), "geno_pos"]
  
  ## input mutation data
  maf <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/", maf_files[grepl(x = maf_files, pattern = cancer)]), data.table = F)
  maf$geno_pos <- paste0(str_split_fixed(string = maf$Chromosome, pattern = "chr", 2)[,2], ":", maf$Start_Position, "-", maf$End_Position)
  
  for (cluster in c("756.0")) {
    hotpho_cluster <- hotpho_mut_ptm[hotpho_mut_ptm$Cluster == cluster,]
    hotpho_cluster_ptm <- hotpho_cluster[hotpho_cluster$Alternate == "ptm",]
    hotpho_cluster_mut <- hotpho_cluster[hotpho_cluster$Alternate != "ptm",]
    for (i in 1:nrow(hotpho_cluster_ptm)) {
      gene_ptm <- as.character(hotpho_cluster_ptm[i, "Gene_Drug"])
      rsd_ptm <- as.character(hotpho_cluster_ptm[i, "SUB_MOD_RSD"])
      geno_ptm <- as.character(hotpho_cluster_ptm[i, "geno_pos"])
      for (j in 1:nrow(hotpho_cluster_mut)) {
        gene_mut <- as.character(hotpho_cluster_mut[j, "Gene_Drug"])
        aa_mut <- as.character(hotpho_cluster_mut[j, "Mutation_Gene"])
        geno_mut <- as.character(hotpho_cluster_mut[j, "geno_pos"])
        
        pho_ptm <- pho_data[pho_head$SUBSTRATE == gene_ptm & !is.na(pho_head$geno_pos) & pho_head$geno_pos == geno_ptm,]
        pho_ptm.m <- melt(pho_ptm)
        
        gene_sub <- "GSK3B"
        rsd_sub <- "S9"
        pho_sub <- pho_data[pho_head$SUBSTRATE == gene_sub & pho_head$SUB_MOD_RSD == rsd_sub,]
        pho_sub.m <- melt(pho_sub)
        
        pho_ptm.m <- merge(pho_ptm.m, pho_sub.m, by = c("variable"), all.x = T, suffixes = c("", ".sub"))
        pho_ptm.m$partID <- sampID2partID(sampleID_vector = as.vector(pho_ptm.m$variable), sample_map =  clinical)
        
        # maf_mut <- maf[maf$Hugo_Symbol == gene_mut & maf$geno_pos == geno_mut,]
        # maf_mut$partID <- str_split_fixed(maf_mut$Tumor_Sample_Barcode, "_", 2)[,1]
        # pho_ptm.m$mutation <- ifelse(test = (pho_ptm.m$partID %in% maf_mut$partID), yes = aa_mut, no = "")
        
        maf_mut <- maf[maf$Hugo_Symbol == gene_mut,]
        maf_mut$partID <- str_split_fixed(maf_mut$Tumor_Sample_Barcode, "_", 2)[,1]
        pho_ptm.m <- merge(pho_ptm.m, maf_mut[, c("partID", "HGVSp_Short")], by = c("partID"), all.x = T)
        colnames(pho_ptm.m)[ncol(pho_ptm.m)] <- "mutation"
        pho_ptm.m$mutation[is.na(pho_ptm.m$mutation)] <- ""
        pho_ptm.m$cancer <- cancer
        
        tab2p <- pho_ptm.m
        p = ggplot(tab2p, aes(x=value, y=value.sub, color = mutation, label= as.character(mutation)))
        p = p + geom_point(aes(shape = ifelse(mutation !=aa_mut, "b", "a")), 
                           position = position_jitterdodge(), stroke = 0, alpha = 0.8)
        p = p + geom_text_repel(aes(segment.color = mutation),
                                force = 1,
                                segment.size = 0.5, segment.alpha = 1,
                                size=2,alpha=0.8)
        p = p + labs(x=paste0(gene_ptm, " ", rsd_ptm, " phosphorylation abundance(log2 ratio"),
                     y=paste0(gene_sub, " ", rsd_sub, " phosphorylation abundance(log2 ratio"))
        p = p + theme_nogrid()
        p = p + theme(axis.title = element_text(size=6), legend.position = 'none',
                      axis.text.x = element_text(colour="black", size=8, vjust=0.5),
                      axis.text.y = element_text(colour="black", size=8))
        p = p + theme(title = element_text(size = 8))
        p = p + theme(axis.title.x = element_text(size = 10))
        p = p + theme(axis.title.y = element_text(size = 10))
        fn = paste0(makeOutDir(resultD = resultD), "cluster_", cluster, "_", gene_mut, "_", aa_mut, "-", gene_ptm, "_", rsd_ptm, "~", gene_sub, "_", rsd_sub, ".pdf")
        pdf(file = fn, height=4, width=4)
        print(p)
        dev.off()
      }
    }
  }
}

for (cancer in c("BRCA")) {
  hotpho_mut_ptm <- hotpho[hotpho[, paste0("mutation_and_ptm_in_", cancer)] & !is.na(hotpho[, paste0("mutation_and_ptm_in_", cancer)]),]
  clusters_mut_ptm <- unique(hotpho_mut_ptm$Cluster)
  
  ## input phospho level
  pho_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
  pho_head <- formatPhosphosite(phosphosite_vector = as.vector(pho_data$Phosphosite), gene_vector = as.vector(pho_data$Gene))
  pho_head$id <- paste0(pho_head$SUBSTRATE, ":", pho_head$SUB_MOD_RSD)
  pho_head$geno_pos <- pho_head_geno_pos[as.vector(pho_head$id), "geno_pos"]
  
  ## input mutation data
  maf <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/", maf_files[grepl(x = maf_files, pattern = cancer)]), data.table = F)
  maf$geno_pos <- paste0(str_split_fixed(string = maf$Chromosome, pattern = "chr", 2)[,2], ":", maf$Start_Position, "-", maf$End_Position)
  
  for (cluster in c("756.0")) {
    hotpho_cluster <- hotpho_mut_ptm[hotpho_mut_ptm$Cluster == cluster,]
    hotpho_cluster_ptm <- hotpho_cluster[hotpho_cluster$Alternate == "ptm",]
    hotpho_cluster_mut <- hotpho_cluster[hotpho_cluster$Alternate != "ptm",]
    for (i in 1:nrow(hotpho_cluster_ptm)) {
      gene_ptm <- as.character(hotpho_cluster_ptm[i, "Gene_Drug"])
      rsd_ptm <- as.character(hotpho_cluster_ptm[i, "SUB_MOD_RSD"])
      geno_ptm <- as.character(hotpho_cluster_ptm[i, "geno_pos"])
      for (j in 1:nrow(hotpho_cluster_mut)) {
        gene_mut <- as.character(hotpho_cluster_mut[j, "Gene_Drug"])
        aa_mut <- as.character(hotpho_cluster_mut[j, "Mutation_Gene"])
        geno_mut <- as.character(hotpho_cluster_mut[j, "geno_pos"])
        
        pho_ptm <- pho_data[pho_head$SUBSTRATE == gene_ptm & !is.na(pho_head$geno_pos) & pho_head$geno_pos == geno_ptm,]
        pho_ptm.m <- melt(pho_ptm)
        
        gene_sub <- "BAD"
        rsd_sub <- "S99"
        pho_sub <- pho_data[pho_head$SUBSTRATE == gene_sub & pho_head$SUB_MOD_RSD == rsd_sub,]
        pho_sub.m <- melt(pho_sub)
        
        pho_ptm.m <- merge(pho_ptm.m, pho_sub.m, by = c("variable"), all.x = T, suffixes = c("", ".sub"))
        pho_ptm.m$partID <- sampID2partID(sampleID_vector = as.vector(pho_ptm.m$variable), sample_map =  clinical)
        
        # maf_mut <- maf[maf$Hugo_Symbol == gene_mut & maf$geno_pos == geno_mut,]
        # maf_mut$partID <- str_split_fixed(maf_mut$Tumor_Sample_Barcode, "_", 2)[,1]
        # pho_ptm.m$mutation <- ifelse(test = (pho_ptm.m$partID %in% maf_mut$partID), yes = aa_mut, no = "")
        
        maf_mut <- maf[maf$Hugo_Symbol == gene_mut,]
        maf_mut$partID <- str_split_fixed(maf_mut$Tumor_Sample_Barcode, "_", 2)[,1]
        pho_ptm.m <- merge(pho_ptm.m, maf_mut[, c("partID", "HGVSp_Short")], by = c("partID"), all.x = T)
        colnames(pho_ptm.m)[ncol(pho_ptm.m)] <- "mutation"
        pho_ptm.m$mutation[is.na(pho_ptm.m$mutation)] <- ""
        pho_ptm.m$cancer <- cancer
        
        tab2p <- pho_ptm.m
        p = ggplot(tab2p, aes(x=value, y=value.sub, color = mutation, label= as.character(mutation)))
        p = p + geom_point(aes(shape = ifelse(mutation !=aa_mut, "b", "a")), 
                           position = position_jitterdodge(), stroke = 0, alpha = 0.8)
        p = p + geom_text_repel(aes(segment.color = mutation),
                                force = 1,
                                segment.size = 0.5, segment.alpha = 1,
                                size=2,alpha=0.8)
        p = p + labs(x=paste0(gene_ptm, " ", rsd_ptm, " phosphorylation abundance(log2 ratio"),
                     y=paste0(gene_sub, " ", rsd_sub, " phosphorylation abundance(log2 ratio"))
        p = p + theme_nogrid()
        p = p + theme(axis.title = element_text(size=6), legend.position = 'none',
                      axis.text.x = element_text(colour="black", size=8, vjust=0.5),
                      axis.text.y = element_text(colour="black", size=8))
        p = p + theme(title = element_text(size = 8))
        p = p + theme(axis.title.x = element_text(size = 10))
        p = p + theme(axis.title.y = element_text(size = 10))
        fn = paste0(makeOutDir(resultD = resultD), "cluster_", cluster, "_", gene_mut, "_", aa_mut, "-", gene_ptm, "_", rsd_ptm, "~", gene_sub, "_", rsd_sub, ".pdf")
        pdf(file = fn, height=4, width=4)
        print(p)
        dev.off()
      }
    }
  }
}

for (cancer in c("BRCA")) {
  hotpho_mut_ptm <- hotpho[hotpho[, paste0("mutation_and_ptm_in_", cancer)] & !is.na(hotpho[, paste0("mutation_and_ptm_in_", cancer)]),]
  clusters_mut_ptm <- unique(hotpho_mut_ptm$Cluster)
  
  ## input phospho level
  pho_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
  pho_head <- formatPhosphosite(phosphosite_vector = as.vector(pho_data$Phosphosite), gene_vector = as.vector(pho_data$Gene))
  pho_head$id <- paste0(pho_head$SUBSTRATE, ":", pho_head$SUB_MOD_RSD)
  pho_head$geno_pos <- pho_head_geno_pos[as.vector(pho_head$id), "geno_pos"]
  
  ## input mutation data
  maf <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/", maf_files[grepl(x = maf_files, pattern = cancer)]), data.table = F)
  maf$geno_pos <- paste0(str_split_fixed(string = maf$Chromosome, pattern = "chr", 2)[,2], ":", maf$Start_Position, "-", maf$End_Position)
  
  for (cluster in c("756.0")) {
    hotpho_cluster <- hotpho_mut_ptm[hotpho_mut_ptm$Cluster == cluster,]
    hotpho_cluster_ptm <- hotpho_cluster[hotpho_cluster$Alternate == "ptm",]
    hotpho_cluster_mut <- hotpho_cluster[hotpho_cluster$Alternate != "ptm",]
    for (i in 1:nrow(hotpho_cluster_ptm)) {
      gene_ptm <- as.character(hotpho_cluster_ptm[i, "Gene_Drug"])
      rsd_ptm <- as.character(hotpho_cluster_ptm[i, "SUB_MOD_RSD"])
      geno_ptm <- as.character(hotpho_cluster_ptm[i, "geno_pos"])
      for (j in 1:nrow(hotpho_cluster_mut)) {
        gene_mut <- as.character(hotpho_cluster_mut[j, "Gene_Drug"])
        aa_mut <- as.character(hotpho_cluster_mut[j, "Mutation_Gene"])
        geno_mut <- as.character(hotpho_cluster_mut[j, "geno_pos"])
        
        pho_ptm <- pho_data[pho_head$SUBSTRATE == gene_ptm & !is.na(pho_head$geno_pos) & pho_head$geno_pos == geno_ptm,]
        pho_ptm.m <- melt(pho_ptm)
        
        gene_sub <- "AKT1S1"
        rsd_sub <- "T246"
        pho_sub <- pho_data[pho_head$SUBSTRATE == gene_sub & pho_head$SUB_MOD_RSD == rsd_sub,]
        pho_sub.m <- melt(pho_sub)
        
        pho_ptm.m <- merge(pho_ptm.m, pho_sub.m, by = c("variable"), all.x = T, suffixes = c("", ".sub"))
        pho_ptm.m$partID <- sampID2partID(sampleID_vector = as.vector(pho_ptm.m$variable), sample_map =  clinical)
        pho_ptm.m$subtype <- sampID2pam50(sampleID_vector = as.vector(pho_ptm.m$variable), subtype_map = subtype_map)
        
        # maf_mut <- maf[maf$Hugo_Symbol == gene_mut & maf$geno_pos == geno_mut,]
        # maf_mut$partID <- str_split_fixed(maf_mut$Tumor_Sample_Barcode, "_", 2)[,1]
        # pho_ptm.m$mutation <- ifelse(test = (pho_ptm.m$partID %in% maf_mut$partID), yes = aa_mut, no = "")
        
        maf_mut <- maf[maf$Hugo_Symbol == gene_mut,]
        maf_mut$partID <- str_split_fixed(maf_mut$Tumor_Sample_Barcode, "_", 2)[,1]
        pho_ptm.m <- merge(pho_ptm.m, maf_mut[, c("partID", "HGVSp_Short")], by = c("partID"), all.x = T)
        colnames(pho_ptm.m)[ncol(pho_ptm.m)] <- "mutation"
        pho_ptm.m$mutation[is.na(pho_ptm.m$mutation)] <- ""
        pho_ptm.m$cancer <- cancer
        
        tab2p <- pho_ptm.m
        p = ggplot(tab2p, aes(x=value, y=value.sub))
        p = p + geom_point(aes(shape = ifelse(mutation !=aa_mut, "b", "a"), color = subtype), 
                           position = position_jitterdodge(), stroke = 0, alpha = 0.8)
        p = p + geom_text_repel(aes(segment.color = mutation, label= as.character(mutation)),
                                force = 1,
                                segment.size = 0.5, segment.alpha = 1,
                                size=2,alpha=0.8)
        p = p + labs(x=paste0(gene_ptm, " ", rsd_ptm, " phosphorylation abundance(log2 ratio"),
                     y=paste0(gene_sub, " ", c, " phosphorylation abundance(log2 ratio"))
        p = p + theme_nogrid()
        p = p + theme(axis.text.x = element_text(colour="black", size=8, vjust=0.5),
                      axis.text.y = element_text(colour="black", size=8))
        p = p + theme(title = element_text(size = 8))
        p = p + theme(axis.title.y = element_text(size = 10))
        p = p + theme(axis.title.y = element_text(size = 10))
        fn = paste0(makeOutDir(resultD = resultD), "cluster_", cluster, "_", gene_mut, "_", aa_mut, "-", gene_ptm, "_", rsd_ptm, "~", gene_sub, "_", rsd_sub, ".pdf")
        pdf(file = fn, height=4, width=6)
        print(p)
        dev.off()
        
        tab2p <- pho_ptm.m
        p = ggplot(tab2p, aes(x=cancer, y=value.sub))
        p = p + geom_point(aes(shape = ifelse(mutation !=aa_mut, "b", "a"), color = subtype), 
                           position = position_jitterdodge(), stroke = 0, alpha = 0.8)
        p = p + geom_violin(fill = "grey", color = NA, alpha = 0.2)
        p = p + geom_text_repel(aes(segment.color = mutation, label= as.character(mutation)),
                                force = 1,
                                segment.size = 0.5, segment.alpha = 1,
                                size=2,alpha=0.8)
        p = p + labs(y=paste0(gene_sub, " ", gene_sub, " phosphorylation abundance(log2 ratio"))
        p = p + theme_nogrid()
        p = p + theme(axis.text.x = element_text(colour="black", size=8, vjust=0.5),
                      axis.text.y = element_text(colour="black", size=8))
        p = p + theme(title = element_text(size = 8))
        p = p + theme(axis.title.y = element_text(size = 10))
        p = p + theme(axis.title.y = element_text(size = 10))
        fn = paste0(makeOutDir(resultD = resultD), "cluster_", cluster, "_", gene_mut, "_", aa_mut, "-", gene_sub, "_", rsd_sub, ".pdf")
        pdf(file = fn, height=4, width=6)
        print(p)
        dev.off()
      }
    }
  }
}

