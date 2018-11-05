# Yige Wu @ WashU 2018 July
## check the occurance of hotpho cluster mutation and phosphosite in the cptac prospective dataset

# source ------------------------------------------------------------------
source("Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R")
library(readr)

# set variables -----------------------------------------------------------
cancers_sort <- c("BRCA", "OV", "CO")

# inputs hotpho results (hg38 based)------------------------------------------------------------------
# hotpho <- read_delim("~/Box Sync/cptac2p/hotpho/Data_201807_cc.p0.05.cluster_transcriptSynced.tsv", 
#                      "\t", escape_double = FALSE, col_types = cols(Cluster = col_character()), trim_ws = TRUE)
hotpho <- read_delim("~/Box Sync/Ding_Lab/Projects_Current/hotpho_data/output/Data_201807_cc.p0.05.cluster_transcriptSynced_hybrid_SourceAnnotated.tsv", 
                     "\t", escape_double = FALSE, col_types = cols(Cluster = col_character()), trim_ws = TRUE)
hotpho$geno_pos <- paste0(hotpho$Chromosome, ":", hotpho$Start, "-", hotpho$Stop)
tmp <- as.vector(hotpho$GenomicPosition)
tmp <- str_split_fixed(tmp, pattern = 'chr|:g.|_|/c', 5)
hotpho$geno_pos[hotpho$Alternate == "ptm"] <- paste0(tmp[hotpho$Alternate == "ptm",2], ":", tmp[hotpho$Alternate == "ptm",3], "-", tmp[hotpho$Alternate == "ptm",4])
hotpho$SUB_MOD_RSD <- str_split_fixed(string = hotpho$Mutation_Gene, pattern = "p.", 2)[,2]
hotpho_hybrid <- hotpho[hotpho$Type == "Hybrid",]


# input sample mapping file -----------------------------------------------
clinical <- read.delim(paste0(cptac_sharedD, "Specimen_Data_20161005_Yige_20180307.txt"), sep = "\t")
clinical <- data.frame(clinical)

# input phosphosite genomic coordinates (hg38) -----------------------------------
pho_head_geno_pos <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid.txt"), data.table = F)
# pho_head_geno_pos <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid_hg38_hg19.txt"), data.table = F)
pho_head_geno_pos$id <- paste0(pho_head_geno_pos$SUBSTRATE, ":", pho_head_geno_pos$SUB_MOD_RSD)
pho_head_geno_pos <- pho_head_geno_pos[!duplicated(pho_head_geno_pos$id),]
rownames(pho_head_geno_pos) <- pho_head_geno_pos$id
pho_head_geno_pos$geno_pos <- paste0(pho_head_geno_pos$seq_region_name, ":", pho_head_geno_pos$start, "-", pho_head_geno_pos$end)


# check mutations ---------------------------------------------------------
## go from each cluster whoese mutated gene is found in the maf
## check the exact mutation to narrow down the hotpho mutations present in CPTAC data
## a list for recording the patient IDs and sample IDs with hotpho mutations
genes_hotpho_maf <- list()
cptac_samples <- list()
for (cancer in cancers_sort){
  ## a column for numbering CPTAC data availibility
  num_cptac_samples <- vector(mode = "numeric", length = nrow(hotpho))

  ## a column for logging CPTAC data availibility
  notes_cptac_samples <- vector(mode = "character", length = nrow(hotpho))
  
  cptac_samples[[cancer]] <- list()

  maf <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/", maf_files[grepl(x = maf_files, pattern = cancer)]), data.table = F)
  maf$geno_pos <- paste0(str_split_fixed(string = maf$Chromosome, pattern = "chr", 2)[,2], ":", maf$Start_Position, "-", maf$End_Position)
  
  genes_hotpho_maf[[cancer]] <- intersect(unique(hotpho$Gene_Drug[hotpho$Type == "Hybrid" & hotpho$Alternate != "ptm"]), unique(maf$Hugo_Symbol))
  
  ## get the cluster with mutated gene name in CPTAC MAF
  clusters_can <- unique(hotpho_hybrid$Cluster[hotpho_hybrid$Gene_Drug %in% genes_hotpho_maf[[cancer]]])
  notes_cptac_samples[hotpho$Type == "Hybrid" & !(hotpho$Cluster %in% clusters_can) & hotpho$Alternate != "ptm"] <- "maf_not_match_gene_name"
  
  ## 342 clusters
  for (cluster in clusters_can) {
    hotpho_cluster <- hotpho[hotpho$Cluster == cluster,]
    cptac_samples[[cancer]][[cluster]] <- list()
    cptac_samples[[cancer]][[cluster]][["mutation"]] <- list()
    ## extract only mutations
    hotpho_cluster_muts <- hotpho_cluster[hotpho_cluster$Alternate != "ptm",]
    for (i in 1:nrow(hotpho_cluster_muts)) {
      gene_hotpho <- as.character(hotpho_cluster_muts[i, "Gene_Drug"])
      hgvsp_hotpho <- as.character(hotpho_cluster_muts[i, "Mutation_Gene"])
      transcript_hotpho <- as.character(hotpho_cluster_muts[i, "Transcript"])
      geno_pos_hotpho <- as.character(hotpho_cluster_muts[i, "geno_pos"])
      
      row_hotpho <- (hotpho$Cluster == cluster & hotpho$Gene_Drug == gene_hotpho & hotpho$Mutation_Gene == hgvsp_hotpho)
      maf_gene <- maf[maf$Hugo_Symbol == gene_hotpho,]
      
      if (nrow(maf_gene) > 0) {
        maf_gene_pos <- maf_gene[maf_gene$geno_pos == geno_pos_hotpho,]
        if (nrow(maf_gene_pos) > 0) {
          num_cptac_samples[row_hotpho] <- nrow(maf_gene_pos)
          cptac_samples[[cancer]][[cluster]][["mutation"]][[paste0(gene_hotpho, ":", hgvsp_hotpho)]] <- list()
          cptac_samples[[cancer]][[cluster]][["mutation"]][[paste0(gene_hotpho, ":", hgvsp_hotpho)]][["Gene_Drug"]] <- gene_hotpho
          cptac_samples[[cancer]][[cluster]][["mutation"]][[paste0(gene_hotpho, ":", hgvsp_hotpho)]][["Mutation_Gene"]] <- hgvsp_hotpho
          cptac_samples[[cancer]][[cluster]][["mutation"]][[paste0(gene_hotpho, ":", hgvsp_hotpho)]][["partIDs"]] <- str_split_fixed(maf_gene_pos$Tumor_Sample_Barcode, "_", 2)[,1]
          if (transcript_hotpho == unique(maf_gene_pos$Transcript_ID)) {
            notes_cptac_samples[row_hotpho] <- "maf_match_genomic_coordinates_and_Transcript_ID"
          } else {
            notes_cptac_samples[row_hotpho] <- "maf_match_genomic_coordinates_but_not_Transcript_ID"
          }
        } else {
          notes_cptac_samples[row_hotpho] <- "maf_match_gene_but_not_genomic_coordinates"
        }
      } else {
        notes_cptac_samples[row_hotpho] <- "maf_not_match_gene_name"
      }
    }
  }
  hotpho[, paste0("num_cptac_samples_", cancer)] <- num_cptac_samples
  hotpho[, paste0("notes_cptac_samples_", cancer)] <- notes_cptac_samples
}

# check ptm ---------------------------------------------------------------
## input genomic coordinates of phosphosites
for (cancer in cancers_sort){
  ## a column for numbering CPTAC data availibility
  num_cptac_samples <- unlist(as.data.frame(hotpho[, paste0("num_cptac_samples_", cancer)]))
  
  ## a column for logging CPTAC data availibility
  notes_cptac_samples <- unlist(as.data.frame(hotpho[, paste0("notes_cptac_samples_", cancer)]))
  
  pho <- fread(input = paste0(cptac_sharedD, cancer, '/', prefix[cancer], '_PHO_formatted_normalized_noControl.txt'), data.table = F)
  pho_head <- formatPhosphosite(phosphosite_vector = as.vector(pho$Phosphosite), gene_vector = as.vector(pho$Gene))
  pho_head$id <- paste0(pho_head$SUBSTRATE, ":", pho_head$SUB_MOD_RSD)
  pho_head$geno_pos <- pho_head_geno_pos[as.vector(pho_head$id), "geno_pos"]
  ## get the cluster with mutated gene name in CPTAC MAF
  clusters_can <- unique(hotpho_hybrid$Cluster[hotpho_hybrid$Gene_Drug %in% pho$Gene & hotpho_hybrid$Alternate == "ptm"])
  notes_cptac_samples[hotpho$Type == "Hybrid" & !(hotpho$Cluster %in% clusters_can) & hotpho$Alternate == "ptm"] <- "phospho_not_match_gene_name"
  
  for (cluster in clusters_can) {
    hotpho_cluster <- hotpho[hotpho$Cluster == cluster,]
    if (!(cluster %in% names(cptac_samples[[cancer]]))) {
      cptac_samples[[cancer]][[cluster]] <- list()
      cptac_samples[[cancer]][[cluster]][["ptm"]] <- list()
    }
    
    ## extract only mutations
    hotpho_cluster_ptms <- hotpho_cluster[hotpho_cluster$Alternate == "ptm",]
    for (i in 1:nrow(hotpho_cluster_ptms)) {
      gene_hotpho <- as.character(hotpho_cluster_ptms[i, "Gene_Drug"])
      rsd_hotpho <- as.character(hotpho_cluster_ptms[i, "SUB_MOD_RSD"])
      hgvsp_hotpho <- as.character(hotpho_cluster_ptms[i, "Mutation_Gene"])
      geno_pos_hotpho <- as.character(hotpho_cluster_ptms[i, "geno_pos"])
      
      row_hotpho <- (hotpho$Cluster == cluster & hotpho$Gene_Drug == gene_hotpho & hotpho$Mutation_Gene == hgvsp_hotpho)
      pho_gene <- pho[pho$Gene == gene_hotpho,]
      
      if (nrow(pho_gene) > 0) {
        pho_gene_pos <- pho[pho$Gene == gene_hotpho & !is.na(pho_head$geno_pos) & pho_head$geno_pos == geno_pos_hotpho, !(colnames(pho) %in% c("Gene", "Phosphosite"))]
        pho_gene_pos <- pho_gene_pos[, !is.na(pho_gene_pos)]
        if (ncol(pho_gene_pos) > 0) {
          num_cptac_samples[row_hotpho] <- ncol(pho_gene_pos)
          cptac_samples[[cancer]][[cluster]][["ptm"]][[paste0(gene_hotpho, ":", hgvsp_hotpho)]] <- list()
          cptac_samples[[cancer]][[cluster]][["ptm"]][[paste0(gene_hotpho, ":", hgvsp_hotpho)]][["Gene_Drug"]] <- gene_hotpho
          cptac_samples[[cancer]][[cluster]][["ptm"]][[paste0(gene_hotpho, ":", hgvsp_hotpho)]][["Mutation_Gene"]] <- hgvsp_hotpho
          cptac_samples[[cancer]][[cluster]][["ptm"]][[paste0(gene_hotpho, ":", hgvsp_hotpho)]][["sampIDs"]] <- colnames(pho_gene_pos)
          cptac_samples[[cancer]][[cluster]][["ptm"]][[paste0(gene_hotpho, ":", hgvsp_hotpho)]][["partIDs"]] <- sampID2partID(sampleID_vector = colnames(pho_gene_pos), sample_map = clinical)
          notes_cptac_samples[row_hotpho] <- "phospho_match_genomic_coordinates"
        } else {
          notes_cptac_samples[row_hotpho] <- "phospho_match_gene_but_not_genomic_coordinates"
        }
      } else {
        notes_cptac_samples[row_hotpho] <- "phospho_not_match_gene_name"
      }
    }
  }
  hotpho[, paste0("num_cptac_samples_", cancer)] <- num_cptac_samples
  hotpho[, paste0("notes_cptac_samples_", cancer)] <- notes_cptac_samples
}

# check ptms & mutations ---------------------------------------------------------
for (cancer in c("BRCA", "OV", "CO")) {
  clusters_can <- intersect(unique(hotpho$Cluster[hotpho[, paste0("num_cptac_samples_", cancer)] > 0 & hotpho$Alternate != "ptm"]),
                            unique(hotpho$Cluster[hotpho[, paste0("num_cptac_samples_", cancer)] > 0 & hotpho$Alternate == "ptm"]))
  hotpho[, paste0("mutation_and_ptm_in_", cancer)] <- FALSE
  
  for (cluster in clusters_can) {
    hotpho_cluster <- hotpho[hotpho$Cluster == cluster,]
    hotpho_cluster_ptms <- hotpho_cluster[hotpho_cluster$Alternate == "ptm",]
    hotpho_cluster_muts <- hotpho_cluster[hotpho_cluster$Alternate != "ptm",]
    partID_mut <- cptac_samples[[cancer]][[cluster]][["mutation-ptm"]] <- list()
    
    for (i in 1:nrow(hotpho_cluster_ptms)) {
      gene_ptm <- as.character(hotpho_cluster_ptms[i, "Gene_Drug"])
      hgvsp_ptm <- as.character(hotpho_cluster_ptms[i, "Mutation_Gene"])
      partID_ptm <- cptac_samples[[cancer]][[cluster]][["ptm"]][[paste0(gene_ptm, ":", hgvsp_ptm)]][["partIDs"]]
      for (j in 1:nrow(hotpho_cluster_muts)) {
        gene_mut <- as.character(hotpho_cluster_muts[j, "Gene_Drug"])
        hgvsp_mut <- as.character(hotpho_cluster_muts[j, "Mutation_Gene"])
        partID_mut <- cptac_samples[[cancer]][[cluster]][["mutation"]][[paste0(gene_mut, ":", hgvsp_mut)]][["partIDs"]]
        
        partID_ptm_mut <- intersect(partID_ptm, partID_mut)
        if (length(partID_ptm_mut) > 0) {
          cptac_samples[[cancer]][[cluster]][["mutation-ptm"]][[paste0(gene_mut, ":", hgvsp_mut, "-", gene_ptm, ":", hgvsp_ptm)]] <- list()
          cptac_samples[[cancer]][[cluster]][["mutation-ptm"]][[paste0(gene_mut, ":", hgvsp_mut, "-", gene_ptm, ":", hgvsp_ptm)]][["gene_ptm"]] <- gene_ptm
          cptac_samples[[cancer]][[cluster]][["mutation-ptm"]][[paste0(gene_mut, ":", hgvsp_mut, "-", gene_ptm, ":", hgvsp_ptm)]][["hgvsp_ptm"]] <- gene_ptm
          cptac_samples[[cancer]][[cluster]][["mutation-ptm"]][[paste0(gene_mut, ":", hgvsp_mut, "-", gene_ptm, ":", hgvsp_ptm)]][["gene_mut"]] <- gene_mut
          cptac_samples[[cancer]][[cluster]][["mutation-ptm"]][[paste0(gene_mut, ":", hgvsp_mut, "-", gene_ptm, ":", hgvsp_ptm)]][["hgvsp_mut"]] <- hgvsp_mut
          cptac_samples[[cancer]][[cluster]][["mutation-ptm"]][[paste0(gene_mut, ":", hgvsp_mut, "-", gene_ptm, ":", hgvsp_ptm)]][["partIDs"]] <- partID_ptm_mut
          
          hotpho[hotpho$Cluster == cluster &  ( (hotpho$Mutation_Gene == hgvsp_mut & hotpho$Gene_Drug == gene_mut ) | (hotpho$Mutation_Gene == hgvsp_ptm & hotpho$Gene_Drug == gene_ptm)), paste0("mutation_and_ptm_in_", cancer)] <- TRUE
        }
      }
    }
  }
}
write.table(x = hotpho, file = paste0(makeOutDir(resultD = resultD), 
                                      "Data_201807_cc.p0.05.cluster_transcriptSynced_CPTAC2p_MAF_availability_annotated.tsv"),
            row.names = F, quote = F, sep = "\t")
saveRDS(object = cptac_samples, file = paste0(makeOutDir(resultD = resultD), "hotpho_hybrid_clusters_in_cptac_samples.RDS"))
