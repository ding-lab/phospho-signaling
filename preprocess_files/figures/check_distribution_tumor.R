# Yige Wu @ WashU Jan 2018
# check the distribution of log ratios per sample (for deciding normalization method)
# load packages -----------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/tables/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/PHOsPHO_network/PHOsPHO_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/PHOsPHO_network/PHOsPHO_network_shared.R')


# inputs ------------------------------------------------------------------
clinical <- loadSampMap()
pam50_map <- loadPAM50Map()
msi_map <- loadMSIMap()

# process proteome and PHOsPHOproteome data ---------------------------------------------
for (cancer in cancers_sort) {
  cancer <- "OV"
  # keep shared peptide log ratio
  Pro.f = fread(paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PRO_formatted.txt",sep=""), data.table=FALSE)
  ## sometimes there are replicates for the same sample, their smaple IDs are irregular
  samples_raw <- colnames(Pro.f)
  samples <- as.vector(str_split_fixed(string = samples_raw, pattern = "\\.", 2)[,1])
  samples_tumor <- samples[samples %in% clinical$Specimen.Label[clinical$tumor_normal == "Tumor" & !is.na(clinical$Specimen.Label)]]
  samples_tumor_unique <- unique(samples_tumor)
  samples_raw_tumor <- samples_raw[samples %in% clinical$Specimen.Label[clinical$tumor_normal == "Tumor" & !is.na(clinical$Specimen.Label)]]
  
  if (cancer == "BRCA") {
    subtypemap <- data.frame(specimen_id = samples_tumor,
                             specimen_id_raw = samples_raw_tumor,
                             subtype = sampID2pam50(sampleID_vector = samples_tumor, subtype_map = pam50_map))
  }
  if (cancer == "CO") {
    subtypemap <- data.frame(specimen_id = samples_tumor, 
                             specimen_id_raw = samples_raw_tumor,
                             subtype = sampID2MSI(sampleID_vector = samples_tumor, subtype_map = msi_map, sample_map = clinical))
  }

  Pro.f <- Pro.f[, samples_raw_tumor]
  
  ## plot the distribution of un-normalized unshared log ratio
  Pro.f.m <- melt(Pro.f)
  rm(Pro.f)
  colnames(Pro.f.m) <- c('specimen_id_raw', "value")
  tab2p <- Pro.f.m
  rm(Pro.f.m)
  
  if (cancer == "OV") {
    tab2p$subtype = ""
  } else {
    subtypes <- as.vector(subtypemap$subtype)
    subtypes[is.na(subtypes)] <- "no_genomic_data"
    subtypes[subtypes == "xx"] <- "no_RNAseq"
    subtype_count <- data.frame(table(subtypes))
    rownames(subtype_count) <- as.vector(subtype_count$subtypes)
    subtypeswcount <- paste0(subtypes, "(n=", subtype_count[subtypes, "Freq"], ")")
    names(subtypeswcount) <- as.vector(subtypemap$specimen_id_raw)
    tab2p$subtype <- subtypeswcount[as.vector(tab2p$specimen_id_raw)]
    
    p <- ggplot(tab2p, aes(x=value))
    p <- p + geom_density(aes(group=specimen_id_raw, colour=subtype), alpha = 0.1)
    p <- p + theme_minimal()
    p <- p + facet_grid(subtype~., scales = "fixed", space = "fixed")
    p <- p + xlim(c(-10, 10))
    p <- p + xlab(paste0("Protein abundance detected in ", cancer))
    p
    ggsave(filename = paste0(makeOutDir(resultD = resultD), cancer, '_Pro_unnormalized_log_ratios_by_sample_by_subtype.pdf'),
           width = 5, height = 10)
    rm(p)
  }
  p <- ggplot(tab2p, aes(x=value))
  p <- p + geom_density(aes(group=specimen_id_raw, colour=subtype), alpha = 0.1)
  p <- p + theme_minimal()
  p <- p + xlim(c(-10, 10))
  p <- p + xlab(paste0("Protein abundance detected in ", cancer))
  p
  ggsave(filename = paste0(makeOutDir(resultD = resultD), cancer, '_Pro_unnormalized_log_ratios_by_sample.pdf'),
         width = 8, height = 5)
  rm(p)
  rm(tab2p)
  
  # normalize by sample
  Pro.n = fread(paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized.txt",sep=""), data.table=FALSE)
  Pro.n <- Pro.n[, samples_raw_tumor]
  ## plot the distribution of un-normalized unshared log ratio
  Pro.n.m <- melt(Pro.n)
  colnames(Pro.n.m) <- c('specimen_id_raw', "value")
  tab2p <- Pro.n.m
  rm(Pro.n.m)
  if (cancer == "OV") {
    tab2p$subtype = ""
  } else {
    subtypes <- as.vector(subtypemap$subtype)
    subtypes[is.na(subtypes)] <- "no_genomic_data"
    subtypes[subtypes == "xx"] <- "no_RNAseq"
    subtype_count <- data.frame(table(subtypes))
    rownames(subtype_count) <- as.vector(subtype_count$subtypes)
    subtypeswcount <- paste0(subtypes, "(n=", subtype_count[subtypes, "Freq"], ")")
    names(subtypeswcount) <- as.vector(subtypemap$specimen_id_raw)
    tab2p$subtype <- subtypeswcount[as.vector(tab2p$specimen_id_raw)]
    
    
    p <- ggplot(tab2p, aes(x=value))
    p <- p + geom_density(aes(group=specimen_id_raw, colour=subtype), alpha = 0.1)
    p <- p + theme_minimal()
    p <- p + facet_grid(subtype~., scales = "fixed", space = "fixed")
    p <- p + xlim(c(-10, 10))
    p <- p + xlab(paste0("Protein abundance detected in ", cancer))
    p
    ggsave(filename = paste0(makeOutDir(resultD = resultD), cancer, '_Pro_normalized_log_ratios_by_sample_by_subtype.pdf'),
           width = 5, height = 10)
    rm(p)
  }
  p <- ggplot(tab2p, aes(x=value))
  p <- p + geom_density(aes(group=specimen_id_raw, colour=subtype), alpha = 0.1)
  p <- p + theme_minimal()
  p <- p + xlim(c(-10, 10))
  p <- p + xlab(paste0("Protein abundance detected in ", cancer))
  p
  ggsave(filename = paste0(makeOutDir(resultD = resultD), cancer, '_Pro_normalized_log_ratios_by_sample.pdf'),
         width = 8, height = 5)
  rm(p)
  rm(tab2p)
  
  
  ## process PHOsPHOproteome
  PHO.f = fread(paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted.txt",sep=""), data.table=FALSE)
  PHO.f <- PHO.f[, samples_raw_tumor]
  PHO.f.m <- melt(PHO.f)
  colnames(PHO.f.m) <- c('specimen_id_raw', "value")
  rm(PHO.f)
  tab2p <- PHO.f.m
  rm(PHO.f.m)
  if (cancer == "OV") {
    tab2p$subtype = ""
  } else {
    subtypes <- as.vector(subtypemap$subtype)
    subtypes[is.na(subtypes)] <- "no_genomic_data"
    subtypes[subtypes == "xx"] <- "no_RNAseq"
    subtype_count <- data.frame(table(subtypes))
    rownames(subtype_count) <- as.vector(subtype_count$subtypes)
    subtypeswcount <- paste0(subtypes, "(n=", subtype_count[subtypes, "Freq"], ")")
    names(subtypeswcount) <- as.vector(subtypemap$specimen_id_raw)
    tab2p$subtype <- subtypeswcount[as.vector(tab2p$specimen_id_raw)]
    
    p <- ggplot(tab2p, aes(x=value))
    p <- p + geom_density(aes(group=specimen_id_raw, colour=subtype), alpha = 0.1)
    p <- p + theme_minimal()
    p <- p + facet_grid(subtype~., scales = "fixed", space = "fixed")
    p <- p + xlim(c(-10, 10))
    p <- p + xlab(paste0("Protein abundance detected in ", cancer))
    p
    ggsave(filename = paste0(makeOutDir(resultD = resultD), cancer, '_PHo_unnormalized_log_ratios_by_sample_by_subtype.pdf'),
           width = 5, height = 10)
    rm(p)
  }
  p <- ggplot(tab2p, aes(x=value))
  p <- p + geom_density(aes(group=specimen_id_raw, colour=subtype), alpha = 0.1)
  p <- p + theme_minimal()
  p <- p + xlim(c(-10, 10))
  p <- p + xlab(paste0("Protein abundance detected in ", cancer))
  p
  ggsave(filename = paste0(makeOutDir(resultD = resultD), cancer, '_PHO_unnormalized_log_ratios_by_sample.pdf'),
         width = 8, height = 5)
  rm(p)

  rm(tab2p)
  
  ## plot the distribution of un-normalized unshared log ratio
  PHO.n = fread(paste(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized.txt",sep=""), data.table=FALSE)
  PHO.n <- PHO.n[, samples_raw_tumor]
  PHO.n.m <- melt(PHO.n)
  colnames(PHO.n.m) <- c('specimen_id_raw', "value")
  rm(PHO.n)
  tab2p <- PHO.n.m
  rm(PHO.n.m)
  if (cancer == "OV") {
    tab2p$subtype = ""
  } else {
    subtypes <- as.vector(subtypemap$subtype)
    subtypes[is.na(subtypes)] <- "no_genomic_data"
    subtypes[subtypes == "xx"] <- "no_RNAseq"
    subtype_count <- data.frame(table(subtypes))
    rownames(subtype_count) <- as.vector(subtype_count$subtypes)
    subtypeswcount <- paste0(subtypes, "(n=", subtype_count[subtypes, "Freq"], ")")
    names(subtypeswcount) <- as.vector(subtypemap$specimen_id_raw)
    tab2p$subtype <- subtypeswcount[as.vector(tab2p$specimen_id_raw)]
    
    p <- ggplot(tab2p, aes(x=value))
    p <- p + geom_density(aes(group=specimen_id_raw, colour=subtype), alpha = 0.1)
    p <- p + theme_minimal()
    p <- p + facet_grid(subtype~., scales = "fixed", space = "fixed")
    p <- p + xlim(c(-10, 10))
    p <- p + xlab(paste0("Protein abundance detected in ", cancer))
    p
    ggsave(filename = paste0(makeOutDir(resultD = resultD), cancer, '_PHo_normalized_log_ratios_by_sample_by_subtype.pdf'),
           width = 5, height = 10)
    rm(p)
  }
  p <- ggplot(tab2p, aes(x=value))
  p <- p + geom_density(aes(group=specimen_id_raw, colour=subtype), alpha = 0.1)
  p <- p + theme_minimal()
  p <- p + xlim(c(-10, 10))
  p <- p + xlab(paste0("Protein abundance detected in ", cancer))
  p
  ggsave(filename = paste0(makeOutDir(resultD = resultD), cancer, '_PHO_normalized_log_ratios_by_sample.pdf'),
         width = 8, height = 5)
  rm(p)
  rm(tab2p)
}