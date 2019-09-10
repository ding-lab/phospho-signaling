# Yige Wu @ WashU 2018 Jan
## shared parameters and functions for phospho_network

# library -----------------------------------------------------------------
library(stringr)
library(reshape)
library(grid)
library(data.table)
library(readr)
library(rstudioapi)
library(ggrepel)
library(dplyr)

# source ------------------------------------------------------------------
baseD = "~/Box/"
setwd(baseD)

# folders <- strsplit(x = rstudioapi::getSourceEditorContext()$path, split = "\\/")[[1]]
# baseD <- paste0(paste0(folders[1:which(folders == "Box Sync")], collapse = "/"), "/")
# setwd(baseD)
code_top_dir <- "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/"
path2cptac_shared <- paste0(code_top_dir, "cptac2p_analysis_shared.R")
source(path2cptac_shared)

# static parameters -------------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
ppnD <- paste0(resultD, "phospho_network/")

## significance level for FDR from regression model
fdr_pp <- 0.1
fdr_pk <- 0.05

reg_sig <- c(fdr_pp, fdr_pk)
names(reg_sig) <- c("phosphatase", "kinase")

# functions ---------------------------------------------------------------
makeOutDir = function() {
  folders <- strsplit(x = rstudioapi::getSourceEditorContext()$path, split = "\\/")[[1]]
  folder_num <- which(folders == "phospho-signaling_analysis") + 1
  resultDnow <- paste(strsplit(paste(folders[folder_num:length(folders)], collapse = "/"), split = "\\.")[[1]][1], sep = "/")
  resultDnow <- paste0(resultD, resultDnow, "/")
  dir.create(resultDnow)
  resultDnow_son <- resultDnow
  dirs2make <- NULL
  while (!dir.exists(resultDnow_son)) {
    tmp <- strsplit(resultDnow_son, split = "\\/")[[1]]
    resultDnow_parent <-paste(tmp[-length(tmp)], collapse = "/")
    dir.create(resultDnow_parent)
    dir.create(resultDnow_son)
    dir.create(resultDnow)
    if (!dir.exists(resultDnow_son)) {
      dirs2make[length(dirs2make) + 1] <- resultDnow_son
    }
    resultDnow_son <- resultDnow_parent
  }
  
  if (length(dirs2make) > 0){
    for (i in 1:length(dirs2make)) {
      dir.create(dirs2make[i])
    }
  } 
  return(resultDnow)
}

sampID2partID = function(sampleID_vector, sample_map) {
  sampleIDs = sampleID_vector
  partIDs <- NULL
  for (sampID in sampleIDs) {
    partID <- as.character(unique(sample_map$Participant.ID[grepl(x = sample_map$Specimen.Label, pattern = sampID) & !is.na(sample_map$Specimen.Label) & !is.na(sample_map$Participant.ID)]))
    if (!is.null(partID)) {
      partIDs <- c(partIDs, partID)
    } else {
      partIDs <- c(partIDs, NA)
    }
  }
  return(partIDs)
}

sampID2clinicalMSI = function(sampleID_vector, sample_map, subtype_map) {
  sampleIDs = sampleID_vector
  partIDs <- sampID2partID(sampleIDs, sample_map)
  subtypes <- NULL
  for (partID in partIDs) {
    subtype <- unique(as.character(subtype_map$`Microsatellite Instability (Abnormal @ >33% loci tested)`[subtype_map$`Participant ID` == partID]))
    subtype <- subtype[!is.na(subtype)]
    if (subtype == 'Yes') {
      subtypes <- c(subtypes, "MSI")
    } else if (subtype == "No") {
      subtypes <- c(subtypes, "MSS")
    } else {
      subtypes <- c(subtypes, 'other')
    }
  }
  return(subtypes)
}

sampID2MSI = function(sampleID_vector, sample_map, subtype_map) {
  sampleIDs = sampleID_vector
  partIDs <- sampID2partID(sampleIDs, sample_map)
  subtypes <- NULL
  for (partID in partIDs) {
    subtype <- unique(as.character(subtype_map$MSI_type[subtype_map$Sample == partID]))
    subtype <- subtype[!is.na(subtype)]
    if (length(subtype) > 0) {
      subtypes <- c(subtypes, subtype)
    } else {
      subtypes <- c(subtypes, "no_genomic_data")
    }
  }
  return(subtypes)
}

sampID2pam50fromclinTab = function(sampleID_vector, subtype_map) {
  sampleIDs = sampleID_vector
  subtypes <- NULL
  for (sampID in sampleIDs) {
    subtype <- unique(as.character(subtype_map$`PAM50 Call`[grepl(x = subtype_map$`Specimen Label`, pattern = sampID)]))
    subtype <- subtype[!is.na(subtype)]
    if (length(subtype) > 0) {
      subtypes <- c(subtypes, subtype)
    } else {
      subtypes <- c(subtypes, NA)
    }
  }
  return(subtypes)
}

sampID2pam50 = function(sampleID_vector, pam50_map, sample_map) {
  sampleIDs = sampleID_vector
  partIDs <- sampID2partID(sampleIDs, sample_map)
  subtypes <- NULL
  for (partID in partIDs) {
    subtype <- unique(as.character(pam50_map$Call[pam50_map$V1 == partID]))
    subtype <- subtype[!is.na(subtype)]
    if (length(subtype) > 0) {
      subtypes <- c(subtypes, subtype)
    } else {
      subtypes <- c(subtypes, "no_RNA-seq")
    }
  }
  return(subtypes)
}


#function for normalize the variables for regression

markVadRetro = function(regression, sig_thres, enzyme_type) {
  ## Usage: input the regression result table and the significance level
  ## mark table with retrospective validated results (FDR < significance level and coefficient > 0)
  retro_table <- fread(input = "~/Box Sync/pan3can_shared_data/analysis_results/regression/table/kinase_substrate_regression_trans_edited.txt",
                       data.table = F)
  retro_table <- markSigKS(retro_table, sig_thres = sig_thres, enzyme_type = enzyme_type)
  retro_table_vad <- retro_table[retro_table$fdr_sig & retro_table$coef_sig,]
  table <- markSigKS(regression, sig_thres = sig_thres, enzyme_type = enzyme_type)
  table[, "vad_retro_pairwise"] <- (paste(table$KINASE, table$SUBSTRATE, sep = ":") %in% paste(retro_table_vad$KINASE, retro_table_vad$SUBSTRATE, sep = ":"))
  table[, "vad_retro_sitewise"] <- (table$pair %in% retro_table_vad$pair)
  return(table)
}


# evaluate regression significance ----------------------------------------
markSigKS = function(regression, sig_thres, enzyme_type) {
  ## Usage: input the regression result table and the significance level
  ## mark table with significant cis and trans pairs
  if (any(is.na(regression$enzyme_type))) {
    stop("annotate the enzyme type first!")
  }
  
  table <- regression
  table_cis <- table[table$SELF == "cis",]
  table_cis$fdr_sig <- (table_cis$FDR_pro_kin < sig_thres)
  table_trans <- table[table$SELF == "trans",]
  table_trans$fdr_sig <- (table_trans$FDR_pho_kin < sig_thres)
  if (enzyme_type == 'kinase') {
    table_cis$coef_sig <- (table_cis$coef_pro_kin > 0)
    table_trans$coef_sig <- (table_trans$coef_pho_kin > 0)
  } else if (enzyme_type == "phosphatase") {
    table_cis$coef_sig <- (table_cis$coef_pro_kin < 0)
    table_trans$coef_sig <- (table_trans$coef_pho_kin < 0)
  }
  table <- as.data.frame(rbind(table_cis, table_trans))
  table$regulated <- (table$fdr_sig & table$coef_sig)
  
  return(table)
}

markSigCan = function(regression, sig_thres, enzyme_type) {
  ## Usage: input the regression result table and the significance level
  ## mark table with whether showing significant in other cancer types
  table <- regression
  table <- markSigKS(table, sig_thres = sig_thres, enzyme_type = enzyme_type)
  table$ks_pair <- paste(table$KINASE, table$SUBSTRATE, sep = ':')
  
  sig_pairs <- vector('list', 3)
  sig_ks_pairs <- vector('list', 3)
  for (cancer in unique(table$Cancer)) {
    table_can <- table[table$Cancer == cancer & table$fdr_sig,]
    sig_pairs[[cancer]] <- as.vector(table_can$pair)
    sig_ks_pairs[[cancer]] <- as.vector(table_can$ks_pair)
  }
  for (cancer in unique(table$Cancer)) {
    table[, paste0('sig_', cancer)] <- (table$pair %in% sig_pairs[[cancer]])
    table[, paste0('ks_sig_', cancer)] <- (table$ks_pair %in% sig_ks_pairs[[cancer]])
  }
  
  table$shared3can <- (table$ks_sig_BRCA & table$ks_sig_OV & table$ks_sig_CO)
  table$uniq_BRCA <- (table$ks_sig_BRCA & !table$ks_sig_OV & !table$ks_sig_CO)
  table$uniq_OV <- (!table$ks_sig_BRCA & table$ks_sig_OV & !table$ks_sig_CO)
  table$uniq_CO <- (!table$ks_sig_BRCA & !table$ks_sig_OV & table$ks_sig_CO)
  return(table)
}


markSigSiteCan = function(regression, sig_thres, enzyme_type) {
  ## Usage: input the regression result table and the significance level
  ## mark table with whether showing significant in other cancer types
  if ('Cancer' %in% colnames(regression)) {
    table <- regression
    table$pair <- paste0(table$GENE, ":", table$SUB_GENE, ":", table$SUB_MOD_RSD)
    table <- markSigKS(table, sig_thres = sig_thres, enzyme_type = enzyme_type)
    sig_pairs <- vector('list', 3); examined_pairs <- vector('list', 3)
    for (cancer in unique(table$Cancer)) {
      table_can <- table[table$Cancer == cancer,]
      examined_pairs[[cancer]] <- as.vector(table_can$pair[!is.na(table_can$fdr_sig)])
      table_can_sig <- table[table$Cancer == cancer & table$fdr_sig & table$coef_sig,]
      sig_pairs[[cancer]] <- as.vector(table_can_sig$pair)
    }
    for (cancer in unique(table$Cancer)) {
      table[, paste0('sig_', cancer)] <- NA
      examined <- table$pair %in% examined_pairs[[cancer]]
      table[examined, paste0('sig_', cancer)] <- (table$pair[examined] %in% sig_pairs[[cancer]])
    }
    
    ## ks pairs need to be examined in all 3 cancer before marking them as shared or unique
    examined3can <- (table$pair %in% examined_pairs[["BRCA"]]) & (table$pair %in% examined_pairs[["OV"]]) & (table$pair %in% examined_pairs[["CO"]])
    table$shared3can <-  NA; table$shared3can[examined3can] <- (table$sig_BRCA[examined3can] & table$sig_OV[examined3can] & table$sig_CO[examined3can])
    table$uniq_BRCA <-  NA; table$uniq_BRCA[examined3can] <- (table$sig_BRCA[examined3can] & !table$sig_OV[examined3can] & !table$sig_CO[examined3can])
    table$uniq_OV <- NA; table$uniq_OV[examined3can] <- (!table$sig_BRCA[examined3can] & table$sig_OV[examined3can] & !table$sig_CO[examined3can])
    table$uniq_CO <-NA; table$uniq_CO[examined3can]  <- (!table$sig_BRCA[examined3can]  & !table$sig_OV[examined3can]  & table$sig_CO[examined3can] )
    return(table)
  } else {
    print('no column named Cancer!')
  }
}



formatPhosphosite = function(phosphosite_vector, gene_vector) {
  data <- phosphosite_vector
  gene <- gene_vector
  # pho_rsd_split <- data.frame(SUB_MOD_RSD = sapply(as.vector(data), FUN = function(p) strsplit(x = p, split = ":")[[1]][1]))
  pho_rsd_split <- as.data.frame(str_split_fixed(data, ":", 2))
  colnames(pho_rsd_split) <- c("transcript", "SUB_MOD_RSD")
  pho_rsd_split$SUB_MOD_RSD <- toupper(pho_rsd_split$SUB_MOD_RSD)
  pho_rsd_split$SUBSTRATE <- gene
  return(pho_rsd_split)
}


# load data ---------------------------------------------------------------
load_mut_impact_proteome <- function() {
  mut_impact_tab <- fread(input = paste0(ppnD, "genoalt/tables/adjustFDR_mut_impact_proteome/mut_impact_proteome_RNA_cptac2p_cptac3_tab_wFDR.txt"), data.table = F, sep = "\t")
  return(mut_impact_tab)
}


load_psp <- function() {
  # k_s_table = read.delim(paste(dir2cptac2_retrospective,"Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt",sep=""))
  k_s_table = read.delim(paste(dir2phospho_signaling,"resources/Phospho_databases/PhosphositePlus/May_28_2019/Kinase_Substrate_Dataset",sep=""), skip = 3)
  k_s_table <- k_s_table %>%
    filter(KIN_ORGANISM == "human" & SUB_ORGANISM == "human")
  k_s_table$pair_pro <- paste0(k_s_table$GENE, ":", k_s_table$SUB_GENE)
  k_s_table$pair <- paste0(k_s_table$pair_pro, ":", k_s_table$SUB_MOD_RSD)
  return(k_s_table)
}

load_omnipath <- function() {
  ## Usage: k_s_table <- load_ks_table({kinase/phosphatase})
  k_s_table <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
  k_s_table$pair_pro <- paste0(k_s_table$GENE, ":", k_s_table$SUB_GENE)
  k_s_table$pair <- paste0(k_s_table$pair_pro, ":", k_s_table$SUB_MOD_RSD)
  return(k_s_table)
}

omnipath_tab_new <- fread(input = "./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/resources/Phospho_databases/OmniPath/July_25_2019/ptms.txt", data.table = F)

load_es_pro_table <- function(protein) {
  omnipath_tab <- load_omnipath()
  omnipath_tab_new <- fread(input = "./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/resources/Phospho_databases/OmniPath/July_25_2019/ptms.txt", data.table = F)
  psp_tab <- load_psp()
  omnipath_tab <- omnipath_tab %>%
    filter(!(Source %in% c("PhosphoNetworks", "NetKIN", "MIMP"))) %>%
    filter(GENE %in% c(kinases, phosphatases))
  
  omnipath_tab_new <- omnipath_tab_new %>%
    filter(!(sources %in% c("PhosphoNetworks", "MIMP"))) %>%
    filter(GENE %in% c(kinases, phosphatases))
  
  omnipath_tab <- data.frame(omnipath_tab)
  
  k_s_table <- rbind(omnipath_tab %>%
                       select(GENE, SUB_GENE),
                     psp_tab %>%
                       select(GENE, SUB_GENE),
                     omnipath_tab_new %>%
                       mutate(GENE = enzyme_genesymbol, SUB_GENE = substrate_genesymbol) %>%
                       select(GENE, SUB_GENE)) %>%
    unique()
  k_s_table <- rbind(k_s_table,
                     data.frame(GENE = c("ATM", "ATR", 
                                         "MAPK1", "MAPK8", "MAPK9", "MAPK14"),
                                SUB_GENE = c("BAP1", "BAP1",
                                             "NFE2L2", "NFE2L2", "NFE2L2", "NFE2L2")))
  
  return(k_s_table)
}

load_es_pro_with_predicted_table <- function(protein) {
  omnipath_tab <- load_omnipath()
  omnipath_tab_new <- fread(input = "./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/resources/Phospho_databases/OmniPath/July_25_2019/ptms.txt", data.table = F)
  psp_tab <- load_psp()
  
  k_s_table <- rbind(omnipath_tab %>%
                       filter(GENE %in% c(kinases, phosphatases)) %>%
                       select(GENE, SUB_GENE),
                     psp_tab %>%
                       select(GENE, SUB_GENE),
                     omnipath_tab_new %>%
                       mutate(GENE = enzyme_genesymbol, SUB_GENE = substrate_genesymbol) %>%
                       filter(GENE %in% c(kinases, phosphatases)) %>%
                       select(GENE, SUB_GENE)) %>%
    unique()
  k_s_table <- rbind(k_s_table,
                     data.frame(GENE = c("ATM", "ATR", 
                                         "MAPK1", "MAPK8", "MAPK9", "MAPK14"),
                                SUB_GENE = c("BAP1", "BAP1",
                                             "NFE2L2", "NFE2L2", "NFE2L2", "NFE2L2")))
  
  return(k_s_table)
}

# list of kinase gene symbols ---------------------------------------------
omnipath_tab <- load_omnipath()
psp_tab <- load_psp()
kinases <- unique(c(as.vector(omnipath_tab$GENE[omnipath_tab$enzyme_type == "kinase"]), as.vector(psp_tab$GENE)))
phosphatases <- unique(omnipath_tab$GENE[omnipath_tab$enzyme_type == "phosphatase"])
## clean up
kinases <- kinases[!(kinases %in% c("APC" ,"AXIN1", "CTNNB1", "CCNE1", "CCNE2", "CCNA2", "RPTOR", "RICTOR", "AKT1S1", "DEPTOR", "IL6ST", "PRKAR2B", 
                                    "CCND1", "PIK3R1", "CDKN2A", "CCNB1"))]


# list of tyrosine kinases ------------------------------------------------
library(dplyr)
psp_tab$aa <- substr(x = psp_tab$SUB_MOD_RSD, start = 1, stop = 1)
psp_tab <- data.frame(psp_tab)

psp_aa_long <- psp_tab %>% 
  filter(KIN_ORGANISM == "human" & SUB_ORGANISM == "human") %>%
  dplyr::select(GENE, aa) %>%
  unique() %>%
  table() %>%
  as.data.frame()

psp_aa_wide <- reshape2::dcast(psp_aa_long, GENE ~ aa)
psp_aa_wide <- data.frame(psp_aa_wide)
psp_aa_wide %>%
  filter(Y == 1) %>%
  nrow()

TK_list <- psp_aa_wide %>%
  filter(Y >= 1)

TK_list <- as.vector(TK_list$GENE)


# Annotate regression table -----------------------------------------------
annotate_enzyme_type  <- function(regression, kinases, phosphatases) {
  table2annotate <- regression
  if (!("GENE" %in% colnames(table2annotate))) {
    table2annotate$GENE <- table2annotate$KINASE
  }
  
  table2annotate$enzyme_type <- NULL
  
  table2annotate$enzyme_type <- ""
  table2annotate$enzyme_type[table2annotate$GENE %in% kinases] <- "kinase"
  table2annotate$enzyme_type[table2annotate$GENE %in% phosphatases] <- "phosphatase"
  return(table2annotate)
}

pairs_cis_direct <- c(psp_tab$pair[psp_tab$pair_pro == "ERBB2:ERBB2"], 
                      "ERBB2:ERBB2:Y1109", "ERBB2:ERBB2:Y1218", 
                      psp_tab$pair[psp_tab$pair_pro == "AKT1:AKT1"], 
                      "AKT1:AKT1:S184", "AKT1:AKT1:T10", "AKT1:AKT1:Y315", "AKT1:AKT1:S129", "AKT1:AKT1:T308", "AKT1:AKT1:T450",
                      psp_tab$pair[psp_tab$pair_pro == "EGFR:EGFR"],
                      psp_tab$pair[psp_tab$pair_pro == "MTOR:MTOR"], 
                      "MTOR:MTOR:S2448", "MTOR:MTOR:S2478S2481", "MTOR:MTOR:S2481",
                      psp_tab$pair[psp_tab$pair_pro == "MAPK1:MAPK1"],
                      "MAPK1:MAPK1:T185", "MAPK1:MAPK1:Y187",
                      psp_tab$pair[psp_tab$pair_pro == "MAP2K1:MAP2K1"],
                      "MAP2K1:MAP2K1:S218", "MAP2K1:MAP2K1:S218S222", "MAP2K1:MAP2K1:S222", "MAP2K1:MAP2K1:S286", "MAP2K1:MAP2K1:S212",
                      psp_tab$pair[psp_tab$pair_pro == "BRAF:BRAF"],
                      "BRAF:BRAF:S151", "BRAF:BRAF:S365", "BRAF:BRAF:S446", "BRAF:BRAF:S750", "BRAF:BRAF:S750T753", "BRAF:BRAF:T401",
                      psp_tab$pair[psp_tab$pair_pro == "FGFR2:FGFR2"],
                      psp_tab$pair[psp_tab$pair_pro == "INSR:INSR"], "INSR:INSR:S1354",
                      psp_tab$pair[psp_tab$pair_pro == "MET:MET"], 
                      psp_tab$pair[psp_tab$pair_pro == "PDGFRA:PDGFRA"], 
                      psp_tab$pair[psp_tab$pair_pro == "PIK3CG:PIK3CG"])

annotate_ks_source <- function(regression) {
  table2annotate <- regression
  table2annotate$Source <- NULL
  table2annotate$is.direct <- NULL
  
  omnipath_tab <- load_omnipath()
  psp_tab <- load_psp()
  
  source_tab <- data.frame(pair = psp_tab$pair[!(psp_tab$pair %in% omnipath_tab$pair)], Source = "PhosphoSite")
  source_tab <- rbind(source_tab, omnipath_tab[, c("pair", "Source")])
  source_tab <- rbind(source_tab, data.frame(pair = pairs_cis_direct, Source = "PhosphoSite"))
  source_tab <- unique(source_tab)
  table2annotate <- merge(table2annotate, source_tab, by = c("pair"), all.x = T)
  table2annotate$is.direct <- ifelse(!is.na(table2annotate$Source) & !(table2annotate$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP")), T, F)
  return(table2annotate)
}


# filter regression and adjust results ------------------------------------
change_regression_nonNA <- function(regression, reg_nonNA, reg_sig) {
  table2filter <- regression
  table2filter$enzyme_type <- NULL
  table2filter <- annotate_enzyme_type(regression = table2filter, kinases = kinases, phosphatases = phosphatases)
  table2filter <- table2filter[table2filter$Size >= reg_nonNA,]
  ## clean up
  table2filter <- table2filter %>%
    filter(enzyme_type != "")
  
  ## adjust p-values to FDR
  name = c("pro_kin","pro_sub","pho_kin")
  
  for (cancer_tmp in unique(table2filter$Cancer)) {
    for (enzyme_type_tmp in unique(table2filter$enzyme_type)) {
      for(self_tmp in unique(table2filter$SELF)) {
        for(coln_tmp in name) {#adjust pvalues for each variable
          row <- (table2filter$self==self_tmp) & (table2filter$Cancer==cancer_tmp) & (table2filter$enzyme_type==enzyme_type_tmp)
          table2filter[row, paste("FDR_",coln_tmp,sep = "")] <- p.adjust(table2filter[row,paste("P_",coln_tmp,sep = "")], method = "fdr")
        }
      }
    }
  }
  
  table2mark <- NULL
  ## mark significance
  for (enzyme_type_tmp in unique(table2filter$enzyme_type)) {
    tab_tmp <- table2filter[table2filter$enzyme_type == enzyme_type_tmp,]
    tab_tmp <- markSigKS(regression = tab_tmp, sig_thres = reg_sig[enzyme_type_tmp], enzyme_type = enzyme_type_tmp)
    table2mark <- rbind(table2mark, tab_tmp)
  }
  
  return(table2mark)
}

adjust_regression_by_nonNA <- function(regression, reg_nonNA, reg_sig) {
  table2filter <- regression
  table2filter$enzyme_type <- NULL
  table2filter <- annotate_enzyme_type(regression = table2filter, kinases = kinases, phosphatases = phosphatases)
  table2filter <- table2filter[table2filter$Size >= reg_nonNA,]
  ## clean up
  table2filter <- table2filter %>%
    filter(enzyme_type != "")
  
  ## adjust p-values to FDR
  name = c("pro_kin","pro_sub","pho_kin")
  
  for(coln_tmp in name) {#adjust pvalues for each variable
    table2filter[, paste("FDR_",coln_tmp,sep = "")] <- FDR_by_id_columns(p_vector = table2filter[,paste("P_",coln_tmp,sep = "")], id_columns = c("SELF", "Cancer", "enzyme_type"), df = table2filter)
  }
  
  table2mark <- NULL
  ## mark significance
  for (enzyme_type_tmp in unique(table2filter$enzyme_type)) {
    tab_tmp <- table2filter[table2filter$enzyme_type == enzyme_type_tmp,]
    tab_tmp <- markSigKS(regression = tab_tmp, sig_thres = reg_sig[enzyme_type_tmp], enzyme_type = enzyme_type_tmp)
    table2mark <- rbind(table2mark, tab_tmp)
  }
  
  return(table2mark)
}

omnipath_tab <- annotate_ks_source(regression = omnipath_tab)
