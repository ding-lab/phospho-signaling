# Yige Wu @ WashU 2018 Jan
## shared parameters and functions for phospho_network

# library -----------------------------------------------------------------
library(stringr)
library(reshape)
library(grid)
require(plyr)
library(data.table)
library(WGCNA)


# static parameters -------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
inputD <- "/Users/yigewu/Box Sync/cptac2p/cptac_shared/"
resultD <- paste0(inputD, "analysis_results/phospho_network/")
prefix <- NULL
prefix["BRCA"] <- "BRCA_BI"
prefix["OV"] <- "OV_PNNL"
prefix["CO"] <- "CO_PNNL"
sig <- 0.05 # significance level


# functions ---------------------------------------------------------------
setwd(baseD)
makeOutDir = function(...) {
  resultDnow <- paste0(gsub('cptac2p_analysis/','cptac2p/cptac_shared/analysis_results/', 
                            str_split_fixed(dirname(rstudioapi::getSourceEditorContext()$path), baseD, 2)[,2]),
                       '/')
  system(paste0('mkdir ', resultDnow))
  return(resultDnow)
}

makeOutDir = function(...) {
  folders <- strsplit(x = rstudioapi::getSourceEditorContext()$path, split = "\\/")[[1]]
  folder_num <- which(folders == "cptac2p_analysis") + 1
  resultDnow <- paste("cptac2p/cptac_shared/analysis_results", paste(folders[folder_num:length(folders)], collapse = "/"), sep = "/")
  resultDnow <- paste0("./", strsplit(resultDnow, split = "\\.")[[1]][1])
  system(paste0('mkdir ', resultDnow))
  return(resultDnow)
}

bubble_heatmap = function(cancer, protein, self, id, table_bubble, table_substrate, table_kinase, outDir) {
  ## Usage: plut bubble plot and heatmap side by side
  ## input:
  ##        identifier(string): cancer type of the dataset, kinase/phophotase substrate pairs; cis/trans; additional ID
  ##        table with columns: "KINASE"      "SUBSTRATE"   "SUB_MOD_RSD" "fdr"         "coef"        "Cancer"      "pair"        "sig"   
  ##        subtype mean of each phosphosite, protein in input table
  pho_rsd_split <- formatPhosphosite(table_substrate$Phosphosite, table_substrate$Gene)
  cohorts <- colnames(table_substrate)[!(colnames(table_substrate) %in% c('Gene', 'Phosphosite'))]
  
  lim = max(abs(max(table_bubble$coef, na.rm = T)),abs(min(table_bubble$coef, na.rm = T)))
  p = ggplot(table_bubble,aes(x='Assoc', y=pair))# make this the original ethni
  p = p + facet_grid(KINASE~., drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
  p = p + geom_point(aes(color=coef, size =-log10(fdr)),pch=16) 
  p = p + scale_color_gradientn(name= "Coefficient", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
  p = p + theme_bw() + theme_nogrid()
  p = p + labs(x = "", y="kinase:substrate:SUB_MOD_RSD")
  p = p + theme(axis.title = element_blank(), axis.text.y = element_text(colour="black", size=10))
  p = p + theme(axis.text.x = element_text(angle = -30, vjust = -0.1, hjust = 0))
  p = p + theme(axis.ticks=element_blank(),legend.position="bottom")
  p = p + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
#  p
  # fn = paste(baseD,'pan3can_shared_data/analysis_results/bubble_charts/',protein,'/',c,'_sorted_cis_',self,'_fdr',sig,'_top_',top,'.pdf',sep ="")
  # ggsave(file=fn, height=10, width=6, useDingbats=FALSE)
  
  # extract the subtype phosphorylation level according to bubble ch --------
  pho_table <- c()
  for (cohort in cohorts) {
    temp <- table_bubble[,c("KINASE","SUBSTRATE","SUB_MOD_RSD","Cancer")]
    temp$pair <- paste(temp$SUBSTRATE,temp$SUB_MOD_RSD,sep = ":")
    temp$cohort <- cohort
    temp$pho_subtype <- NA
    for (i in 1:nrow(table_bubble)) {
      pho_temp <- table_substrate[pho_rsd_split$SUBSTRATE==as.character(temp$SUBSTRATE[i]) & pho_rsd_split$SUB_MOD_RSD==as.character(temp$SUB_MOD_RSD[i]),cohort]
      if (length(pho_temp) > 0) {
        temp$pho_subtype[i] <- pho_temp
      }
    }
    pho_table <- rbind(pho_table,temp)
  }
  
  # plot heatmap for subtype phosphosite level --------------------------
  lim = max(abs(max(pho_table$pho_subtype, na.rm = T)),abs(min(pho_table$pho_subtype, na.rm = T)))
  cap <- min(2, lim)
  pho_table$pho_capped <- pho_table$pho_subtype
  pho_table$pho_capped[pho_table$pho_capped > cap] <- cap
  pho_table$pho_capped[pho_table$pho_capped < (-cap)] <- (-cap)
  plot2 = ggplot(pho_table)
  plot2 = plot2 + geom_tile(aes(x=cohort, y=pair, fill=pho_capped), color=NA)#, linetype="blank") 
  plot2 = plot2 + facet_grid(KINASE~., drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
  plot2 = plot2 + scale_fill_gradientn(name= "phospho_level", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
  plot2 = plot2 + theme_bw() + theme_nogrid()
  plot2 = plot2 + theme(axis.text.x = element_text(angle = -30, vjust = 0, hjust = 0))
  plot2 = plot2 + theme(axis.title = element_blank(),  axis.text.y = element_blank())
  plot2 = plot2 + theme(axis.ticks=element_blank(),legend.position="bottom")
  plot2 = plot2 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
  
  # extract the subtype protein level according to bubble ch --------
  pro_table <- c()
  for (cohort in cohorts) {
    temp <- table_bubble[,c("KINASE","SUBSTRATE","SUB_MOD_RSD","Cancer","pair")]
    temp$cohort <- cohort
    temp$pro_subtype <- NA
    for (i in 1:nrow(table_bubble)) {
      temp$pro_subtype[i] <- table_kinase[table_kinase$Gene == as.character(temp$KINASE[i]) ,cohort]
    }
    pro_table <- rbind(pro_table,temp)
  }
  
  # plot heatmap for subtype kinase expression level --------------------------
  lim = max(abs(max(pro_table$pro_subtype)),abs(min(pro_table$pro_subtype)))
  cap <- min(2, lim)
  pro_table$pro_capped <- pro_table$pro_subtype
  pro_table$pro_capped[pro_table$pro_capped > cap] <- cap
  pro_table$pro_capped[pro_table$pro_capped < (-cap)] <- (-cap)
  plot3 = ggplot(pro_table)
  plot3 = plot3 + geom_tile(aes(x=cohort, y=pair, fill=pro_capped), color=NA)#, linetype="blank") 
  plot3 = plot3 + facet_grid(KINASE~., drop=T, space = "free_y",scales = "free_y")#, space = "free", scales = "free")
  plot3 = plot3 + scale_fill_gradientn(name= "pro_level", na.value=NA, colours=RdBu1024, limit=c(-lim,lim))
  plot3 = plot3 + theme_bw() + theme_nogrid()
  plot3 = plot3 + theme(axis.title = element_blank(), axis.text.y = element_blank())
  plot3 = plot3 + theme(axis.text.x = element_text(angle = -30, vjust = 0, hjust = 0))
  plot3 = plot3 + theme(axis.ticks=element_blank(), legend.position="bottom")
  plot3 = plot3 + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
  
  # plot together --------------------------------------------------------------------
  fn = paste(outDir, "cptac2p_", cancer,"_", protein,'_regrssion_',self, "_", id, '_heatmapWsubtype.pdf',sep ="")
  grid.newpage()
  # pdf(fn, height = 14*(len+10)/(top+10), width = 8*(len+10)/(top+10))
  pdf(fn, height = 12*(len+10)/(top+10), width = 5.5)
  
  plot123 <- cbind(ggplotGrob(p), ggplotGrob(plot2), ggplotGrob(plot3), size = "last")
  title <- textGrob(paste0(cancer, " ", self,"-regulated ",protein, "-substrate pairs(", id ,")"),gp=gpar(fontsize=16))
  padding <- unit(5,"mm")
  plottgt <- gtable_add_rows(plot123, 
                             heights = grobHeight(title) + padding,
                             pos = 0)
  plottgt <- gtable_add_grob(plottgt, title, 1, 1, 1, ncol(plottgt))
  grid.draw(plottgt)
  dev.off()
}

sampID2partID = function(sampleID_vector, sample_map) {
  sampleIDs = sampleID_vector
  partIDs <- NULL
  for (sampID in sampleIDs) {
    partID <- as.character(sample_map$Participant.ID[sample_map$Specimen.Label == sampID])
    if (!is.null(partID)) {
      partIDs <- c(partIDs, partID)
    } else {
      partIDs <- c(partIDs, NA)
    }
  }
  return(partIDs)
}

sampID2msi = function(sampleID_vector, sample_map, subtype_map) {
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

sampID2pam50 = function(sampleID_vector, subtype_map) {
  sampleIDs = sampleID_vector
  subtypes <- NULL
  for (sampID in sampleIDs) {
    subtype <- unique(as.character(subtype_map$`PAM50 Call`[subtype_map$`Specimen Label` == sampID]))
    subtype <- subtype[!is.na(subtype)]
    if (length(subtype) > 0) {
      subtypes <- c(subtypes, subtype)
    } else {
      subtypes <- c(subtypes, NA)
    }
  }
  return(subtypes)
}

#function for normalize the variables for regression
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

markSigKS = function(regression, sig_thres) {
  ## Usage: input the regression result table and the significance level
  ## mark table with significant cis and trans pairs
  table <- regression
  table_cis <- table[table$SELF == "cis",]
  table_cis$fdr_sig <- (table_cis$FDR_pro_kin < sig_thres)
  table_cis$coef_pos <- (table_cis$coef_pro_kin > 0)
  table_trans <- table[table$SELF == "trans",]
  table_trans$fdr_sig <- (table_trans$FDR_pho_kin < sig_thres)
  table_trans$coef_pos <- (table_trans$coef_pho_kin > 0)
  table <- as.data.frame(rbind(table_cis, table_trans))
  return(table)
}

markSigCan = function(regression, sig_thres) {
  ## Usage: input the regression result table and the significance level
  ## mark table with whether showing significant in other cancer types
  table <- regression
  table <- markSigKS(table, sig_thres = sig_thres)
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


markSigSiteCan = function(regression, sig_thres) {
  ## Usage: input the regression result table and the significance level
  ## mark table with whether showing significant in other cancer types
  table <- regression
  table <- markSigKS(table, sig_thres = sig_thres)
  sig_pairs <- vector('list', 3)
  for (cancer in c("BRCA", "OV", "CO")) {
    table_can <- table[table$Cancer == cancer & table$fdr_sig,]
    sig_pairs[[cancer]] <- as.vector(table_can$pair)
  }
  for (cancer in unique(table$Cancer)) {
    table[, paste0('sig_', cancer)] <- (table$pair %in% sig_pairs[[cancer]])
  }
  table$shared3can <- (table$sig_BRCA & table$sig_OV & table$sig_CO)
  table$uniq_BRCA <- (table$sig_BRCA & !table$sig_OV & !table$sig_CO)
  table$uniq_OV <- (!table$sig_BRCA & table$sig_OV & !table$sig_CO)
  table$uniq_CO <- (!table$sig_BRCA & !table$sig_OV & table$sig_CO)
  return(table)
}

remove_outliers <- function(x, out_thres, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- out_thres * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

formatPhosphosite = function(phosphosite_vector, gene_vector) {
  data <- phosphosite_vector
  gene <- gene_vector
  pho_rsd_split <- as.data.frame(str_split_fixed(data, ":", 2))
  colnames(pho_rsd_split) <- c("transcript", "SUB_MOD_RSD")
  pho_rsd_split$SUB_MOD_RSD <- toupper(pho_rsd_split$SUB_MOD_RSD)
  pho_rsd_split$SUBSTRATE <- gene
  return(pho_rsd_split)
}

load_ks_table <- function(protein_type) {
  ## Usage: k_s_table <- load_ks_table({kinase/phosphotase})
  if ( protein == "kinase" ) {
    ### read in the kinase/substrate table/ phosphorylation data ###
    k_s_table = read.delim(paste(baseD,"pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt",sep=""))
  }
  
  if ( protein == "phosphotase" ) {
    ### read in the phosphotase/substrate table/ phosphorylation data ### 
    k_s_table <- read.csv(paste(baseD,"pan3can_shared_data/Phospho_databases/DEPOD/DEPOD_201612_human_phosphatase-protein_substrate_to_Kuan-lin.csv",sep = ""))
    colnames(k_s_table) <- c("Phosphatase_UniProtAC_human","GENE","Substrate_UniProtAC_ref","SUB_GENE","Substrate_Type","DephosphoSite","BioassayType", "PubMed_ID_rev")
  }
  return(k_s_table)
}

