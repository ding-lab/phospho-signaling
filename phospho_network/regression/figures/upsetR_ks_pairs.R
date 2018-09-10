# Yige Wu @ WashU 2018 Jan
# look at the intersections of regulated kinase-substrate pairs in BRCA, OV and CO datasets


# input -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
resultDnow <- makeOutDir()
source(paste0(baseD, "cptac2p_analysis/pan3can_aes.R")) # aes for general purposes
library(UpSetR)
plotshared <- function(mydata, x, y){
  shared_or_uniq = y
  top <- 10
  mydata2 <- mydata[mydata[,shared_or_uniq],]
  if (nrow(mydata2) > 0) {
    mydata3 <- data.frame(table(mydata2$KINASE))
    colnames(mydata3)[1] <- 'KINASE'
    mydata3 <- mydata3[order(mydata3$Freq, decreasing = T),]
    mydata4 <- mydata3[1:(min(nrow(mydata3), top)),]
    mydata4[,x] <- reorder(mydata4$KINASE, -mydata4$Freq)
    myplot <- (ggplot(mydata4, aes_string(x= x, y = 'Freq')) + 
                 geom_bar(stat="identity",position='stack') +
                 theme_bw() + theme_nogrid() + scale_y_log10() +
                 ggtitle(label = paste0("Kinases of ks pairs(", shared_or_uniq, ")")) +
                 theme(axis.text.x = element_text(colour="black", size = 5, angle= 30, vjust=0.5), 
                       axis.text.y = element_text(colour="black", size=3),
                       axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       plot.title = element_text(size=5)))
  } else {
    print('not enough data!')
  }
}

# find shared/unique ks pairs to specific phosphosites ------------------------------------------
for (protein in c("kinase", "phosphotase")) {
  tn = paste0(resultD,"regression/tables/",protein,"_substrate_regression_cptac2p_3can.txt")
  table_3can <- fread(tn, data.table = F)
  table_3can <- markSigSiteCan(table_3can, sig_thres = sig, enzyme_type = protein)
  ## only look at ks pairs that are examined in all 3 cancer types
  table4upset <- table_3can[table_3can$fdr_sig & !is.na(table_3can$shared3can),]
  table4upset <- cbind(table4upset[,c('KINASE', 'SUBSTRATE', 'SUB_MOD_RSD', 'SELF', 'shared3can', 'uniq_BRCA', 'uniq_OV','uniq_CO')], 
                       table4upset[,c('sig_BRCA', 'sig_OV', 'sig_CO')]*1)
  table4upset <- unique(table4upset)
  for (self in c('cis', 'trans')) {
    table4upset_tmp <- table4upset[table4upset$SELF == self,]
    if (nrow(table4upset_tmp) > 0) {
      ## upSet subplot
      subplots <- list()
      for (geneset in c('shared3can', 'uniq_BRCA', 'uniq_OV', 'uniq_CO')) {
        if (length(which(table4upset_tmp[, geneset])) > 0) {
          subplots[[length(subplots) + 1]] <- list(plot=plotshared, x=protein, y = geneset)
        }
      }
      
      ## main upset plot
      fn = paste(resultDnow, "cptac2p_3can_", protein,'_regrssion_',self, "_FDR_", sig, '_SitewiseSharedorUnique.pdf',sep ="")
      pdf(fn, height = 8, width = 6, onefile = F)
      upset(table4upset_tmp, sets = c('sig_BRCA', 'sig_OV', 'sig_CO'), sets.bar.color = "#56B4E9",
            order.by = "freq", empty.intersections = "on",
            attribute.plots=list(gridrows=60,plots=subplots, nrow = 2))
      dev.off()
    }
  }
}

# find shared/unique ks pairs to just kinase and substrate protein ------------------------------------------
for (protein in c("kinase", "phosphotase")) {
  tn = paste0(resultD,"regression/tables/",protein,"_substrate_regression_cptac2p_3can.txt")
  table_3can <- read.delim(file = tn)
  table_3can <- markSigCan(table_3can, sig_thres = sig)
  table4upset <- table_3can[table_3can$fdr_sig,]
  table4upset <- cbind(table4upset[,c('KINASE', 'SUBSTRATE', 'SUB_MOD_RSD', 'SELF', 'shared3can', 'uniq_BRCA', 'uniq_OV','uniq_CO')], 
                       table4upset[,c('ks_sig_BRCA', 'ks_sig_OV', 'ks_sig_CO')]*1)
  table4upset <- unique(table4upset)
  for (self in c('cis', 'trans')) {
    table4upset_tmp <- table4upset[table4upset$SELF == self,]
    ## set up main upset plot
    fn = paste(resultDnow, "cptac2p_3can_", protein,'_regrssion_',self, "_FDR_", sig, '_KSpairSharedorUnique.pdf',sep ="")
    pdf(fn, height = 8, width = 6)
    upset(table4upset_tmp, sets = c('ks_sig_BRCA', 'ks_sig_OV', 'ks_sig_CO'), sets.bar.color = "#56B4E9",
          order.by = "freq", empty.intersections = "on",
          attribute.plots=list(gridrows=60,plots=list(list(plot=plotshared, x=protein, y = 'shared3can'),
                                                      list(plot=plotshared, x=protein, y = 'uniq_BRCA'),
                                                      list(plot=plotshared, x=protein, y = 'uniq_OV'),
                                                      list(plot=plotshared, x=protein, y = 'uniq_CO')), nrow = 2))
    dev.off()
  }
}
