# Yige Wu @ WashU 2018 Jan
# look at the intersections of consistently high/low kinase-substrate pairs in BRCA, OV and CO datasets


# source -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source(paste0(baseD, "cptac2p_analysis/pan3can_aes.R")) # aes for general purposes
library(UpSetR)
plotshared <- function(mydata, x, y){
  shared_or_uniq = y
  top <- 10
  mydata2 <- mydata[mydata[,shared_or_uniq],]
  if (nrow(mydata2) > 0) {
    mydata3 <- data.frame(table(mydata2$GENE))
    colnames(mydata3)[1] <- 'KINASE'
    mydata3 <- mydata3[order(as.vector(mydata3[,"Freq"]), decreasing = T),]
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
resultDnow <- makeOutDir()
cancers <- c("BRCA", "OV", "CO")

# mark cancer specific or shared kinase-substrate pairs -------------------
## formatting
cons_tabs <- list()
pairs_uniq <- NULL
for (cancer in c("BRCA", "OV", "CO")) {
# for (cancer in c("BRCA")) {
  ksea_diffexp <- fread(input = paste0("cptac2p/cptac_shared/analysis_results/phospho_network/kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_regression/", 
                                       cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)

  cons_tab <- ksea_diffexp[!is.na(ksea_diffexp$consistent) ,]
  cons_tab$pair <- paste0(cons_tab$GENE, ":", cons_tab$SUB_GENE, ":", cons_tab$SUB_MOD_RSD)
  cons_tab <- cons_tab[!duplicated(cons_tab$pair),]
  cons_tab$id <- paste0(cons_tab$pair, ":", cons_tab$enzyme_direction, ":", cons_tab$substrate_direction)
  rownames(cons_tab) <- cons_tab$id
  cons_tabs[[cancer]] <- cons_tab
  pairs_uniq <- unique(rbind(pairs_uniq, cons_tab[cons_tab$consistent, c("GENE", "SUB_GENE", "SUB_MOD_RSD", "pair", "enzyme_type", "enzyme_direction" ,"substrate_direction", "id")]))
}

## mark consistently high/low pairs in each cancer types
table4upset <- data.frame(pairs_uniq)
table4upset <- table4upset[order(table4upset$GENE, table4upset$SUB_GENE, table4upset$SUB_MOD_RSD),]
for (cancer in c("BRCA", "OV", "CO")) {
  cons_tab <- cons_tabs[[cancer]]
  sig_pairs <- cons_tab$id[cons_tab$consistent]
  table4upset[, paste0("sig_", cancer)] <- (table4upset$id %in% sig_pairs)
}

## mark unique or shared across cancer types
for (cancer in cancers) {
  rest_cancers <- cancers[cancers != cancer]
  table4upset[, paste0("uniq_", cancer)] <- (table4upset[, paste0("sig_", cancer)] & rowSums(table4upset[, paste0("sig_", rest_cancers)]) == 0)
}
table4upset[, "shared"] <- (rowSums(table4upset[, paste0("sig_", cancers)]) == length(cancers))
table4upset[, paste0("sig_", cancers)] <- table4upset[, paste0("sig_", cancers)]*1
table4upset <- data.frame(table4upset)
# find shared/unique ks pairs to specific phosphosites ------------------------------------------
for (enzyme_type in c("kinase", "phosphotase")) {
# for (enzyme_type in c("kinase")) {
  for (sub_direction in c('up', 'down')) {
    table4upset_tmp <- table4upset[table4upset$enzyme_type == enzyme_type & table4upset$substrate_direction == sub_direction,]
    if (nrow(table4upset_tmp) > 0) {
      ## upSet subplot
      subplots <- list()
      for (geneset in c('shared', 'uniq_BRCA', 'uniq_OV', 'uniq_CO')) {
        if (length(which(table4upset_tmp[, geneset])) > 0) {
          subplots[[length(subplots) + 1]] <- list(plot=plotshared, x=enzyme_type, y = geneset)
        }
      }
      
      ## main upset plot
      fn = paste(resultDnow, "cptac2p_3can_", enzyme_type,'_substrate_',sub_direction, "_FDR_", sig, '_SitewiseSharedorUnique.pdf',sep ="")
      pdf(fn, height = 8, width = 6, onefile = F)
      upset(table4upset_tmp, sets = c('sig_BRCA', 'sig_OV', 'sig_CO'), sets.bar.color = "#56B4E9",
            order.by = "freq", empty.intersections = "on",
            attribute.plots=list(gridrows=60,plots=subplots, nrow = 2))
      dev.off()
    }
  }
}

