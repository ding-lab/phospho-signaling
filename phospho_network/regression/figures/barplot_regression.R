# Yige Wu @ WashU 2018 Jan
# barplot the validation statistics for regression results

# choose kinase/phosphotase, cancer , significance level -----------------------------------------------
protein <- "kinase"
top <- 10

# input -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')

resultDnow <- makeOutDir()
source(paste0(baseD, "cptac2p_analysis/pan3can_aes.R")) # aes for general purposes; it should be one directory out of the working directory
library(readr)
RTKs <- loadGeneList("RTK", cancer = "BRCA")

# loop around 3 cancer types individually -------------------------------------------------------
for (cancer in c('BRCA', 'OV', 'CO')) {
  for (self in c('cis', 'trans')) {
    tn = paste(resultD,"regression/tables/cptac2p_",cancer,"_",protein,"_", self, "_regression_validation_statistics.txt", sep="")
    valid_stat <- read_delim(tn, "\t", escape_double = FALSE, trim_ws = TRUE)
    valid_stat$Cancer <- cancer
    valid_stat$none <- valid_stat$all_count-valid_stat$pos_sig
    
    # barplot ranking validation count for trans----------------------------------------
    top <- 10
    valid_count_thres <- 1
    valid_stat0 <- valid_stat[!is.na(valid_stat$valid_ratio) & valid_stat$valid_count > 1 ,c("kinase","valid_count","Cancer","pos_sig","none")]
    valid_stat0 <- valid_stat0[order(valid_stat0$valid_count, decreasing = T),]
    valid_stat0$kinase_print <- valid_stat0$kinase
    driver_tag <- valid_stat0$kinase %in% driver_table$Gene[driver_table$`Cancer type` == driver_table_cancers[cancer]]
    valid_stat0$kinase_print[driver_tag] <- paste0('*', valid_stat0$kinase[driver_tag])
    valid_stat0 <- data.frame(valid_stat0)
    valid_stat_top <- data.frame(valid_stat0[1:(min(top,nrow(valid_stat0))),])
    
    ## plot individual top regulated kinases
    table_stat <- melt(valid_stat_top,id=c("kinase","valid_count","Cancer","kinase_print"))
    colnames(table_stat) <- c("kinase","valid_count","Cancer","kinase_print","coef_FDR","count")
    table_stat$KINASE <- reorder(table_stat$kinase_print, -table_stat$valid_count)
    table_stat$coef_FDR = as.character(table_stat$coef_FDR)
    table_stat$SELF <- self
    
    p <- ggplot()
    p <- p + geom_bar(data=table_stat, aes(y = count, x = KINASE, fill = coef_FDR ), stat="identity",
                      position='stack')
    p <- p + facet_grid(SELF~cancer, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p <- p + scale_fill_manual(values = c("pos_sig" = "red","none" = "grey"))
    p <- p + theme_bw() + theme_nogrid() + scale_y_log10()
    p <- p + xlab(protein)+ylab("number of substrate phosphosites")
    p <- p + theme(axis.title=element_text(size=10))
    p <- p + theme(axis.text.x = element_text(colour="black", size=10,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
    p
    fn = paste0(resultD,'regression/figures/cptac2p_', cancer, '_', self, '_',protein,'_substrate_validation_top_count_ranking_per_gene.pdf')
    ggsave(file=fn, height=3, width=4)
    
    ## plot regulated kinases landscape
    table_stat <- melt(valid_stat0,id=c("kinase","valid_count","Cancer","kinase_print"))
    colnames(table_stat) <- c("kinase","valid_count","Cancer","kinase_print","coef_FDR","count")
    table_stat$KINASE <- reorder(table_stat$kinase_print, -table_stat$valid_count)
    table_stat$coef_FDR = as.character(table_stat$coef_FDR)
    table_stat$SELF <- self
    
    p <- ggplot()
    p <- p + geom_bar(data=table_stat, aes(y = count, x = KINASE, fill = coef_FDR ), stat="identity",
                      position='stack')
    p <- p + facet_grid(SELF~Cancer, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p <- p + scale_fill_manual(values = c("pos_sig" = "red","none" = "grey"))
    p <- p + theme_bw() + theme_nogrid() + scale_y_log10()
    p <- p + xlab(protein)+ylab("number of substrate phosphosites")
    p <- p + theme(axis.title=element_text(size=10))
    p <- p + theme(axis.text.x = element_text(colour="black", size=10,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
    p
    fn = paste0(resultD,'regression/figures/cptac2p_', cancer, '_', self, '_',protein,'_substrate_validation_count_per_gene.pdf')
    ggsave(file=fn, height=6, width=10)
  }
}

# plot 3 cancers together -------------------------------------------------
for (self in c('cis', 'trans')) {
  valid_stat_top_3can <- NULL
  for (cancer in c('BRCA', 'OV', 'CO')) {
    tn = paste(resultD,"regression/tables/calculate_regression_validation/cptac2p_",cancer,"_",protein,"_", self, "_regression_validation_statistics.txt", sep="")
    valid_stat <- read_delim(tn, "\t", escape_double = FALSE, trim_ws = TRUE)
    valid_stat$Cancer <- cancer
    valid_stat$none <- valid_stat$all_count-valid_stat$valid_count
    valid_count_thres <- 1
    valid_stat0 <- valid_stat[!is.na(valid_stat$valid_ratio) & valid_stat$valid_count > 1,]
    valid_stat0 <- valid_stat0[order(valid_stat0$valid_count, decreasing = T),]
    valid_stat0$kinase_print <- valid_stat0$kinase
    drivers <- loadGeneList("driver", cancer)
    driver_tag <- (valid_stat0$kinase %in% drivers)
    rtk_tag <- (valid_stat0$kinase %in% RTKs)
    valid_stat0$kinase_print[driver_tag] <- paste0('*', valid_stat0$kinase_print[driver_tag])
    valid_stat0$kinase_print[rtk_tag] <- paste0('*', valid_stat0$kinase_print[rtk_tag])
    valid_stat0 <- data.frame(valid_stat0)
    valid_stat_top <- data.frame(valid_stat0[1:(min(top,nrow(valid_stat0))),])
    valid_stat_top_3can <- rbind(valid_stat_top_3can, valid_stat_top)
  }

  ## plot individual top regulated kinases
  valid_stat_top_3can$valid_not_retro <- valid_stat_top_3can$valid_count - valid_stat_top_3can$vad_retro_sitewise_count
  table_stat <- melt(valid_stat_top_3can,id=c("kinase", "valid_count", "all_count", "vad_retro_pairwise_count", "valid_ratio", "Cancer", "kinase_print"))
  table_stat$Phosphosite <- "Other"
  table_stat$Phosphosite[table_stat$variable == "valid_not_retro"] <- "regulated"
  table_stat$Phosphosite[table_stat$variable == "vad_retro_sitewise_count"] <- "regulated_also_in_retrospective"
  table_stat$KINASE <- reorder(table_stat$kinase_print, -table_stat$valid_count)
  table_stat$Cancer_f <- factor(table_stat$Cancer, levels = c("BRCA", "OV", "CO"))
  p <- ggplot()
  p <- p + geom_bar(data=table_stat, aes(y = value, x = KINASE, fill = Phosphosite ), stat="identity",
                    position='stack')
  p <- p + facet_grid(.~Cancer_f, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p <- p + scale_fill_manual(values = c("regulated_also_in_retrospective" = "#FF7F00", "regulated" = "#FF7F00","Other" = "grey"))
  p <- p + theme_bw()
  # p <- p + theme_nogrid()
  p <- p + scale_y_log10()
  p <- p + xlab(protein)+ylab("number of substrate phosphosites")
  p <- p + theme(axis.title=element_text(size=10))
  p <- p + theme(axis.text.x = element_text(colour="black", size=10,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
  p
  fn = paste0(resultDnow,'cptac2p_3can_', self, '_',protein,'_substrate_validation_top_count_ranking_per_gene.pdf')
  ggsave(file=fn, height=4, width=10)
  
  for (cancer in c("BRCA", "OV", "CO")) {
    table_can <- table_stat[table_stat$Cancer == cancer,]
    tmp <- as.vector(table_can$kinase_print)[order(-table_can$valid_count)]
    tmp <- tmp[!duplicated(tmp)]
    table_can$KINASE <- factor(table_can$kinase_print, levels = tmp)
    p <- ggplot()
    p <- p + geom_bar(data=table_can, aes(y = value, x = KINASE, fill = Phosphosite ), stat="identity",
                      position='stack')
    p <- p + facet_grid(.~Cancer_f, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p <- p + scale_fill_manual(values = c("regulated_also_in_retrospective" = "#FF7F00", "regulated" = "#FF7F00","Other" = "grey"))
    p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    p <- p + scale_y_log10()
    p <- p + xlab(protein)+ylab("number of substrate phosphosites")
    p <- p + theme(axis.title=element_text(size=10)) + theme(legend.position="none")
    p <- p + theme(axis.text.x = element_text(colour="black", size=10,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
    p
    fn = paste0(resultDnow,'cptac2p_', cancer, '_', self, '_',protein,'_substrate_validation_top_count_ranking_per_gene.pdf')
    ggsave(file=fn, height=4, width=4)
  }
}


