# Yige Wu @ WashU 2019 Mar
## plot volcano plot for known partners affected by TP53 mutations


# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/p53/TP53_shared.R")


# input result table ------------------------------------------------------
results<- fread("./cptac2p/analysis_results/p53/tables/test_mut_impact_proteome_TP53/TP53_mut_impact_proteome_RNA_cptac2p_cptac3_tab.txt", data.table = F)


# set variables -----------------------------------------------------------
num_genoalt_thres <- 5
text_cutoff <- -log10(0.05)
num_top2show <- 20

for (mut_type in c('missense', 'truncation', 'not_silent')) {
  print(mut_type)
  for (type in c( 'PRO', 'RNA', 'PHO')) {
    mut_cnv_cans <- results
    mut_cnv_cans <- mut_cnv_cans[mut_cnv_cans$num >= num_genoalt_thres,]
    mut_cnv_cans <- mut_cnv_cans[order(mut_cnv_cans$p),]
    mut_cnv_cans <- subset(mut_cnv_cans, variant_class==mut_type)
    RdBu1024 <- colorRampPalette(rev(brewer.pal(10, "RdBu")))(1024)
    mut_cnv_cans <- subset(mut_cnv_cans, affected_exp_type==type )
    
    for (self in c("trans")) {
      print(self)
      tab2p <- mut_cnv_cans
      tab2p <- subset(tab2p, SELF == self)
      tab2p$y <- -log10(tab2p$fdr)
      tab2p$x <- tab2p$meddiff
      tab2p_text <- tab2p[tab2p$fdr > 0 & tab2p$fdr < sig_p_thres,]
      #tab2p_text <- tab2p_text[!duplicated(tab2p_text$pair_pro),]
      tab2p_text <- tab2p_text[1:num_top2show,]
      if (type == "PHO") {
        tab2p_text$text <- str_split_fixed(tab2p_text$pair, ':', 2)[,2]
      } else {
        tab2p_text$text <- str_split_fixed(tab2p_text$pair, ':', 3)[,2]
      }
      
      cap <- min(max(abs(tab2p$meddiff), na.rm = T), 1)
      tab2p$meddiff_capped <- tab2p$meddiff
      tab2p$meddiff_capped[tab2p$meddiff > cap] <- cap
      tab2p$meddiff_capped[tab2p$meddiff < (-cap)] <- (-cap)
      
      
      p = ggplot()
      p = p + geom_point(data = tab2p, mapping = aes(x= x, y= y, fill = meddiff_capped), 
                         alpha = 1, shape = 21, size = 2, stroke = 0, color = "white")
      # p = p + geom_point(data = subset(tab2p, p>0.05), mapping = aes(x= x, y= y), 
      #                    alpha = 0.5, shape = 16, size = 2, stroke = 0, color = 'lightgrey')
      p = p + scale_fill_gradientn(name= paste0(type, " change"), na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
      #p = p + scale_color_brewer(palette = 'Dark2')
      p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
      p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
      p <- p + geom_text_repel(data = tab2p_text, mapping = aes(x=x, y = y, label = text, color = is_uniq.variant_class),
                               size = 3, force = 2)
      p = p + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), name = paste0("only associated in ", mut_type))
      p = p + facet_grid(~cancer, drop=T, space = "free",scales = "free")
      p = p + labs(x = paste0("log2(Median Fold Change ", mut_type, " vs WT)"), y="-log10(FDR P-value)")
      p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
      p <- p + theme_classic() #+ theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
      p
      
      fn <- paste0(makeOutDir(resultD = resultD), '/Volcano_plot_', self ,'_tp53_', mut_type,'_impact_on_', type,'_global_search_FDR.pdf')
      ggsave(filename = fn, width = 15, height = 4, useDingbats = F)
      
    }
  }
}