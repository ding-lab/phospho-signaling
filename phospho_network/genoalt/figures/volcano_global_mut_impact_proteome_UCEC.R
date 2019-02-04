# source ------------------------------------------------------------------
setwd("~/Box Sync/")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R")
library(ggrepel)

# set variables -----------------------------------------------------------
show_p_thres <- 0.2
sig_p_thres <- 0.05
cancer <- "UCEC"
num_genoalt_thres <- 5
smg_cancer <- c("FLNA",
                "MAP3K4",
                "HUWE1",
                "NSD1",
                "JAK1",
                "RPL22",
                "SCAF4",
                "FBXW7",
                "KMT2D",
                "INPPL1",
                "ZFHX3",
                "KMT2B",
                "TP53",
                "CTCF",
                "CTNNB1",
                "KRAS",
                "PIK3R1",
                "ARID1A",
                "PIK3CA",
                "PTEN")
text_cutoff <- -log10(0.05)
num_top2show <- 10

for (variant_class in c("not_silent")) {
  # plot protein section ----------------------------------------------------
  for (affected_exp_type in c("PRO", "PHO")) {
    for (genoalt_type in c("mut")) {
      mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_impact_proteome/", cancer, "_", variant_class, "_mut_impact_", affected_exp_type , "_tab.txt"), data.table = F, sep = "\t")
      mut_cnv_cans <- mut_cnv_cans[mut_cnv_cans$num >= num_genoalt_thres,]
      mut_cnv_cans$pair_pro <- paste(mut_cnv_cans$GENE, ":", mut_cnv_cans$SUB_GENE)
      mut_cnv_cans <- mut_cnv_cans[order(mut_cnv_cans$p),]
      
      for (SELF in c("cis", "trans")) {
        df <- mut_cnv_cans
        df <- df[df$SELF == SELF,]
        df$y <- -log10(df$p)
        df$x <- df$meddiff
        df_text <- df[df$p > 0 & df$p < sig_p_thres & df$GENE %in% smg_cancer,]
        df_text <- df_text[!duplicated(df_text$pair_pro),]
        df_text <- df_text[1:num_top2show,]
        cap <- min(max(abs(df$meddiff), na.rm = T), 1.5)
        df$meddiff_capped <- df$meddiff
        df$meddiff_capped[df$meddiff > cap] <- cap
        df$meddiff_capped[df$meddiff < (-cap)] <- (-cap)
        
        
        p = ggplot()
        p = p + geom_point(data = df, mapping = aes(x= x, y= y, color = meddiff_capped), 
                            alpha = 0.8, shape = 16, size = 2, stroke = 0)
        p = p + scale_color_gradientn(name= "protein\nchange", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
        p <- p + geom_text_repel(data = df_text, mapping = aes(x=x, y = y, label = pair), 
                                 size = 3, force = 2, color = "black")
        # p = p + scale_color_manual(values = c("up" = "#E31A1C", "down" = "#1F78B4"))
        p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
        p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
        p = p + labs(x = paste0("Median (Y-mutated samples - Y-unmutated samples) at protein X level"), y="-log10(P-value)")
        p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))

        p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
        p
        fn = paste(makeOutDir(resultD = resultD), variant_class, "_", genoalt_type, "_", SELF, "_impact_", affected_exp_type, "_num_genoalt_thres_", num_genoalt_thres, "_volcano.pdf",sep = "")
        ggsave(filename = fn, width = 6, height = 4, useDingbats = F)
        
      }
    }
  }
}