# source ------------------------------------------------------------------
setwd("~/Box Sync/")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R")
library(ggrepel)
library(dplyr)

# set variables -----------------------------------------------------------
show_p_thres <- 0.2
sig_p_thres <- 0.05
cancer <- "UCEC"
num_genoalt_thres <- 5
gene_alt <- "TP53"
text_cutoff <- -log10(0.05)
num_top2show <- 10
cancers2process <- c("BRCA", "UCEC")
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC", "LIHC")

# inputs ------------------------------------------------------------------
## input enzyme-substrate table
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)
### take out pairs that are purely computationally predicted
ptms_site_pairs_sup <- ptms_site_pairs_sup[!(ptms_site_pairs_sup$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP")) ,]
ptms_site_pairs_sup$pair <- paste0(ptms_site_pairs_sup$pair_pro, ":", ptms_site_pairs_sup$SUB_MOD_RSD)
ptms_site_pairs_sup$pair_pro_rev <- paste0(ptms_site_pairs_sup$SUB_GENE, ":",  ptms_site_pairs_sup$GENE)

mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_impact_proteome/mut_impact_proteome_RNA_cptac2p_cptac3_tab.txt"), data.table = F, sep = "\t")
mut_cnv_cans <- mut_cnv_cans[mut_cnv_cans$num >= num_genoalt_thres,]
mut_cnv_cans %>%
  select(cancer) %>%
  unique()
mut_cnv_cans$pair %>% head()
mut_cnv_cans <- merge(mut_cnv_cans, ptms_site_pairs_sup[, c("pair", "Source")], by = c("pair"), all.x = T)
mut_cnv_cans$pair_pro <- paste0(mut_cnv_cans$GENE, ":",  mut_cnv_cans$SUB_GENE)

# ``` {r annotate-mut-impact-table-with-kinase/phosphatase}
mut_cnv_cans$SUB_GENE_role <- "complex_partner"
mut_cnv_cans$SUB_GENE_role[mut_cnv_cans$pair_pro %in% ptms_site_pairs_sup$pair_pro] <- "substrate"
mut_cnv_cans$SUB_GENE_role[mut_cnv_cans$pair_pro %in% ptms_site_pairs_sup$pair_pro_rev[ptms_site_pairs_sup$enzyme_type == "kinase"]] <- "kinase"
mut_cnv_cans$SUB_GENE_role[mut_cnv_cans$pair_pro %in% ptms_site_pairs_sup$pair_pro_rev[ptms_site_pairs_sup$enzyme_type == "phosphatase"]] <- "phosphatase"
# ```

for (SUB_GENE_role in c("kinase_phosphatase")) {
# for (SUB_GENE_role in c("kinase", "phosphatase", "complex_partner")) {
  for (variant_class in c("missense", "truncation", "not_silent")) {
  # for (variant_class in c("not_silent")) {
      
    # plot protein section ----------------------------------------------------
    # for (affected_exp_type in c("PRO", "PHO", "RNA")) {
    for (affected_exp_type in c("PRO", "PHO")) {
        
      for (genoalt_type in c("mut")) {
        for (SELF in c("trans")) {
          df <- mut_cnv_cans
          if (affected_exp_type %in% c("PRO", "RNA")) {
            df <- df[ df$SUB_MOD_RSD == affected_exp_type,]
            df <- df[df$SELF == SELF & df$num >= num_genoalt_thres & (df$GENE ==  gene_alt) & (df$cancer %in% cancers2process) & df$variant_class == variant_class & df$SUB_MOD_RSD == affected_exp_type,]
          } else {
            df <- df[df$SELF == SELF & df$num >= num_genoalt_thres & df$GENE ==  gene_alt & df$cancer %in% cancers2process & df$variant_class == variant_class & !(df$SUB_MOD_RSD %in% c("PRO", "RNA")) ,]
          }
          
          if (SUB_GENE_role == "kinase_phosphatase") {
            df <- df[df$SUB_GENE_role %in% c("kinase", "phosphatase"),]
            
          } else {
            df <- df[df$SUB_GENE_role == SUB_GENE_role,]
            
          }
          df$y <- -log10(df$p)
          df$x <- df$meddiff
          df_text <- df[df$p > 0 & df$p < sig_p_thres,]
          df_text_order <- data.frame(table(df$pair_pro))
          df_text_order <- df_text_order[order(df_text_order$Freq, decreasing = T),]
          df_text_order_rev <- df_text_order[order(df_text_order$Freq, decreasing = F),]
          head(df_text_order)
          df_text <- df_text[!duplicated(df_text$pair_pro),]
          # df_text <- df_text[df_text$pair_pro %in% c(as.vector(df_text_order$Var1[1:num_top2show]), as.vector(df_text_order_rev$Var1[1:num_top2show])),]
          cap <- min(max(abs(df$meddiff), na.rm = T), 1.5)
          df$meddiff_capped <- df$meddiff
          df$meddiff_capped[df$meddiff > cap] <- cap
          df$meddiff_capped[df$meddiff < (-cap)] <- (-cap)
          df$cancer <- factor(as.vector(df$cancer), levels = rev(cancers2process))
          
          p = ggplot()
          p = p + geom_point(data = df, mapping = aes(x= x, y= y, color = meddiff_capped), 
                             alpha = 0.8, shape = 16, size = 2, stroke = 0)
          p <- p + geom_text_repel(data = df_text, mapping = aes(x=x, y = y, label = pair), 
                                   size = 2, force = 2, color = "black")
          p = p + scale_color_gradientn(name= paste0(affected_exp_type, "\nchange"), na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
          # p = p + scale_color_manual(values = c("up" = "#E31A1C", "down" = "#1F78B4"))
          p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
          p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
          # p = p + facet_grid(cancer~., scales = "free_y", space = "fixed")
          p = p + facet_wrap(c("cancer"), scales = "free_y")
          p = p + labs(x = paste0("Median (Y-mutated samples - Y-unmutated samples) at X " , affected_exp_type, " level"), y="-log10(P-value)")
          p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))
          
          p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
          p
          fn = paste(makeOutDir(resultD = resultD), gene_alt, "_", variant_class, "_", genoalt_type, "_", SELF, "_impact_its_", SUB_GENE_role, "_", affected_exp_type, "_num_genoalt_thres_", num_genoalt_thres, "_volcano.pdf",sep = "")
          if (length(unique(df$cancer)) == 2) {
            ggsave(filename = fn, width = 6, height = 4, useDingbats = F)
            
          } else {
            # ggsave(filename = fn, width = 5, height = 10, useDingbats = F)
            ggsave(filename = fn, width = 9, height = 5, useDingbats = F)
            
          }
        }
      }
    }
  }
}

