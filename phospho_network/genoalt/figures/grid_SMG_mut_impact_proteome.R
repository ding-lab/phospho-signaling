# Yige Wu @ WashU 2018 July
## overview figure for mutational impact of kinases on substrates

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R")


# set variables -----------------------------------------------------------
show_p_thres <- 0.2
fdr_thres <- 0.05
num_genoalt_thres <- 5
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")
affected_exp_type <- "PRO"
do.smg <- TRUE
pathway_order <- c("TP53", "Cell Cycle","WNT", "PI3K","RTK RAS", "MAPK", "SWI/SNF", "NOTCH", "MMR", "HIPPO", "TGF-Beta", "other")

# input pair annotation table ---------------------------------------------
pair_tab_annotated <- fread(input = "./cptac2p/analysis_results/dependencies/tables/compile_protein_pair_table/protein_pair_table.txt", data.table = F)


# input mutation impact table ---------------------------------------------
mut_impact_tab <- fread(input = "./cptac2p/analysis_results/phospho_network/genoalt/tables/test_SMG_mut_impact_proteome/SMG_mut_impact_tab.txt", data.table = F)
mut_impact_tab$GENE.is_SMG <- get_SMG_by_cancer(gene_vector = mut_impact_tab$GENE, cancer_vector = mut_impact_tab$cancer)
mut_impact_tab <- mut_impact_tab %>%
  mutate(GENE_cancer = paste0(GENE, ":", cancer)) %>%
  mutate(SUB_GENE_cancer = paste0(SUB_GENE, ":", cancer))

# get the SMG-kinase/phosphatase to show ----------------------------------
df <- mut_impact_tab
df <- df %>%
  filter(GENE.is_SMG == T) %>%
  filter(affected_exp_type == "PHO") %>%
  filter(cancer %in% cancers2process) %>%
  filter(num >= num_genoalt_thres)

gene_altered2show_kinase <- df %>%
  filter(fdr_by_gene < fdr_thres) %>%
  filter(GENE %in% c(kinases, phosphatases) & (SUB_GENE.is_kinase_substrate == T | SUB_GENE.is_phosphatase_substrate == T)) %>%
  select(GENE_cancer) %>%
  unique()
gene_altered2show_kinase <- gene_altered2show_kinase$GENE_cancer

gene_affected2show_kinase <- df %>%
  filter(fdr_by_gene < fdr_thres) %>%
  filter(GENE_cancer %in% gene_altered2show_kinase & (SUB_GENE.is_kinase_substrate == T | SUB_GENE.is_phosphatase_substrate == T)) %>%
  select(SUB_GENE_cancer) %>%
  unique()
gene_affected2show_kinase <- gene_affected2show_kinase$SUB_GENE_cancer

# get the SMG-TF to show ----------------------------------
df <- mut_impact_tab
df <- df %>%
  filter(GENE.is_SMG == T) %>%
  filter(affected_exp_type == "PRO") %>%
  filter(cancer %in% cancers2process) %>%
  filter(num >= num_genoalt_thres)

gene_altered2show_TF <- df %>%
  filter(fdr_by_gene < fdr_thres) %>%
  filter(SUB_GENE.is_TF_downstream == T) %>%
  filter(SUB_GENE %in% c(kinases, phosphatases)) %>%
  select(GENE_cancer) %>%
  unique()
gene_altered2show_TF
gene_altered2show_TF <- gene_altered2show_TF$GENE_cancer

gene_affected2show_TF <- df %>%
  filter(fdr_by_gene < fdr_thres) %>%
  filter(GENE_cancer %in% gene_altered2show_TF & (SUB_GENE.is_TF_downstream == T)) %>%
  filter(SUB_GENE %in% c(kinases, phosphatases)) %>%
  select(SUB_GENE_cancer) %>%
  unique()
gene_affected2show_TF <- gene_affected2show_TF$SUB_GENE_cancer

# get the SMG-other to show ----------------------------------
gene_altered2show_other <- mut_impact_tab %>% 
  filter(GENE.is_SMG == T) %>%
  filter(cancer %in% cancers2process) %>%
  filter(!(GENE %in% c(kinases, phosphatases))) %>%
  filter(!(GENE %in% pair_tab_annotated$GENE[pair_tab_annotated$SUB_GENE.is_TF_downstream])) %>%
  filter((SUB_GENE %in% c(kinases, phosphatases))) %>%
  filter(fdr_by_gene < fdr_thres) %>%
  filter(num >= num_genoalt_thres) %>%
  filter(affected_exp_type == "PHO" | affected_exp_type == "PRO") %>%
  filter(SELF == "trans") %>%
  select(GENE_cancer) %>%
  unique()
gene_altered2show_other <- gene_altered2show_other$GENE_cancer

gene_affected2show_other <- mut_impact_tab %>% 
  filter(GENE.is_SMG == T) %>%
  filter(cancer %in% cancers2process) %>%
  filter(fdr_by_gene < fdr_thres) %>%
  filter(GENE_cancer %in% gene_altered2show_other) %>%
  filter(SUB_GENE %in% c(kinases, phosphatases)) %>%
  filter(affected_exp_type == "PHO" | affected_exp_type == "PRO") %>%
  select(SUB_GENE_cancer) %>%
  unique()
gene_affected2show_other <- gene_affected2show_other$SUB_GENE_cancer

# plot by RNA/PRO/PHO -----------------------------------------------------
for (affected_exp_type_tmp in c("RNA", "PRO", "PHO")) {
  df <- mut_impact_tab
  df <- df %>%
    filter(affected_exp_type == affected_exp_type_tmp) %>%
    filter(cancer %in% cancers2process) %>%
    filter(num >= num_genoalt_thres) %>%
    filter(GENE_cancer %in% c(gene_altered2show_kinase, gene_altered2show_TF, gene_altered2show_other)) %>%
    filter(SUB_GENE_cancer %in% c(gene_affected2show_kinase, gene_affected2show_TF, gene_affected2show_other) | SUB_GENE_cancer %in% c(gene_altered2show_kinase, gene_altered2show_TF, gene_altered2show_other))
    
  df$sig <- (df$fdr_by_gene < fdr_thres)
  df$pair_can  <- paste0(df$pair, ":", df$cancer)
  df <- df[order(df$p, decreasing = T),]
  # df <- df[!duplicated(df$pair_can),]
  df$log10_FDR <- -log10(df$fdr_by_gene)
  
  cap <- min(max(abs(df$meddiff), na.rm = T), 2)
  df$meddiff_capped <- df$meddiff
  df$meddiff_capped[df$meddiff > cap] <- cap
  df$meddiff_capped[df$meddiff < (-cap)] <- (-cap)
  
  sig_colors <- c("black", "w hite")
  names(sig_colors) <- c(paste0('FDR<', fdr_thres), paste0('FDR>', fdr_thres))
  
  # tmp <- as.vector(unique(df$GENE))
  # df$GENE <- factor(df$GENE, levels = c(c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"), tmp[!(tmp %in% c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"))]))
  
  df$shape <- ifelse(df$SELF == "cis", "a", "b")
  df$shape <- as.factor(df$shape)
  pos_up <- position_nudge(x = 0, y = 0.3)
  pos_down <- position_nudge(x = 0, y = -0.3)
  if (affected_exp_type == "PRO" | affected_exp_type == "RNA") {
    p <- ggplot()
    p <- p + geom_point(data = df[df$variant_class == "not_silent",], mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_FDR, shape = variant_class,
                                                 colour = ifelse(sig == TRUE, paste0('FDR<', fdr_thres), paste0('FDR>', fdr_thres))), alpha = 0.5, stroke = 1)
    p <- p + geom_point(data = df[df$variant_class == "missense",], mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_FDR, shape = variant_class,
                                                                                   colour = ifelse(sig == TRUE, paste0('FDR<', fdr_thres), paste0('FDR>', fdr_thres))), alpha = 0.5, stroke = 1, position = pos_up)
    p <- p + geom_point(data = df[df$variant_class == "truncation",], mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_FDR, shape = variant_class,
                                                                                   colour = ifelse(sig == TRUE, paste0('FDR<', fdr_thres), paste0('FDR>', fdr_thres))), alpha = 0.5, stroke = 1, position = pos_down)
    p <- p  + scale_shape_manual(values = c("not_silent" = 21, "missense" = 22, "truncation" = 23))
    p <- p + scale_color_manual(values = sig_colors)
    p = p + scale_fill_gradientn(name= "protein\nchange", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
    p <- p + facet_grid(cancer ~ ., drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p <- p + xlab("Affected Protein") + ylab("Mutated Gene")
    p <- p + guides(colour = guide_legend(title = "FDR"))
    p <- p + theme_bw()
    p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
    p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
    p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
    p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
    p <- p + theme(strip.text.y = element_text(angle = 0))
    p <- p + theme(strip.background = element_rect(fill = "white", color = "white"))
    p
    fn = paste(makeOutDir(resultD = resultD), "mut_impact_", affected_exp_type_tmp, "_num_genoalt_thres_", num_genoalt_thres,'.pdf',sep = "")
    ggsave(filename = fn, width = 8, height = 6)
  }
  
  if (affected_exp_type == "PHO") {
    p <- ggplot()
    p <- p + geom_point(data = df[df$variant_class == "not_silent",], mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_FDR, shape = variant_class,
                                                                                    colour = ifelse(sig == TRUE, paste0('FDR<', fdr_thres), paste0('FDR>', fdr_thres))), alpha = 0.5, stroke = 1)
    p <- p + geom_point(data = df[df$variant_class == "missense",], mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_FDR, shape = variant_class,
                                                                                  colour = ifelse(sig == TRUE, paste0('FDR<', fdr_thres), paste0('FDR>', fdr_thres))), alpha = 0.5, stroke = 1, position = pos_up)
    p <- p + geom_point(data = df[df$variant_class == "truncation",], mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_FDR, shape = variant_class,
                                                                                    colour = ifelse(sig == TRUE, paste0('FDR<', fdr_thres), paste0('FDR>', fdr_thres))), alpha = 0.5, stroke = 1, position = pos_down)
    p <- p  + scale_shape_manual(values = c("not_silent" = 21, "missense" = 22, "truncation" = 23))
    p <- p + scale_color_manual(values = sig_colors)
    p = p + scale_fill_gradientn(name= "phosphoprotein\nchange", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
    p <- p + facet_grid(cancer ~ ., drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
    p <- p + xlab("Affected Phosphorotein") + ylab("Mutated Gene")
    p <- p + guides(colour = guide_legend(title = "FDR"))
    p <- p + theme_bw()
    p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
    p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
    p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
    p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
    p <- p + theme(strip.text.y = element_text(angle = 0))
    p <- p + theme(strip.background = element_rect(fill = "white", color = "white"))
    p
    fn = paste(makeOutDir(resultD = resultD), "mut_impact_", affected_exp_type_tmp, "_num_genoalt_thres_", num_genoalt_thres,'.pdf',sep = "")
    ggsave(filename = fn, width = 8, height = 6)
  }
}

# Plot RNA/PRO/PHO tgt ----------------------------------------------------
df <- mut_impact_tab
df <- df %>%
  filter(cancer %in% cancers2process) %>%
  filter(num >= num_genoalt_thres) %>%
  filter(GENE_cancer %in% c(gene_altered2show_kinase, gene_altered2show_TF, gene_altered2show_other)) %>%
  filter(SUB_GENE_cancer %in% c(gene_affected2show_kinase, gene_affected2show_TF, gene_affected2show_other) | SUB_GENE_cancer %in% c(gene_altered2show_kinase, gene_altered2show_TF, gene_altered2show_other))

df$sig <- (df$fdr_by_gene < fdr_thres)
df$pair_can  <- paste0(df$pair, ":", df$cancer)
df <- df[order(df$p, decreasing = T),]
# df <- df[!duplicated(df$pair_can),]
df$log10_FDR <- -log10(df$fdr_by_gene)

cap <- min(max(abs(df$meddiff), na.rm = T), 2)
df$meddiff_capped <- df$meddiff
df$meddiff_capped[df$meddiff > cap] <- cap
df$meddiff_capped[df$meddiff < (-cap)] <- (-cap)

sig_colors <- c("black", "w hite")
names(sig_colors) <- c(paste0('FDR<', fdr_thres), paste0('FDR>', fdr_thres))

# tmp <- as.vector(unique(df$GENE))
# df$GENE <- factor(df$GENE, levels = c(c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"), tmp[!(tmp %in% c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"))]))

df$shape <- ifelse(df$SELF == "cis", "a", "b")
df$shape <- as.factor(df$shape)
pos_up <- position_nudge(x = 0, y = 0.3)
pos_down <- position_nudge(x = 0, y = -0.3)

df$affected_exp_type <- factor(df$affected_exp_type, levels = c("RNA", "PRO", "PHO"))
p <- ggplot()
p <- p + geom_point(data = df[df$variant_class == "not_silent",], mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_FDR, shape = variant_class,
                                                                                colour = ifelse(sig == TRUE, paste0('FDR<', fdr_thres), paste0('FDR>', fdr_thres))), alpha = 0.5, stroke = 1)
p <- p + geom_point(data = df[df$variant_class == "missense",], mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_FDR, shape = variant_class,
                                                                              colour = ifelse(sig == TRUE, paste0('FDR<', fdr_thres), paste0('FDR>', fdr_thres))), alpha = 0.5, stroke = 1, position = pos_up)
p <- p + geom_point(data = df[df$variant_class == "truncation",], mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_FDR, shape = variant_class,
                                                                                colour = ifelse(sig == TRUE, paste0('FDR<', fdr_thres), paste0('FDR>', fdr_thres))), alpha = 0.5, stroke = 1, position = pos_down)
p <- p  + scale_shape_manual(values = c("not_silent" = 21, "missense" = 22, "truncation" = 23))
p <- p + scale_color_manual(values = sig_colors)
p = p + scale_fill_gradientn(name= "protein\nchange", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
p <- p + facet_grid(cancer ~ affected_exp_type, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
p <- p + xlab("Affected Protein") + ylab("Mutated Gene")
p <- p + guides(colour = guide_legend(title = "FDR"))
p <- p + theme_bw()
p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
p <- p + theme(strip.text.y = element_text(angle = 0))
p <- p + theme(strip.background = element_rect(fill = "white", color = "white"))
p
fn = paste(makeOutDir(resultD = resultD), "mut_impact_", "RNA_PRO_PHO", "_num_genoalt_thres_", num_genoalt_thres,'.pdf',sep = "")
ggsave(filename = fn, width = 1, height = 6)
