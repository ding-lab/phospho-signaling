# Yige Wu @ WashU 2019 Jan
## test mutation impact on protein/phosphorylation within kinase-substrate pairs or protein complex pairs


# source ------------------------------------------------------------------
setwd(dir = "~/Box Sync/")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/p53/TP53_shared.R")

if(!("pheatmap" %in% installed.packages()[,"Package"])) {
  install.packages("pheatmap")
}
library(pheatmap)
save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf <- function(x, filename, width=6, height=6) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
library(RColorBrewer)

# set variables -----------------------------------------------------------
gene <- "TP53"

# inputs ------------------------------------------------------------------
for (cancer in c("BRCA", "OV", "CO", "CCRCC", "UCEC", "LIHC")) {
# for (cancer in c("CO")) {
  sup_tab <- NULL
 
  maf <- loadMaf(cancer = cancer, maf_files = maf_files)
  
  if (cancer %in% c("CCRCC", "UCEC", "OV")) {
    maf_gene <- maf %>%
      filter(Hugo_Symbol == gene)
    
    maf_gene <- maf %>%
      filter(Hugo_Symbol == gene) %>%
      mutate(VAF = as.numeric(t_alt_count)/as.numeric(t_depth))
  }
  
  if (cancer %in% c("BRCA", "CO")) {
    cancer_tmp <- cancer
    source("./cptac2p_analysis/p53/tables/extract_VAF_per_gene.R")
    maf_gene <- maf_gene_new
  }
  
  if (cancer %in% c("LIHC")) {
    maf_gene <- maf %>%
      filter(Hugo_Symbol == gene) %>%
      mutate(VAF = as.numeric(n_alt_count)/as.numeric(n_depth))
  }
  
  vaf_part <- maf_gene %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(VAF = max(VAF, na.rm = T)) %>%
    mutate(partID = str_split_fixed(string = Tumor_Sample_Barcode, pattern = "_", n = 2)[,1])
  vaf_part$VAF[is.infinite(vaf_part$VAF)] <- NA
  
  mut_mat <- generate_somatic_mutation_matrix(pair_tab = gene, maf = maf)
  
  ## input CNA matrix
  # cna_tab <- loadCNAstatus(cancer = cancer)
  cna_tab <- loadTP53Deletion(cancer = cancer)
  
  cna_tab <- cna_tab[cna_tab$gene == gene,]
  cna_tab %>% head()
  cna_tab.m <- melt(cna_tab, id.vars = "gene")
  cna_tab.m %>% head()

  ## input protein data
  if (cancer %in% c("BRCA", "OV", "CO")) {
    pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
    pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")

  } else if (cancer == "UCEC") {
    pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
    pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
  } else if (cancer == "CCRCC") {
    pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
    pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
  } else if (cancer == "LIHC"){
    pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
    pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
  }
  pro_tab <- pro_tab[!is.na(pro_tab$Gene) & pro_tab$Gene == gene,]
  pho_tab <- pho_tab[pho_tab$Gene == gene,]
  
  if (nrow(pro_tab) > 0) {
    pro_tab.m <- melt(pro_tab, id.vars = "Gene")
    pro_tab.m %>% head()
    colnames(pro_tab.m) <- c("Gene", "partID", "exp_value")
    pro_tab.m$Phosphosite <- "PRO"
    sup_tab <- rbind(sup_tab, pro_tab.m[, c("Gene", "Phosphosite", "partID", "exp_value")])
  }
  
  if (nrow(pho_tab) > 0) {
    pho_tab <- pho_tab[!(colnames(pho_tab) %in% c("Peptide_ID"))]
    pho_tab.m <- melt(pho_tab, id.vars = c("Gene", "Phosphosite"))
    pho_tab.m %>% head()
    colnames(pho_tab.m) <- c("Gene", "Phosphosite", "partID", "exp_value")
    pho_tab.m <- pho_tab.m[,c("Gene", "Phosphosite", "partID", "exp_value")]
    sup_tab <- rbind(sup_tab, pho_tab.m[, c("Gene", "Phosphosite", "partID", "exp_value")])
  }

  ## load RNA
  rna_tab <- loadRNA(cancer = cancer)
  rna_tab <- rna_tab[rna_tab$gene == gene,]
  rna_tab.m <- melt(rna_tab, id.vars = "gene")
  colnames(rna_tab.m) <- c("Gene", "partID", "exp_value")
  rna_tab.m$Phosphosite <- "RNA"
  if (nrow(rna_tab.m) > 0) {
    sup_tab <- rbind(sup_tab, rna_tab.m[,c("Gene", "Phosphosite", "partID", "exp_value")])
  }
  
  if (nrow(mut_mat) > 0){
    mut_mat.m <- melt(mut_mat, id.vars = "Hugo_Symbol")
    mut_mat.m %>% head()
    colnames(mut_mat.m) <- c("Gene", "partID", "variant_class")

    ## distinguish by missense and truncation
    tmp <- as.vector(mut_mat.m$variant_class)
    tmp[is.na(tmp)] <- ""
    mut_mat.m$variant_class <- tmp
    mut_mat.m$variant_class_sim <- "other_mutation"
    mut_mat.m$variant_class_sim[mut_mat.m$variant_class == ""] <- "wild_type"
    mut_mat.m$variant_class_sim[mut_mat.m$variant_class  == "Silent"] <- "silent"
    mut_mat.m$variant_class_sim[grepl(x = mut_mat.m$variant_class, pattern = "Missense_Mutation")] <- "missense"
    mut_mat.m$variant_class_sim[grepl(x = mut_mat.m$variant_class, pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del")] <- "truncation"
    mut_mat.m$variant_class_sim[sapply(X = mut_mat.m$variant_class, FUN = function(v) (grepl(pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del", x = v) & grepl(pattern = "Missense_Mutation", x = v)))] <- "missense&truncation"
    
    sup_tab <- merge(sup_tab, mut_mat.m, by.x = c("Gene", "partID"), by.y = c("Gene", "partID"), all.x = T)
    sup_tab %>% tail() 
    sup_tab$variant_class[is.na(sup_tab$variant_class)] <- ""
    sup_tab$variant_class_sim[is.na(sup_tab$variant_class_sim)] <- "wild_type"
  } else {
    stop("no mutation!")
    sup_tab$variant_class_sim <- NA
  }
  
  if (nrow(cna_tab.m) > 0) {
    colnames(cna_tab.m) <- c("gene", "partID", "CNA")
    cna_tab.m$is.deleted <- (cna_tab.m$CNA == "deletion")
    sup_tab <- merge(sup_tab, cna_tab.m[, c("gene", "partID", "is.deleted")], all.x = T)
    head(sup_tab)
    sup_tab$is.deleted[is.na(sup_tab$is.deleted)] <- FALSE
  } else {
    stop("no CNA!")
    sup_tab$is.deleted <- FALSE
  }
  
  if (nrow(vaf_part) > 0) {
    sup_tab <- merge(sup_tab, vaf_part[, c("partID", "VAF")], by = c("partID"), all.x = T)
    head(sup_tab)
  } else {
    stop("no VAF!")
    sup_tab$VAF <- NA
  }
  
  head(sup_tab)
  
  ## bind with super table
  sup_tab$cancer <- cancer
  
  ## delete phosphosites with too much NAs
  sup_tab <- sup_tab[sup_tab$Phosphosite %in% c("RNA", "PRO", "S315", "S392"),]
  # sup_tab <- sup_tab[sup_tab$Phosphosite %in% c("PRO", "S315"),]
  expression2plot <- c("RNA", "PRO", "S315", "S392")
  
  # Plot heatmap ------------------------------------------------------------
  sup_tab$id_row <- paste0(sup_tab$Gene, "_", sup_tab$Phosphosite)
  cap <- 3
  sup_tab$exp_value_capped <- sup_tab$exp_value
  sup_tab$exp_value_capped[sup_tab$exp_value > cap] <- cap
  sup_tab$exp_value_capped[sup_tab$exp_value < (-cap)] <- (-cap)
  breaks = seq(-(cap),cap, by=0.2)
  
  ## make the matrix for the heatmap body
  # df_value <- dcast(data = sup_tab, id_row ~ partID, value.var = "exp_value_capped")
  df_value <- dcast(data = sup_tab, id_row ~ partID, value.var = "exp_value")
  
  df_value %>% head()
  mat_value <- as.matrix(df_value[,-1])
  rownames(mat_value) <- df_value$id_row
  head(mat_value)
  # mat_value <- mat_value[,colSums(!is.na(mat_value)) > 0]
  mat_value %>% head()
  
  ## make the annotation columns
  col_anno <- sup_tab[sup_tab$Phosphosite == "PRO", c("partID", "variant_class_sim", "is.deleted", "VAF", "exp_value")]
  col_anno %>% head()
  col_anno <- col_anno[order(col_anno$exp_value, decreasing = T),]
  col_anno <- col_anno[order(col_anno$is.deleted, decreasing = T),]
  col_anno <- col_anno[order(col_anno$variant_class_sim),]
  col_anno %>% head()
  rownames(col_anno) <- col_anno$partID
  col_anno$partID <- NULL
  col_anno$cancer <- NULL
  col_anno$exp_value <- NULL
  col_anno$is.deleted <- as.character(col_anno$is.deleted)
  colnames(col_anno) <- c("mutation", "deletion", "VAF")
  
  ## order the matrix by variant class and RNA-PRO-Phosphosite
  mat_value <- mat_value[, unique(rownames(col_anno))]
  mat_value <- mat_value[intersect(paste0(gene, "_", expression2plot), rownames(mat_value)),]
  
  ## add color palette
  color.palette <- colorRampPalette(rev(brewer.pal(10, "RdBu")))(length(breaks))
  ann_colors = list(mutation = c(missense = "#E41A1C", truncation = "#377EB8", wild_type = "#33A02C", "missense&truncation" = "#6A3D9A", other_mutation = "#FF7F00", silent = "#B2DF8A"),
                    deletion = c("TRUE" = "black", "FALSE" = "white"))
  vaf_ann_breaks <- seq(0.1, 0.8, 0.1)
  vaf_ann_color <- c(col_paletteR(length(vaf_ann_breaks)))
  names(vaf_ann_color) <- vaf_ann_breaks
  ann_colors[["VAF"]] <- vaf_ann_color
  
  ## actual plotting
  if (cancer == "LIHC") {
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno,
                           scale = "row", legend = FALSE, annotation_legend = F,
                           na_col = "white", 
                           cellwidth = 4,
                           cellheight = 9,
                           cluster_rows=F, cluster_cols=F, show_colnames = F, annotation_colors = ann_colors)
    save_pheatmap_pdf(x = my_heatmap, 
                      filename = paste0(makeOutDir(resultD = resultD), "TP53_", cancer, ".pdf"), 
                      width = 12, height = 1.5)
    
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno,
                           scale = "row", legend = FALSE, annotation_legend = T,
                           na_col = "white", 
                           cellwidth = 4,
                           cellheight = 9,
                           cluster_rows=F, cluster_cols=F, show_colnames = F, annotation_colors = ann_colors)
    save_pheatmap_pdf(x = my_heatmap, 
                      filename = paste0(makeOutDir(resultD = resultD), "TP53_", cancer, "_with_legend.pdf"), 
                      width = 12, height = 3)
  } else {
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno,
                           scale = "row", legend = FALSE, annotation_legend = F,
                           na_col = "white", 
                           cluster_rows=F, cluster_cols=F, show_colnames = F, annotation_colors = ann_colors)
    save_pheatmap_pdf(x = my_heatmap, 
                      filename = paste0(makeOutDir(resultD = resultD), "TP53_", cancer, ".pdf"), 
                      width = 10, height = 1.2)
    
    my_heatmap <- pheatmap(mat_value, 
                           color = color.palette,
                           annotation_col = col_anno,
                           scale = "row", legend = FALSE, annotation_legend = T,
                           na_col = "white", 
                           cluster_rows=F, cluster_cols=F, show_colnames = F, annotation_colors = ann_colors)
    save_pheatmap_pdf(x = my_heatmap, 
                      filename = paste0(makeOutDir(resultD = resultD), "TP53_", cancer, "_with_legend.pdf"), 
                      width = 8, height = 3)
  }
}
