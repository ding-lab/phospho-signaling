# Yige Wu @ WashU Feb 2018
# check how particular phosphosite abundance change after normalization


# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/get_tumor_normal_pair.R')
clinical <- read.delim(paste0(inputD, "Specimen_Data_20161005_Yige_20180127.txt"), sep = "\t")
clinical <- data.frame(clinical)
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(ggplot2)
subtype_map <- read_excel("~/Box Sync/cptac2p/cptac_shared/5_CPTAC2_Breast_Prospective_Collection_BI/proteome-data-v1.01011/data/20171106_CPTAC2_ProspectiveBreastCancer_Sample_Annotations_Table_v79.xlsx")


cancer <- "BRCA"
Pro <- fread(raw[cancer, "protein"], data.table=FALSE)

# make density plot for phsphoprotein data for normals and tumor subtypes ------------------------
Pho <- fread(raw[cancer, "phosphoprotein"], data.table=FALSE)
Pho.f2p = unshared_pho(Pho, Pho[,c("Gene","Phosphosite")])
Pho_head <- formatPhosphosite(phosphosite_vector = Pho$Phosphosite, gene_vector = Pho$Gene)
Pho.n <- get_normal(expression = Pho.f2p, clinical.m = clinical)
Pho.t <- get_tumor(expression = Pho.f2p, clinical.m = clinical)
Pho.scale.n = normalize_by_sample(Pho.n)
Pho.scale.t = normalize_by_sample(Pho.t)

## get PAM50 subtypes
samples.t <- colnames(Pho.t)
subtypes.t <- sampID2pam50(sampleID_vector = samples.t, subtype_map = subtype_map)
samples.n <- colnames(Pho.n)

## melt phosphoprotein data
Pho <- cbind(Pho_head, Pho.n, Pho.t)
# Pho2m <- Pho[rowSums(!is.na(Pho)) > (ncol(Pho)*(2/3)),]
Pho2m <- Pho[sample(1:nrow(Pho), 2000),]
Pho.m <- melt(Pho2m)

Pho.m <- merge(Pho.m, rbind(data.frame(variable = samples.t,
                                       subtype = subtypes.t),
                            data.frame(variable = samples.n,
                                       subtype = "normal")), 
              all.x = T)
Pho.m$tumor_normal <- ifelse(Pho.m$variable %in% samples.t, "tumor", "normal")
Pho.m$tumor_normal <- factor(Pho.m$tumor_normal, levels = c("tumor", "normal"))
Pho.m$method <- "CDAP_original"
Pho.scale <- cbind(Pho_head, Pho.scale.n, Pho.scale.t)
# Pho.scale2m <- Pho.scale[rowSums(!is.na(Pho.scale)) > (ncol(Pho.scale)*(2/3)),]
Pho.scale2m <- Pho.scale[sample(1:nrow(Pho.scale), 2000),]
Pho.scale.m <- melt(Pho.scale2m)

Pho.scale.m <- merge(Pho.scale.m, rbind(data.frame(variable = samples.t,
                                                   subtype = subtypes.t),
                                        data.frame(variable = samples.n,
                                                   subtype = "normal")), 
                     all.x = T)
Pho.scale.m$tumor_normal <- ifelse(Pho.scale.m$variable %in% samples.t, "tumor", "normal")
Pho.scale.m$tumor_normal <- factor(Pho.scale.m$tumor_normal, levels = c("tumor", "normal"))
Pho.scale.m$method <- "CDAP_scaled"

# process broad un-normalized data ----------------------------------------
## input broad data
Pho.broad.scale <- fread(input = "cptac2p/cptac_shared/5_CPTAC2_Breast_Prospective_Collection_BI/phosphoproteome-data-v1.02011/normalized-data/phosphoproteome-ratio-norm.gct", 
                         data.table = F)
Pho.broad <- fread(input = "cptac2p/cptac_shared/5_CPTAC2_Breast_Prospective_Collection_BI/phosphoproteome-data-v1.02011/parsed-data/phosphoproteome-ratio.gct", 
                   data.table = F)

## format to numeric
row_start <- which(Pho.broad[,2] != "na")[1]
col_start <- which(Pho.broad[2,2:ncol(Pho.broad)] != "na")[1]+1

## get tumor/normal and subtype info
Pho.broad_head <- data.frame(SUBSTRATE = Pho.broad[row_start:nrow(Pho.broad), "geneSymbol"],
                             SUB_MOD_RSD = Pho.broad[row_start:nrow(Pho.broad), 1 ])
tumor_normal <- as.character(Pho.broad[Pho.broad[,1] == "Sample.Type", col_start:ncol(Pho.broad)])
pam50.broad <- as.character(Pho.broad[Pho.broad[,1] == "PAM50", col_start:ncol(Pho.broad)])

## melt
Pho.broad2m <- Pho.broad[row_start:nrow(Pho.broad), col_start:ncol(Pho.broad)]
Pho.broad2m <- data.frame(as.matrix(Pho.broad2m))
partIDs.broad <- colnames(Pho.broad2m)
Pho.broad2m <- cbind(Pho.broad_head, Pho.broad2m)
# Pho.broad2m <- Pho.broad2m[rowSums(!is.na(Pho.broad2m)) > ncol(Pho.broad2m)*(2/3),]
Pho.broad2m <- Pho.broad2m[sample(1:nrow(Pho.broad2m), 1000),]
Pho.broad.m <- melt(Pho.broad2m, id.vars = c("SUBSTRATE", "SUB_MOD_RSD"))

## annotate with tumor/normal and subtype info
subtype_map.broad <- rbind(data.frame(variable = partIDs.broad[tumor_normal == "Tumor"],
                                      subtype = pam50.broad[tumor_normal == "Tumor"],
                                      tumor_normal = "tumor"),
                           data.frame(variable = partIDs.broad[tumor_normal == "Adjacent_Normal"],
                                      subtype = "normal",
                                      tumor_normal = "normal"))

Pho.broad.m <- merge(Pho.broad.m, subtype_map.broad, by = c("variable"), all.x = T)
Pho.broad.m$method <- "broad_original"
Pho.broad.m$value <- as.numeric(as.vector(Pho.broad.m$value))
# process broad normalized data ----------------------------------------


## format to numeric
row_start <- which(Pho.broad.scale[,2] != "na")[1]
col_start <- which(Pho.broad.scale[2,2:ncol(Pho.broad.scale)] != "na")[1]+1

## get tumor/normal and subtype info
Pho.broad.scale_head <- data.frame(SUBSTRATE = Pho.broad.scale[row_start:nrow(Pho.broad.scale), "geneSymbol"],
                                   SUB_MOD_RSD = Pho.broad.scale[row_start:nrow(Pho.broad.scale), 1 ])
tumor_normal <- as.character(Pho.broad.scale[Pho.broad.scale[,1] == "Sample.Type", col_start:ncol(Pho.broad.scale)])
pam50.broad <- as.character(Pho.broad.scale[Pho.broad.scale[,1] == "PAM50", col_start:ncol(Pho.broad.scale)])

## melt
Pho.broad.scale2m <- Pho.broad.scale[row_start:nrow(Pho.broad.scale), col_start:ncol(Pho.broad.scale)]
Pho.broad.scale2m <- data.frame(as.matrix(Pho.broad.scale2m))
partIDs.broad <- colnames(Pho.broad.scale2m)
Pho.broad.scale2m <- cbind(Pho.broad.scale_head, Pho.broad.scale2m)
# Pho.broad.scale2m <- Pho.broad.scale2m[rowSums(!is.na(Pho.broad.scale2m)) > ncol(Pho.broad.scale2m)*(2/3),]
Pho.broad.scale2m <- Pho.broad.scale2m[sample(1:nrow(Pho.broad.scale2m), 1000),]
Pho.broad.scale.m <- melt(Pho.broad.scale2m, id.vars = c("SUBSTRATE", "SUB_MOD_RSD"))

## annotate with tumor/normal and subtype info
subtype_map.broad <- rbind(data.frame(variable = partIDs.broad[tumor_normal == "Tumor"],
                                      subtype = pam50.broad[tumor_normal == "Tumor"],
                                      tumor_normal = "tumor"),
                           data.frame(variable = partIDs.broad[tumor_normal == "Adjacent_Normal"],
                                      subtype = "normal",
                                      tumor_normal = "normal"))

Pho.broad.scale.m <- merge(Pho.broad.scale.m, subtype_map.broad, by = c("variable"), all.x = T)
Pho.broad.scale.m$method <- "broad_normalized"
Pho.broad.scale.m$value <- as.numeric(as.vector(Pho.broad.scale.m$value))

# merge broad data with CDAP ----------------------------------------------
col2merge <- intersect(colnames(Pho.broad.m), colnames(Pho.m))
df <- rbind(Pho.m[,col2merge], Pho.scale.m[,col2merge])


# density plot ------------------------------------------------------------
p <- ggplot(data = df)
p <- p + geom_density(mapping = aes(x = value, group = method, color = method, fill = method), alpha = 0.4)
p <- p + facet_grid(tumor_normal + subtype~.)
p <- p + theme_bw()
p
resultDnow <- makeOutDir()
ggsave(filename = paste0(resultDnow, cancer, "_phosphosite_abundance_normalized_density_CDAP_partial.pdf"), height = 6, width = 6)
# ggsave(filename = paste0(resultDnow, cancer, "_phosphosite_abundance_normalized_density.pdf"), height = 6, width = 6)

# df <- rbind(Pho.broad.m[,col2merge], Pho.broad.scale.m[,col2merge])
df <- rbind(Pho.broad.m, Pho.broad.scale.m)
df <- df[!is.na(df$tumor_normal),]
resultDnow <- makeOutDir()
p <- ggplot(data = df)
p <- p + geom_density(mapping = aes(x = value, group = method, color = method, fill = method), alpha = 0.4)
p <- p + facet_grid(tumor_normal + subtype~.)
p <- p + theme_bw()
p
ggsave(filename = paste0(resultDnow, cancer, "_phosphosite_abundance_normalized_density_broad_partial.pdf"), 
       height = 6, width = 6)


# compare column mean and median between tumors and normals ---------------
df_center <- rbind(data.frame(tumor_normal = "normal",
                              value = colMeans(Pho.n, na.rm = T),
                              type = "colMean"),
                   data.frame(tumor_normal = "tumor",
                              value = colMeans(Pho.t, na.rm = T),
                              type = "colMean"),
                   data.frame(tumor_normal = "normal",
                              value = sapply(1:ncol(Pho.n), FUN = function(n, m) median(m[,n], na.rm = T), m = Pho.n),
                              type = "colMedian"),
                   data.frame(tumor_normal = "tumor",
                              value = sapply(1:ncol(Pho.t), FUN = function(n, m) median(m[,n], na.rm = T), m = Pho.t),
                              type = "colMedian"))

p <- ggplot(data = df_center)
p <- p + geom_boxplot(mapping = aes(x = type, y = value, color = tumor_normal))
p <- p + theme_bw()
# p <- p + ggtitle(label = paste0(gene, " ", site, " phosphosite abundance normalized values"))
p = p + theme(axis.text.x = element_text(angle = -15))
p
resultDnow <- makeOutDir()
ggsave(filename = paste0(resultDnow, cancer, "_phosphosite_abundance_colMean_colMedian_tumorORnormal.pdf"), height = 4, width = 5)

# make density plot for CDAP protein data for normals and tumor subtypes ------------------------
Pro_raw <- fread(raw[cancer, "protein"], data.table=FALSE)
Pro.f <- unshared_pro(Pro_raw)
Pro.f2p = data.frame(Gene = rownames(Pro.f))
Pro.f2p <- cbind(Pro.f2p, Pro.f)
Pro_head <- data.frame(Gene = Pro.f2p$Gene)
Pro.n <- get_normal(expression = Pro.f2p, clinical.m = clinical)
Pro.t <- get_tumor(expression = Pro.f2p, clinical.m = clinical)
Pro.scale.n = normalize_by_sample(Pro.n)
Pro.scale.t = normalize_by_sample(Pro.t)

## get PAM50 subtypes
samples.t <- colnames(Pro.t)
subtypes.t <- sampID2pam50(sampleID_vector = samples.t, subtype_map = subtype_map)
samples.n <- colnames(Pro.n)

## melt ProsProprotein data
Pro <- cbind(Pro_head, Pro.n, Pro.t)
Pro2m <- Pro[rowSums(!is.na(Pro)) > (ncol(Pro)/3),]
Pro.m <- melt(Pro2m)

Pro.m <- merge(Pro.m, rbind(data.frame(variable = samples.t,
                                       subtype = subtypes.t),
                            data.frame(variable = samples.n,
                                       subtype = "normal")), 
               all.x = T)
Pro.m$tumor_normal <- ifelse(Pro.m$variable %in% samples.t, "tumor", "normal")
Pro.m$tumor_normal <- factor(Pro.m$tumor_normal, levels = c("tumor", "normal"))
Pro.m$method <- "CDAP_original"
Pro.scale <- cbind(Pro_head, Pro.scale.n, Pro.scale.t)
Pro.scale2m <- Pro.scale[rowSums(!is.na(Pro.scale)) > (ncol(Pro.scale)/3),]
Pro.scale.m <- melt(Pro.scale2m)

Pro.scale.m <- merge(Pro.scale.m, rbind(data.frame(variable = samples.t,
                                                   subtype = subtypes.t),
                                        data.frame(variable = samples.n,
                                                   subtype = "normal")), 
                     all.x = T)
Pro.scale.m$tumor_normal <- ifelse(Pro.scale.m$variable %in% samples.t, "tumor", "normal")
Pro.scale.m$tumor_normal <- factor(Pro.scale.m$tumor_normal, levels = c("tumor", "normal"))
Pro.scale.m$method <- "CDAP_scaled"
df <- rbind(Pro.m, Pro.scale.m)
## density plot
p <- ggplot(data = df)
p <- p + geom_density(mapping = aes(x = value, group = method, color = method, fill = method), alpha = 0.5)
p <- p + facet_grid(tumor_normal + subtype~.)
p <- p + theme_bw()
p
resultDnow <- makeOutDir()
ggsave(filename = paste0(resultDnow, cancer, "_Protein_abundance_normalized_density.pdf"), height = 6, width = 6)


# make density plot for Broad Protein data for normals and tumor subtypes ------------------------
## melt Pro.broadtein data
Pro.broad.m <- melt(Pro.broad2m[sample(1:nrow(Pro.broad2m), 1000),], id.vars = "SUBSTRATE")

Pro.broad.m <- merge(Pro.broad.m, rbind(data.frame(variable = colnames(Pro.broad2m)[-1][tumor_normal == "Tumor"],
                                                   tumor_normal = "tumor",
                                                   subtype = pam50.broad[tumor_normal == "Tumor"]),
                                        data.frame(variable = colnames(Pro.broad2m)[-1][tumor_normal == "Adjacent_Normal"],
                                                   tumor_normal = "normal",
                                                   subtype = "normal")),
                     all.x = T)

Pro.broad.m$method <- "Broad_original"
Pro.broad.m$value <- as.numeric(as.vector(Pro.broad.m$value))

Pro.broad.scale.m <- melt(Pro.broad.scale2m[sample(1:nrow(Pro.broad.scale2m), 1000),], id.vars = "SUBSTRATE")
Pro.broad.scale.m <- merge(Pro.broad.scale.m, rbind(data.frame(variable = colnames(Pro.broad2m)[-1][tumor_normal == "Tumor"],
                                                               tumor_normal = "tumor",
                                                               subtype = pam50.broad[tumor_normal == "Tumor"]),
                                                    data.frame(variable = colnames(Pro.broad2m)[-1][tumor_normal == "Adjacent_Normal"],
                                                               tumor_normal = "normal",
                                                               subtype = "normal")),
                           all.x = T)
Pro.broad.scale.m$method <- "Broad_normalized"
Pro.broad.scale.m$value <- as.numeric(as.vector(Pro.broad.scale.m$value))

df <- rbind(Pro.broad.m, Pro.broad.scale.m)
df <- df[!is.na(df$tumor_normal),]
## density plot
resultDnow <- makeOutDir()
p <- ggplot(data = df)
p <- p + geom_density(mapping = aes(x = value, group = method, color = method, fill = method), alpha = 0.5)
p <- p + facet_grid(tumor_normal + subtype~.)
p <- p + theme_bw()
p
ggsave(filename = paste0(resultDnow, "Broad", "_Protein_abundance_normalized_density.pdf"), height = 6, width = 6)





# plot EGFR phosphosite abundance for CDAP data --------------------------------
gene <- "EGFR"
resultDnow <- makeOutDir()

## try normalizeQuantiles
library(limma)
Pho.nq <- normalizeQuantiles(as.matrix(cbind(Pho.n, Pho.t)))
Pho.nq.n <- Pho.nq[, colnames(Pho.n)]
Pho.nq.t <- Pho.nq[, colnames(Pho.t)]

## try correcting for proportional of each gene in the mixture

for (site in Pho_head$SUB_MOD_RSD[Pho_head$SUBSTRATE == gene]) {
  df_orig <- rbind(data.frame(tumor_normal = "normal",
                              phospho = as.numeric(Pho.n[Pho_head$SUBSTRATE == gene & Pho_head$SUB_MOD_RSD == site,]),
                              transcript = unique(Pho_head$transcript[Pho_head$SUBSTRATE == gene & Pho_head$SUB_MOD_RSD == site])),
                   data.frame(tumor_normal = "tumor",
                              phospho = as.numeric(Pho.t[Pho_head$SUBSTRATE == gene & Pho_head$SUB_MOD_RSD == site,]),
                              transcript = unique(Pho_head$transcript[Pho_head$SUBSTRATE == gene & Pho_head$SUB_MOD_RSD == site])))
  df_orig$method <- "CDAP_original"
  
  values.n <- as.numeric(Pho.scale.n[Pho_head$SUBSTRATE == gene & Pho_head$SUB_MOD_RSD == site,])
  if (length(values.n) >0 ) {
    df_scale <- rbind(data.frame(tumor_normal = "normal",
                                 phospho = values.n,
                                 transcript = unique(Pho_head$transcript[Pho_head$SUBSTRATE == gene & Pho_head$SUB_MOD_RSD == site])),
                      data.frame(tumor_normal = "tumor",
                                 phospho = as.numeric(Pho.scale.t[Pho_head$SUBSTRATE == gene & Pho_head$SUB_MOD_RSD == site,]),
                                 transcript = unique(Pho_head$transcript[Pho_head$SUBSTRATE == gene & Pho_head$SUB_MOD_RSD == site])))
    
    df_scale$method <- "CDAP_scaled_by_colMean&SD"
    
    df_nq <- rbind(data.frame(tumor_normal = "normal",
                              phospho = as.numeric(Pho.nq.n[Pho_head$SUBSTRATE == gene & Pho_head$SUB_MOD_RSD == site,]),
                              transcript = unique(Pho_head$transcript[Pho_head$SUBSTRATE == gene & Pho_head$SUB_MOD_RSD == site])),
                   data.frame(tumor_normal = "tumor",
                              phospho = as.numeric(Pho.nq.t[Pho_head$SUBSTRATE == gene & Pho_head$SUB_MOD_RSD == site,]),
                              transcript = unique(Pho_head$transcript[Pho_head$SUBSTRATE == gene & Pho_head$SUB_MOD_RSD == site])))
    
    
    df_nq$method <- "CDAP_normalizeQuantiles"
    
    cf <- sum(Pho.t[Pho_head$SUBSTRATE == gene & Pho_head$SUB_MOD_RSD == site,], na.rm = T)*nrow(Pho.t)/sum(rowSums(Pho.t, na.rm = T), na.rm = T)
    
    df_nq <- rbind(data.frame(tumor_normal = "normal",
                              phospho = as.numeric(Pho.nq.n[Pho_head$SUBSTRATE == gene & Pho_head$SUB_MOD_RSD == site,]),
                              transcript = unique(Pho_head$transcript[Pho_head$SUBSTRATE == gene & Pho_head$SUB_MOD_RSD == site])),
                   data.frame(tumor_normal = "tumor",
                              phospho = as.numeric(Pho.nq.t[Pho_head$SUBSTRATE == gene & Pho_head$SUB_MOD_RSD == site,]),
                              transcript = unique(Pho_head$transcript[Pho_head$SUBSTRATE == gene & Pho_head$SUB_MOD_RSD == site])))
    
    
    df_nq$method <- "CDAP_normalizeQuantiles"
    
    df <- rbind(df_orig, df_scale, df_nq)
    # df <- rbind(df_orig, df_scale)
    df$tumor_normal <- factor(df$tumor_normal, levels = c("tumor", "normal"))
    df$method <- factor(df$method, levels = c("CDAP_original", "CDAP_scaled_by_colMean&SD", "CDAP_normalizeQuantiles"))
    
    p <- ggplot(data = df)
    p <- p + geom_boxplot(mapping = aes(x = method, y = phospho, color = tumor_normal))
    p <- p + theme_bw()
    p <- p + ggtitle(label = paste0(gene, " phosphosite abundance"), subtitle = paste0(unique(df$transcript), "_", site))
    p <- p + ylab("phosphorylation abundance(log2 ratio)")
    p = p + theme(axis.text.x = element_text(angle = -15)) + theme(plot.title = element_text(size = 10, face = "bold"))
    p
    ggsave(filename = paste0(resultDnow, "CDAP_breast_", gene, "_", unique(df$transcript), "_", site, "_phosphosite_abundance_normalized.pdf"), height = 4, width = 5)
    
  }

}

# plot EGFR protein abundance for CDAP data --------------------------------
df_orig <- rbind(data.frame(tumor_normal = "normal",
                            protein = as.numeric(Pro.n[Pro_head$Gene == gene,])),
                 data.frame(tumor_normal = "tumor",
                            protein = as.numeric(Pro.t[Pro_head$Gene == gene,])))
df_orig$method <- "CDAP_original"

values.n <- as.numeric(Pro.scale.n[Pro_head$Gene == gene,])
if (length(values.n) >0 ) {
  df_scale <- rbind(data.frame(tumor_normal = "normal",
                               protein = values.n),
                    data.frame(tumor_normal = "tumor",
                               protein = as.numeric(Pro.scale.t[Pro_head$Gene == gene,])))
  
  df_scale$method <- "CDAP_scaled_by_colMean&SD"
  
  # df <- rbind(df_orig, df_sweepmed, df_scale)
  df <- rbind(df_orig, df_scale)
  df$tumor_normal <- factor(df$tumor_normal, levels = c("tumor", "normal"))
  p <- ggplot(data = df)
  p <- p + geom_boxplot(mapping = aes(x = method, y = protein, color = tumor_normal))
  p <- p + theme_bw()
  p <- p + ggtitle(label = paste0(gene, " protein abundance"))
  p <- p + ylab("Protein abundance(log2 ratio)")
  p = p + theme(axis.text.x = element_text(angle = -15)) + theme(plot.title = element_text(size = 10, face = "bold"))
  p
  ggsave(filename = paste0(resultDnow, "CDAP_breast_", gene, "_", unique(df$transcript), "_", "Protein", "_abundance_normalized.pdf"), height = 4, width = 5)
}



# plot EGFR phosphosite abundance in broad dataset---------------------------------------------------------------
row_start <- which(Pho.broad[,2] != "na")[1]
col_start <- which(Pho.broad[2,2:ncol(Pho.broad)] != "na")[1]+1
Pho.broad_head <- data.frame(SUBSTRATE = Pho.broad[row_start:nrow(Pho.broad), "geneSymbol"],
                             SUB_MOD_RSD = Pho.broad[row_start:nrow(Pho.broad), 1 ])
Pho.broad2m <- Pho.broad[row_start:nrow(Pho.broad), col_start:ncol(Pho.broad)]
Pho.broad2m <- data.frame(as.matrix(Pho.broad2m))
partIDs.broad <- colnames(Pho.broad2m)
Pho.broad2m <- cbind(Pho.broad_head, Pho.broad2m)

row_start <- which(Pho.broad.scale[,2] != "na")[1]
col_start <- which(Pho.broad.scale[2,2:ncol(Pho.broad.scale)] != "na")[1]+1

## get tumor/normal and subtype info
Pho.broad.scale_head <- data.frame(SUBSTRATE = Pho.broad.scale[row_start:nrow(Pho.broad.scale), "geneSymbol"],
                                   SUB_MOD_RSD = Pho.broad.scale[row_start:nrow(Pho.broad.scale), 1 ])
## melt
Pho.broad.scale2m <- Pho.broad.scale[row_start:nrow(Pho.broad.scale), col_start:ncol(Pho.broad.scale)]
Pho.broad.scale2m <- data.frame(as.matrix(Pho.broad.scale2m))
partIDs.broad <- colnames(Pho.broad.scale2m)
Pho.broad.scale2m <- cbind(Pho.broad.scale_head, Pho.broad.scale2m)


gene <- "EGFR"
Pho.broad2m_gene <- Pho.broad2m[Pho.broad2m$SUBSTRATE == "EGFR",]
Pho.broad.scale2m_gene <- Pho.broad.scale2m[Pho.broad.scale2m$SUBSTRATE == "EGFR",]
# site <- "T693"

for (site in Pho.broad2m_gene$SUB_MOD_RSD) {
  values.t <- as.vector(Pho.broad2m_gene[grepl(pattern = site, x = Pho.broad2m_gene$SUB_MOD_RSD, ignore.case = T), colnames(Pho.broad2m_gene)[which(tumor_normal == "Tumor") + 2]])
  values.t <- sapply(1:length(values.t), FUN = function(n, v) as.numeric(as.vector(v[1,n])), v = values.t)
  values.n <- as.vector(Pho.broad2m_gene[grepl(pattern = site, x = Pho.broad2m_gene$SUB_MOD_RSD, ignore.case = T), colnames(Pho.broad2m_gene)[which(tumor_normal == "Adjacent_Normal") + 2]])
  values.n <- sapply(1:length(values.n), FUN = function(n, v) as.numeric(as.vector(v[1,n])), v = values.n)
  
  df_orig <- rbind(data.frame(tumor_normal = "tumor",
                              phospho = values.t),
                   data.frame(tumor_normal = "normal",
                              phospho = values.n))
  df_orig$method <- "broad_original"
  
  values.t <- as.vector(Pho.broad.scale2m_gene[grepl(pattern = site, x = Pho.broad.scale2m_gene$SUB_MOD_RSD, ignore.case = T), colnames(Pho.broad.scale2m_gene)[which(tumor_normal == "Tumor") + 2]])
  values.t <- sapply(1:length(values.t), FUN = function(n, v) as.numeric(as.vector(v[1,n])), v = values.t)
  values.n <- as.vector(Pho.broad.scale2m_gene[grepl(pattern = site, x = Pho.broad.scale2m_gene$SUB_MOD_RSD, ignore.case = T), colnames(Pho.broad.scale2m_gene)[which(tumor_normal == "Adjacent_Normal") + 2]])
  values.n <- sapply(1:length(values.n), FUN = function(n, v) as.numeric(as.vector(v[1,n])), v = values.n)
  
  df_scale <- rbind(data.frame(tumor_normal = "tumor",
                               phospho = values.t),
                    data.frame(tumor_normal = "normal",
                               phospho = values.n))
  df_scale$method <- "broad_normalized"
  
  df <- rbind(df_orig, df_scale)
  df$method <- factor(df$method, levels = c("broad_original", "broad_normalized"))
  resultDnow <- makeOutDir()
  p <- ggplot(data = df)
  p <- p + geom_boxplot(mapping = aes(x = method, y = phospho, color = tumor_normal))
  p <- p + theme_bw()
  p <- p + ggtitle(label = paste0(gene, " phosphosite abundance"), subtitle = site)
  p <- p + ylab("phosphorylation abundance(log2 ratio)")
  p = p + theme(axis.text.x = element_text(angle = -15)) + theme(plot.title = element_text(size = 10, face = "bold"))
  p
  ggsave(filename = paste0(resultDnow, "Broad_", gene, "_", site, "_phosphosite_abundance_normalized.pdf"), height = 4, width = 5)
}

# plot EGFR protein abundance in broad dataset---------------------------------------------------------------
## input broad protein data
Pro.broad.scale <- fread(input = "cptac2p/cptac_shared/5_CPTAC2_Breast_Prospective_Collection_BI/proteome-data-v1.01011/normalized-data/proteome-ratio-norm.gct", 
                         data.table = F)
Pro.broad <- fread(input = "cptac2p/cptac_shared/5_CPTAC2_Breast_Prospective_Collection_BI/proteome-data-v1.01011/parsed-data/proteome-ratio.gct", 
                   data.table = F)

row_start <- which(Pro.broad[,2] != "na")[1]
col_start <- which(Pro.broad[2,2:ncol(Pro.broad)] != "na")[1]+1
Pro.broad_head <- data.frame(SUBSTRATE = Pro.broad[row_start:nrow(Pro.broad), "geneSymbol"])
Pro.broad2m <- Pro.broad[row_start:nrow(Pro.broad), col_start:ncol(Pro.broad)]
Pro.broad2m <- data.frame(as.matrix(Pro.broad2m))
partIDs.broad <- colnames(Pro.broad2m)
Pro.broad2m <- cbind(Pro.broad_head, Pro.broad2m)

row_start <- which(Pro.broad.scale[,2] != "na")[1]
col_start <- which(Pro.broad.scale[2,2:ncol(Pro.broad.scale)] != "na")[1]+1

## get tumor/normal and subtype info
Pro.broad.scale_head <- data.frame(SUBSTRATE = Pro.broad.scale[row_start:nrow(Pro.broad.scale), "geneSymbol"])
## melt
Pro.broad.scale2m <- Pro.broad.scale[row_start:nrow(Pro.broad.scale), col_start:ncol(Pro.broad.scale)]
Pro.broad.scale2m <- data.frame(as.matrix(Pro.broad.scale2m))
partIDs.broad <- colnames(Pro.broad.scale2m)
Pro.broad.scale2m <- cbind(Pro.broad.scale_head, Pro.broad.scale2m)
tumor_normal <- as.character(Pho.broad[Pho.broad[,1] == "Sample.Type", col_start:ncol(Pho.broad)])


gene <- "EGFR"
Pro.broad2m_gene <- Pro.broad2m[Pro.broad2m$SUBSTRATE == "EGFR",]
Pro.broad.scale2m_gene <- Pro.broad.scale2m[Pro.broad.scale2m$SUBSTRATE == "EGFR",]
# site <- "T693"

values.t <- as.vector(Pro.broad2m_gene[, colnames(Pro.broad2m_gene)[which(tumor_normal == "Tumor") + 1]])
values.t <- sapply(1:length(values.t), FUN = function(n, v) as.numeric(as.vector(v[1,n])), v = values.t)
values.n <- as.vector(Pro.broad2m_gene[, colnames(Pro.broad2m_gene)[which(tumor_normal == "Adjacent_Normal") + 1]])
values.n <- sapply(1:length(values.n), FUN = function(n, v) as.numeric(as.vector(v[1,n])), v = values.n)

df_orig <- rbind(data.frame(tumor_normal = "tumor",
                            phospho = values.t),
                 data.frame(tumor_normal = "normal",
                            phospho = values.n))
df_orig$method <- "broad_original"

values.t <- as.vector(Pro.broad.scale2m_gene[, colnames(Pro.broad.scale2m_gene)[which(tumor_normal == "Tumor") + 1]])
values.t <- sapply(1:length(values.t), FUN = function(n, v) as.numeric(as.vector(v[1,n])), v = values.t)
values.n <- as.vector(Pro.broad.scale2m_gene[, colnames(Pro.broad.scale2m_gene)[which(tumor_normal == "Adjacent_Normal") + 1]])
values.n <- sapply(1:length(values.n), FUN = function(n, v) as.numeric(as.vector(v[1,n])), v = values.n)

df_scale <- rbind(data.frame(tumor_normal = "tumor",
                             phospho = values.t),
                  data.frame(tumor_normal = "normal",
                             phospho = values.n))
df_scale$method <- "broad_normalized"

df <- rbind(df_orig, df_scale)
df$method <- factor(df$method, levels = c("broad_original", "broad_normalized"))
resultDnow <- makeOutDir()
p <- ggplot(data = df)
p <- p + geom_boxplot(mapping = aes(x = method, y = phospho, color = tumor_normal))
p <- p + theme_bw()
p <- p + ggtitle(label = paste0(gene, " protein abundance"))
p <- p + ylab("protein abundance(log2 ratio)")
p = p + theme(axis.text.x = element_text(angle = -15)) + theme(plot.title = element_text(size = 10, face = "bold"))
p
ggsave(filename = paste0(resultDnow, "Broad_", gene, "_protein_abundance_normalized.pdf"), height = 4, width = 5)
