# Yige Wu @ WashU 2018 Feb
# make heatmap showing kinase/phosphotase pairs with consistent phosphoprotein expression in a cascade with substrate in different pathways


# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
clinical <- read.delim(paste0(inputD, "Specimen_Data_20161005_Yige_20180127.txt"), sep = "\t")
clinical <- data.frame(clinical)

# input -----------------------------------------------------------------
us_diffexp <- fread(input = "cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/union_diffexp_regresson_genoalt/kinase_phosphotase_substrate_diffexp.txt")
us_diffexp$Cancer <- factor(us_diffexp$Cancer, levels = c("BRCA", "OV", "CO"))
us_diffexp$consistent <- FALSE
us_diffexp$consistent[(us_diffexp$upstream_direction == us_diffexp$downstream_direction & us_diffexp$upstream == "kinase")] <- TRUE
us_diffexp$consistent[(us_diffexp$upstream_direction != us_diffexp$downstream_direction & us_diffexp$upstream == "phosphotase")] <- TRUE
us_diffexp_nona <- us_diffexp[!is.na(us_diffexp$upstream_direction) & !is.na(us_diffexp$downstream_direction),]
us_diffexp_nona <- data.frame(us_diffexp_nona)

## input pathways
load("~/Box Sync/pan3can_shared_data/analysis_results/2015-08-01_Gene_Set.RData")
keywords <- c("adherens","Mismatch","EMT", "hormone", "apoptosis","immunological","stromal","transmembrane","receptors","integrin",
              "TGFβ","LKB1","AMPK","TSC","mTOR","PI3K","Akt","Ras","MAPK","Notch","Wnt", "catenin", "cell cycle", "p53","RTK","erbb")
all_keggs <- names(KEGG)
keggs <- NULL
for (key in keywords) {
  keggs <- c(keggs, all_keggs[grepl(key, all_keggs, ignore.case = T)])
}
keggs <- unique(keggs)

# plot MTOR cascade as a example ------------------------------------------
# for (hub in c("MTOR", "EGFR", "RB1")) {
# for (hub in c("MTOR")) {
for (hub in c("APC")) {

  ## expand proteins to upstream of the hub and limit to KEGG pathway components
  us_hub <- us_diffexp[(us_diffexp$KINASE == hub | us_diffexp$SUBSTRATE == hub) & us_diffexp$SELF == "trans",]
  all_genes <- unique(c(as.vector(us_hub$KINASE), as.vector(us_hub$SUBSTRATE)))
  for (kegg in keggs) {
    gene_kg <- all_genes[all_genes %in% KEGG[[kegg]]]
    if (length(gene_kg) > 1 & (hub %in% gene_kg)) {
      print(paste0(kegg, ":   ",paste(gene_kg, collapse = " ")))
      # df[, kegg] <- ifelse(df$SUBSTRATE %in% gene_kg, TRUE, FALSE)
      us_hub_kegg <- us_hub[us_hub$KINASE %in% KEGG[[kegg]] & us_hub$SUBSTRATE %in% KEGG[[kegg]],]
      
      # for (cancer in unique(us_hub$Cancer)) {
      for (cancer in c("CO")) {
        us_hub_can <- us_hub_kegg[us_hub_kegg$Cancer == cancer,]
        all_genes_kegg <- unique(c(as.vector(us_hub_can$KINASE), as.vector(us_hub_can$SUBSTRATE)))
        
        # hub %in% us_hub_can$KINASE
        if (TRUE) {
          ## input phospho data
          pho <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PHO", "_formatted_normalized_replicate_averaged.txt"),
                       data.table = F)
          # pho <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PHO", "_formatted.txt"),
          #              data.table = F)
          pho.n <- get_normal(expression = pho, clinical.m = clinical)
          pho.t <- get_tumor(expression = pho, clinical.m = clinical)
          pho_head <- formatPhosphosite(phosphosite_vector = pho$Phosphosite, pho$Gene)
          phog <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "collapsed_PHO", "_formatted_normalized_replicate_averaged.txt"),
                        data.table = F)
          phog <- phog[,intersect(c("Gene", colnames(pho.n), colnames(pho.t)), colnames(phog))]
          pho <- cbind(pho_head, pho[,intersect(c(colnames(pho.n), colnames(pho.t)), colnames(phog))])
          
          ## input somatic mutations
          vcf <- fread(input = paste0("/Users/yigewu/Box Sync/cptac2p/mutations/somatic/", cancer, ".maf"), data.table = F)
          vcf <- data.frame(vcf)
          vcf_hub <- vcf[vcf$Hugo_Symbol %in% all_genes_kegg,]
          rm(vcf)
          vcf_hub_nonsilent <- data.frame(table(vcf_hub[vcf_hub$Variant_Classification != "Silent", c("Hugo_Symbol", "Tumor_Sample_Barcode")]))
          vcf_hub_nonsilent <- vcf_hub_nonsilent[vcf_hub_nonsilent$Freq > 0,]
          colnames(vcf_hub_nonsilent)[3] <- "mut_Freq"
          
          ## input somatic copy number
          cnv <- fread(input = paste0("/Users/yigewu/Box Sync/cptac2p/copy_number/gatk/v1.1.prospective_only/deliverables/", str_sub(tolower(cancer), 1, 2), "/gene_level_CNV.", str_sub(tolower(cancer), 1, 2), ".v1.1.2017-12-10.tsv"), 
                       data.table = F, header = T)
          cnv_hub <- cnv[cnv$gene %in% all_genes_kegg,]
          rm(cnv)
          cnv_hub_m <- melt(cnv_hub)
          cnv_hub_m <- cnv_hub_m[cnv_hub_m$value < log2(1.5/2) | cnv_hub_m$value > log2(2.5/2),]
          colnames(cnv_hub_m) <- c("SUBSTRATE", "Tumor_Sample_Barcode", "log_copy")
          cnv_hub_m$Tumor_Sample_Barcode <- paste0(cnv_hub_m$Tumor_Sample_Barcode, "_T")
          
          ## merge kinase global phospho and substrate phosphosite phospho tgt
          pho_sub_rows <- sapply(1:nrow(us_hub_can), function(n, pho) which(pho$SUBSTRATE == us_hub_can$SUBSTRATE[n] & pho$SUB_MOD_RSD == us_hub_can$SUB_MOD_RSD[n]), pho = pho)
          pho_sub <- pho[pho_sub_rows,]
          pho_sub$transcript <- NULL
          pho_sub_m <- melt(pho_sub)
          pho_kin <- phog[phog$Gene %in% as.vector(us_hub_can$KINASE),]
          pho_kin_m <- melt(pho_kin)
          pho_kin_f <- data.frame(SUB_MOD_RSD = "global_phospho", SUBSTRATE = pho_kin_m$Gene, variable = pho_kin_m$variable, value = pho_kin_m$value)
          
          ## add order of phosphosites and annotate samples to tumors and normals
          # df <- merge(df, plot_order, by.x = c("SUBSTRATE"), by.y = c("Var1"), all.x = T)
          # df$Freq[is.na(df$Freq)] <- 0
          df <- rbind(pho_kin_f, pho_sub_m)
          df$y_print <- paste0(df$SUBSTRATE, ':', df$SUB_MOD_RSD)
          
          ## cap the value to show readable color
          lim = max(abs(max(df$value, na.rm = T)), abs(min(df$value, na.rm = T)))
          cap <- min(2, lim)
          df$pho_capped <- as.numeric(df$value)
          df$pho_capped[df$pho_capped > cap] <- cap
          df$pho_capped[df$pho_capped < (-cap)] <- (-cap)
          
          ## label samples to tumor/normal
          df$tumor_normal <- NA
          df$tumor_normal[df$variable %in% colnames(pho.t)] <- "tumor"
          df$tumor_normal[df$variable %in% colnames(pho.n)] <- "normal"
          
          ## add somatic mutations
          df$partID <- sampID2partID(sampleID_vector = as.vector(df$variable), sample_map = clinical)
          df$Tumor_Sample_Barcode <- NA
          df$Tumor_Sample_Barcode[df$tumor_normal == "tumor"] <- paste0(df$partID[df$tumor_normal == "tumor"], "_T")
          df <- merge(df, vcf_hub_nonsilent, by.x = c("SUBSTRATE", "Tumor_Sample_Barcode"), by.y = c("Hugo_Symbol", "Tumor_Sample_Barcode"), all.x = T, sort = F)
          
          ## add somatic CNV
          df <- merge(df, cnv_hub_m, all.x = T)
          
          ## add differential expression result
          diffexp_directions <- rbind(data.frame(y_print = paste0(us_hub_can$KINASE, ":global_phospho"), diffexp_type = paste0(ifelse(us_hub_can$KINASE == hub, "hub_", "upstream_"), us_hub_can$upstream_direction)),
                                      data.frame(y_print = paste0(us_hub_can$SUBSTRATE, ":", us_hub_can$SUB_MOD_RSD), diffexp_type = paste0(ifelse(us_hub_can$SUBSTRATE == hub, "hub_", "downstream_"), us_hub_can$downstream_direction)))
          diffexp_directions <- unique(diffexp_directions)
          df <- merge(df, diffexp_directions, all.x = T)
          df$diffexp_type <- factor(df$diffexp_type, levels = c("upstream_up", "upstream_down", "upstream_NA", "hub_up", "hub_down", "downstream_up", "downstream_down", "downstream_NA"))
          ## order samples by hub phospho
          pho_kin_hub <- pho_kin_m[pho_kin_m$Gene == hub,]
          df$variable <- factor(df$variable, levels = as.vector(pho_kin_hub$variable[order(pho_kin_hub$value)]))
          
          for (kegg_tmp in keggs) {
            gene_tmp <- all_genes_kegg[all_genes_kegg %in% KEGG[[kegg_tmp]]]
            if (length(gene_tmp) >0) {
              df[, kegg_tmp] <- ifelse(df$SUBSTRATE %in% gene_tmp, TRUE, FALSE)
            }
          }
          df$kegg_hit <- rowSums(df[,colnames(df) %in% keggs])
          
          ## keep hub upstream on the top
          df$kegg_hit[df$SUBSTRATE %in% as.vector(us_hub_can$KINASE[us_hub_can$SUBSTRATE == hub])] <- df$kegg_hit[df$SUBSTRATE %in% as.vector(us_hub_can$KINASE[us_hub_can$SUBSTRATE == hub])] + 100
          df$y_print <- reorder(df$y_print, df$kegg_hit)
          
          ## heatmap
          p = ggplot(df)
          # p = p + coord_fixed(ratio = 1)
          p = p + geom_tile(aes(x=variable, y=y_print, fill=pho_capped, color = log_copy), width = 0.7, height = 0.7, size=0.5)#, linetype="blank") 
          p = p + facet_grid(diffexp_type~tumor_normal, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
          p = p + scale_fill_gradientn(name= "phospho_level", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
          p = p + geom_text(aes(x=variable, y=y_print, label = mut_Freq), size = 1.5)
          p = p + theme_bw() + theme_nogrid()
          p = p + theme(axis.text.x = element_blank())
          p = p + theme(axis.title = element_blank(),  axis.text.y = element_text(size = 7))
          p = p + theme(axis.ticks=element_blank(),legend.position="bottom")
          p = p + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
          p = p + theme(strip.text.y = element_text(size = 5, angle = 0))
          p = p + scale_colour_gradient2(high = "#E31A1C", low = "#1F78B4", midpoint = 0, na.value = NA)
          p
          resultDnow <- makeOutDir()
          ggsave(filename = paste0(resultDnow, hub, "_substrates_and_upstream_in_", kegg, "_in_", cancer, ".pdf"), 
                 height = 3, width = 12)
        }
      }
    }
  }
}


## plot adjancency matrix
# react_anno <- data.frame(gene = all_genes_kegg)
# for (react in names(REACT)) {
#   gene_rt <- all_genes_kegg[all_genes_kegg %in% REACT[[react]]]
#   if (length(gene_rt) >2) {
#     print(paste0(react, ":   ",paste(gene_rt, collapse = " ")))
#     react_anno[, react] <- ifelse(react_anno$gene %in% gene_rt, TRUE, FALSE)
#   }
# }
# react_anno <- react_anno[rowSums(react_anno == TRUE) > 0,]
# react_anno_m <- melt(react_anno)
# ggplot(react_anno_m, aes(x = variable, y = gene, fill = value)) +
#   geom_raster() +
#   theme_bw() +
#   # Because we need the x and y axis to display every node,
#   # not just the nodes that have connections to each other,
#   # make sure that ggplot does not drop unused factor levels
#   scale_x_discrete(drop = FALSE, position = "top") +
#   scale_y_discrete(drop = FALSE) +
#   theme(
#     # Rotate the x-axis lables so they are legible
#     axis.text.x = element_text(angle = 270, hjust = 0),
#     # Force the plot into a square aspect ratio
#     aspect.ratio = 1,
#     # Hide the legend (optional)
#     legend.position = "none")
# 
# 
