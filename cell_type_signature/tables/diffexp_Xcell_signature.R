# Yige Wu @ WashU 2018 Apr
# take cell type specific mark genes from Xcell and test difference between prospective tumors and normals


# souce ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/cptac2p_analysis_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(readr)
library(metap)
# inputs ------------------------------------------------------------------
## input Xcell cell type signature
xcell_raw <- read_delim("~/Box Sync/cptac2p/deconvolution/xCell_gene_signature_non_redundant.txt",
                        "\t", escape_double = FALSE, trim_ws = TRUE)
xcell_list <- lapply(1:nrow(xcell_raw), FUN = function(n, xcell_raw) str_split(string = as.vector(xcell_raw$signature_genes[xcell_raw$cell_type == xcell_raw$cell_type[n]]), pattern = ",")[[1]],
                     xcell_raw = xcell_raw)
names(xcell_list) <- as.vector(xcell_raw$cell_type)
##

cancers <- c("BRCA", "OV", "CO")
xcell_wilcox <- vector("list", length = length(cancers))
pvalue_fisher <- NULL
# loop around each cancer types -------------------------------------------
for (cancer in cancers) {
  xcell_wilcox[[cancer]] <- list()
  ## input protein expression
  pro_t <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PRO", "_formatted_normalized_noControl.txt"),
               data.table = F)
  samples_t <- colnames(pro_t)[!(colnames(pro_t) %in% c("Gene"))]
  pro_n <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PRO", "_formatted_normalized_Control.txt"),
                 data.table = F)
  samples_n <- colnames(pro_n)[!(colnames(pro_n) %in% c("Gene"))]
  
  ## input differential protein expression result
  diff_pro <- fread(paste0(resultD, "diffexp/tables/samr_paired/", cancer, "_", datatype, "_diffexp_paired_", testType, ".txt"),
                data.table = F)
  diff_pro$log2FC <- log2(diff_pro$foldchange)
  
  for (cell_type in unique(xcell_raw$cell_type)) {
    all_genes <- xcell_list[[cell_type]]
    all_genes <- all_genes[all_genes %in% pro_t$Gene & all_genes %in% pro_n$Gene]
    xcell_wilcox[[cancer]][[cell_type]] <- list()
    xcell_wilcox[[cancer]][[cell_type]][["genes"]] <- all_genes
    xcell_wilcox[[cancer]][[cell_type]][["log2fc_ordered"]] <- list()
    xcell_wilcox[[cancer]][[cell_type]][["fdrs_ordered"]] <- list()
    xcell_wilcox[[cancer]][[cell_type]][["fdrs_fisher"]] <- list()
    
    ## order fdrs
    fdrs <- diff_pro$fdr[diff_pro$Gene %in% all_genes]
    names(fdrs) <- diff_pro$Gene[diff_pro$Gene %in% all_genes]

    ## order fold change by fdrs
    log2fcs <- diff_pro$log2FC[diff_pro$Gene %in% all_genes]
    names(log2fcs) <- diff_pro$Gene[diff_pro$Gene %in% all_genes]
    log2fcs <- log2fcs[order(fdrs)]
    
    fdrs <- fdrs[order(fdrs)]
    xcell_wilcox[[cancer]][[cell_type]][["fdrs_ordered"]] <- fdrs
    xcell_wilcox[[cancer]][[cell_type]][["log2fc_ordered"]] <- log2fcs
    
    fdrs_top <- fdrs[1:min(10, length(fdrs))]
    if (length(fdrs_top) > 1) {
      p.fisher <- sumlog(fdrs_top)
      p.fisher <- p.fisher$p
      xcell_wilcox[[cancer]][[cell_type]][["fdrs_fisher"]] <- p.fisher
    } else {
      xcell_wilcox[[cancer]][[cell_type]][["fdrs_fisher"]] <- fdrs_top
    }
  }
  tmp <- sapply(unique(xcell_raw$cell_type), FUN = function(c, l) l[[c]][["fdrs_fisher"]], l = xcell_wilcox[[cancer]])
  tmp <- data.frame(tmp)
  colnames(tmp) <- cancer
  if (is.null(pvalue_fisher)) {
    pvalue_fisher <- tmp
  } else {
    pvalue_fisher <- cbind(pvalue_fisher, tmp)
  }
}

dir.create("./cptac2p/cptac_shared/analysis_results/cell_type_signature/")
dir.create("./cptac2p/cptac_shared/analysis_results/cell_type_signature/tables")
saveRDS(object = xcell_wilcox, file = paste0(makeOutDir(), "xcell_wilcox.RDS"))
pvalue_fisher$cell_type <- rownames(pvalue_fisher)
write.table(x = pvalue_fisher, file = paste0(makeOutDir(), "pvalue_fisher.txt"), 
            row.names = F, sep = '\t')

