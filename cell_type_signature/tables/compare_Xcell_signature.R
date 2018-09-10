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
  
  for (cell_type in unique(xcell_raw$cell_type)) {
    all_genes <- xcell_list[[cell_type]]
    all_genes <- all_genes[all_genes %in% pro_t$Gene & all_genes %in% pro_n$Gene]
    if (length(all_genes) > 0) {
      xcell_wilcox[[cancer]][[cell_type]] <- list()
      xcell_wilcox[[cancer]][[cell_type]][["pvalues_ordered"]] <- list()
      xcell_wilcox[[cancer]][[cell_type]][["pvalues_fisher"]] <- list()
      
      for (side in c("greater", "less")) {
        xcell_wilcox[[cancer]][[cell_type]][[side]] <- list()
        for (gene in all_genes) {
          pro_t_g <- as.numeric(as.vector(pro_t[pro_t$Gene == gene, samples_t]))
          pro_n_g <- as.numeric(as.vector(pro_n[pro_n$Gene == gene, samples_n]))
          
          stat <- try(wilcox.test(x = pro_t_g, y = pro_n_g, alternative = side))
          if ("p.value" %in% names(stat)) {
            xcell_wilcox[[cancer]][[cell_type]][[side]][[gene]] <- stat$p.value
          } else {
            xcell_wilcox[[cancer]][[cell_type]][[side]][[gene]] <- NA
          }
        }
        ## order genes by p-values
        pvalues <- sapply(all_genes, FUN = function(g, l) l[[g]], l = xcell_wilcox[[cancer]][[cell_type]][[side]])
        names(pvalues) <- all_genes
        pvalues <- pvalues[!is.na(pvalues)]
        pvalues <- pvalues[order(pvalues)]
        xcell_wilcox[[cancer]][[cell_type]][["pvalues_ordered"]][[side]] <- pvalues
        
        # ## choose top 10 significant genes and combine the pvalues
        # pvalues_sig <- pvalues[pvalues < sig]
        # pvalues_top <- pvalues_sig[1:min(10, length(pvalues_sig))]
        # if (length(pvalues_top) > 1) {
        #   p.fisher <- sumlog(pvalues_top)
        #   p.fisher <- p.fisher$p
        #   xcell_wilcox[[cancer]][[cell_type]][["pvalues_fisher"]][[side]] <- p.fisher
        # } else {
        #   xcell_wilcox[[cancer]][[cell_type]][["pvalues_fisher"]][[side]] <- pvalues_top
        # }
        
        ## combine the all the pvalues
        p.fisher <- sumlog(pvalues)
        p.fisher <- p.fisher$p
        xcell_wilcox[[cancer]][[cell_type]][["pvalues_fisher"]][[side]] <- p.fisher
      }
      xcell_wilcox[[cancer]][[cell_type]][["pvalues_fisher"]][["min"]] <- min(xcell_wilcox[[cancer]][[cell_type]][["pvalues_fisher"]][["greater"]],
                                                                    xcell_wilcox[[cancer]][[cell_type]][["pvalues_fisher"]][["less"]])
    }
  }
  pvalue_greater <- sapply(unique(xcell_raw$cell_type), FUN = function(c, l) {
    tmp <- NA
    if ( c %in% names(l) ) {
      tmp <- l[[c]][["pvalues_fisher"]][["greater"]]
    }
    return(tmp)
    }, l = xcell_wilcox[[cancer]])
  
  pvalue_less <- sapply(unique(xcell_raw$cell_type), FUN = function(c, l) {
    tmp <- NA
    if ( c %in% names(l) ) {
      tmp <- l[[c]][["pvalues_fisher"]][["less"]]
    }
    return(tmp)
  }, l = xcell_wilcox[[cancer]])
  
  tmp1<- data.frame(cancer = cancer, side = "tumor>normal"); tmp1 <- cbind(tmp1, t(pvalue_greater))
  tmp2<- data.frame(cancer = cancer, side = "tumor<normal"); tmp2 <- cbind(tmp2, t(pvalue_less))
  tmp <- rbind(tmp1, tmp2)
  pvalue_fisher <- rbind(pvalue_fisher, tmp)
  # pvalue_fisher <- sapply(unique(xcell_raw$cell_type), FUN = function(c, l) l[[c]][["pvalue_fisher"]][["min"]], l = xcell_wilcox[[cancer]])
}

dir.create("./cptac2p/cptac_shared/analysis_results/cell_type_signature/")
dir.create("./cptac2p/cptac_shared/analysis_results/cell_type_signature/tables")
saveRDS(object = xcell_wilcox, file = paste0(makeOutDir(), "xcell_wilcox.RDS"))
write.table(x = pvalue_fisher, file = paste0(makeOutDir(), "pvalue_fisher.txt"), 
            row.names = F, sep = '\t')

CD8+ Tem
