# Yige Wu @ WashU 2018 Aug
## normalize the log2 copy number of CPTAC2 prospective samples

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
library(ggplot2)


# set variables -----------------------------------------------------------

# inputs ------------------------------------------------------------------

# loop each cancer --------------------------------------------------------
cnv <- fread(input = "./cptac2p/copy_number/gatk4wxscnv/v1.3.swap_contamination_fixed/prospective_somatic/deliverables/ov/gene_level_CNV.ov.v1.3.2018-03-19.tsv", data.table = F)
cnv <- fread(input = "./cptac2p/copy_number/gatk4wxscnv/v1.3.swap_contamination_fixed/prospective_somatic/deliverables/br/gene_level_CNV.br.v1.3.2018-03-19.tsv", data.table = F)
int_enzyme_cnv <- intersect(unique(cnv$gene), enzymes)


for (cancer in c("OV")) {
  ## input gene-level CNV
  cnv <- fread(input = "./cptac2p/copy_number/gatk4wxscnv/v1.3.swap_contamination_fixed/prospective_somatic/deliverables/ov/gene_level_CNV.ov.v1.3.2018-03-19.tsv", data.table = F)
  for (gene in cnv$gene) {
    cnv_gene <- cnv[cnv$gene == gene, -1]
    tab2p <- data.frame(log2cn = as.numeric(t(cnv_gene)))
    p <- ggplot()
    p <-p + geom_density(data = tab2p, mapping = aes(x = log2cn))
    p
  }
}
