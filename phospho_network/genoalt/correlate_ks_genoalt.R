# Yige Wu @ WashU 2018 Jan
# correlate substrate expression with kinase mutations/cnv

source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
resultDnow <- makeOutDir()


# inputs ------------------------------------------------------------------
cnv <- vector("list", length = 3)
maf <- vector("list", length = 3)
for (cancer in c("BRCA", "OV", "CO")) {
  cnv <- fread(input = "/Users/yigewu/Box Sync/cptac2p/copy_number/gatk/v1.1.prospective_only/deliverables/co/gene_level_CNV.co.v1.1.2017-12-10.tsv",
                  data.table = F)
  maf_co <- fread(input = "/Users/yigewu/Box Sync/cptac2p/mutations/somatic/")
}
