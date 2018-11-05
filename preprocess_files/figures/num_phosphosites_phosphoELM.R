# Yige Wu @ WashU 2018 July
## show the number of phosphosites detected in prospective collection overlapping with phosphoELM 

# source ------------------------------------------------------------------
source("Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R")

library(readr)
library(eulerr)
library(ggrepel)

# set variables -----------------------------------------------------------
cancers_sort <- c("BRCA", "OV", "CO")
num_nonNA2test <- c(15)
num_sites2test <- c(1)

# input CPTAC prospective phosphosite genomic coordinates(hg19) -----------------------------------
pho_head_geno_pos <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid_hg38_hg19.txt"), data.table = F)
pho_head_geno_pos$id <- paste0(pho_head_geno_pos$SUBSTRATE, "_", pho_head_geno_pos$SUB_MOD_RSD)
pho_head_geno_pos <- pho_head_geno_pos[!duplicated(pho_head_geno_pos$id) & !is.na(pho_head_geno_pos$coordinate_hg19),]
rownames(pho_head_geno_pos) <- pho_head_geno_pos$id
## create the column that contains the formatted genomic coordinates
pho_head_geno_pos$genomic_pos <- paste0(pho_head_geno_pos$chromosome, ":g.", pho_head_geno_pos$start_hg19, "_", pho_head_geno_pos$end_hg19)

# input phosphoELM --------------------------------------------------------
phosphoELM <- read_delim("pan3can_shared_data/Phospho_databases/Phospho.ELM/phosphoELM_all_2015-04.dump_human_transvar_out.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
## clean phoshpoELM
### filter out sites cannot be tranformed to genomic coordinates and duplicate coordinates
phosphoELM <- data.frame(phosphoELM)
phosphoELM <- phosphoELM[!duplicated(phosphoELM$coordinates.gDNA.cDNA.protein.),]
phosphoELM <- phosphoELM[phosphoELM$coordinates.gDNA.cDNA.protein. != "././.",]
### create the column that contains only the genomic coordinates
phosphoELM$chr <- str_split_fixed(string = phosphoELM$coordinates.gDNA.cDNA.protein., pattern = "\\:", 2)[,1]
phosphoELM$start <- str_split_fixed(string = phosphoELM$coordinates.gDNA.cDNA.protein., pattern = "\\:g\\.|\\_|\\/c\\.", 4)[,2]
phosphoELM$end <- str_split_fixed(string = phosphoELM$coordinates.gDNA.cDNA.protein., pattern = "\\:g\\.|\\_|\\/c\\.", 4)[,3]
phosphoELM$genomic_pos <- paste0(phosphoELM$chr, ":g.", phosphoELM$start, "_", phosphoELM$end)
## create a column represent the site
phosphoELM$SUB_MOD_RSD <- str_split_fixed(string = phosphoELM$input, pattern = "\\:", 2)[,2]
phosphoELM$site <- paste0(phosphoELM$gene, "_", phosphoELM$SUB_MOD_RSD)
phosphoELM$phosphoELM <- T
# input disease-associated phosphosites from PSP --------------------------
## input file
psp_disease <- read_delim("Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Disease-associated_sites_term_cancer.oma.leuk_wGPos.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
psp_disease <- data.frame(psp_disease)
## clean up
psp_disease <- psp_disease[psp_disease$genomic_pos != "." & !duplicated(psp_disease$genomic_pos),]
## create the column represent the phosphosite
psp_disease$SUB_MOD_RSD <- str_split_fixed(string = psp_disease$MOD_RSD, "\\-", 2)[,1]
psp_disease$site <- paste0(psp_disease$GENE, "_", psp_disease$SUB_MOD_RSD)
psp_disease$psp_disease <- T

# input phosphosites detected in prospective collection and plot venn diagram -------------------
for (num_nonNA in num_nonNA2test) {
  for (num_sites in num_sites2test) {
    ## get all the phosphosites from prospective collection (only single site phosphosite)
    
    pho_sites_cans <- list()
    pho_sites_all <- NULL
    for (cancer in cancers_sort) {
      pho_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
      samples <- colnames(pho_data)[!(colnames(pho_data) %in% c("Gene", "Phosphosite"))]
      pho_data <- pho_data[rowSums(!is.na(pho_data[, samples])) >= num_nonNA,]
      pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
      ## only take single phosphosite
      pho_head$offset2 <- str_split_fixed(string = pho_head$SUB_MOD_RSD, pattern = '[STY]', 3)[,3]
      pho_head <- pho_head[pho_head$offset2 == "",]
      pho_sites <- unique(pho_head[, c("SUBSTRATE", "SUB_MOD_RSD")])
      pho_sites$site <- paste0(pho_sites$SUBSTRATE, "_", pho_sites$SUB_MOD_RSD)
      pho_sites_cans[[cancer]] <- unique(as.vector(pho_sites$site))
      pho_sites_all <- unique(c(pho_sites_all, as.vector(pho_sites$site)))
    }
    
    ## merge genomic positions into prospective phosphosites
    dat <- data.frame(site = pho_sites_all)
    dat$CPTAC <- TRUE
    dat <- merge(dat, pho_head_geno_pos[, c("id", "genomic_pos")], by.x = c("site"), by.y = c("id"), all.x = T)
    ## take out those without genomic coordinates first
    dat <- dat[!is.na(dat$genomic_pos),]
    
    ## merge phosphosites from phosphoELM
    dat <- merge(dat, phosphoELM[, c("genomic_pos", "phosphoELM", "site")], by = c("genomic_pos"), all = T, suffix = c(".cptac", ".elm"))
    ## merge phosphosites from PSP
    # dat <- merge(dat, psp_disease[, c("genomic_pos", "psp_disease", "site")], by = c("genomic_pos"), all = T, suffix = c(".cptac", ".psp"))
    # dat$psp_disease[is.na(dat$psp_disease)] <- FALSE
    
    ## clean up the logical colunns
    dat$CPTAC[is.na(dat$CPTAC)] <- FALSE
    dat$phosphoELM[is.na(dat$phosphoELM)] <- FALSE
    
    ## make venn diagram for phosphosites
    fn = paste(makeOutDir(resultD = resultD), 'all_phosphosites_', num_nonNA, 'nonNA_venn','.pdf',sep ="")
    grid.newpage()
    pdf(fn, height = 12, width = 12, useDingbats = FALSE)
    fit <- euler(combinations = dat, input = "disjoint", shape = 'circle')
    
    p <-plot(fit, quantities = list(fontsize = 50), labels = F, legend = list(fontsize = 50))
    grid.draw(p)
    dev.off()
  }
}
# check how phosphatase are distributed in detected phosphosites across historical classification -------
genes_pho <- unique(str_split_fixed(string = pho_sites_all, pattern = "_", 2)[,1])
genes_pho <- data.frame(gene = genes_pho)

## input human phosphatases from DEPOD
depod_tab <- read_excel("Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/pan3can_shared_data/Phospho_databases/DEPOD/DEPOD_201410_human_phosphatases.xlsx")

## merge by gene symbols
genes_pho <- merge(genes_pho, depod_tab, by.x = c("gene"), by.y = c("Gene symbol"))

## calculate the number of historical classification
sort(table(genes_pho$`Historical class`))

# check how kinases are distributed aross detected phosphosites -------
genes_pho <- unique(str_split_fixed(string = pho_sites_all, pattern = "_", 2)[,1])

## check how many auto-phosphorylating kinases are detected
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source != "NetKIN" | (ptms_site_pairs_sup$Source == "NetKIN" & ptms_site_pairs_sup$networkin_score >= 5),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)
ptms_site_pairs_sup$pair <- paste0(ptms_site_pairs_sup$pair_pro, ":", ptms_site_pairs_sup$SUB_MOD_RSD)
enzyme_type <- "kinase"
k_s_table <- ptms_site_pairs_sup[ptms_site_pairs_sup$enzyme_type == enzyme_type,]
k_s_psp <- load_ks_table(protein = enzyme_type)
kinase_trans <- as.vector(unique(k_s_table$GENE[as.vector(k_s_table$GENE)!=as.vector(k_s_table$SUB_GENE)]))
kinase_cis <- unique(c(as.vector(k_s_table$GENE[as.vector(k_s_table$GENE) == as.vector(k_s_table$SUB_GENE)]), 
                       as.vector(k_s_psp$GENE[as.vector(k_s_psp$GENE)==as.vector(k_s_psp$SUB_GENE)])))

length(unique(intersect(genes_pho, kinase_cis)))

reg_nonNA2test <- c(20)
fdr_thres <- 0.05
cptac_phase2process <- "cptac2p_3can"
for (reg_nonNA in reg_nonNA2test) {
  for (enzyme_type in c("kinase")) {
    regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                       enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor",
                                       "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
    regression <- markSigKS(regression = regression, sig_thres = fdr_thres, enzyme_type = enzyme_type)
    regression$regulated <- (regression$fdr_sig & regression$coef_sig)
    
    cis_enzymes <- unique(regression$GENE[regression$regulated & regression$SELF == "cis"])
    
  }
}
regression$site <- paste0(regression$SUBSTRATE, "_", regression$SUB_MOD_RSD)

length(cis_enzymes)
length(cis_enzymes)/length(unique(intersect(genes_pho, kinase_cis)))

## get the phosphosites known to be autophosphorylated
# k_s_cis <- k_s_table[as.vector(k_s_table$GENE) == as.vector(k_s_table$SUB_GENE),]
k_s_cis <- k_s_psp[as.vector(k_s_psp$GENE) == as.vector(k_s_psp$SUB_GENE),]
k_s_cis$site <- paste0(k_s_cis$GENE, "_", k_s_cis$SUB_MOD_RSD)
k_s_cis$pair <- paste0(k_s_cis$GENE, ":", k_s_cis$SUB_GENE, ":", k_s_cis$SUB_MOD_RSD)

length(intersect(k_s_cis$site, pho_sites_all))
known_cis_sites <- intersect(k_s_cis$site, pho_sites_all)
known_cis_sites_genes <- unique(str_split_fixed(string = known_cis_sites, pattern = "_", 2)[,1]) 
length(known_cis_sites_genes)

wilcox.test(regression$coef_pro_kin[regression$pair %in% k_s_cis$pair], regression$coef_pro_kin[!(regression$pair %in% k_s_cis$pair) & (regression$GENE %in% known_cis_sites_genes) & regression$SELF == "cis"])

for (cancer in cancers_sort) {
  print(wilcox.test(regression$coef_pro_kin[regression$pair %in% k_s_cis$pair & regression$Cancer == cancer], regression$coef_pro_kin[!(regression$pair %in% k_s_cis$pair) & (regression$GENE %in% known_cis_sites_genes) & regression$SELF == "cis" & regression$Cancer == cancer]))
}

## getting known cis-regulated sites also cis-regulated in our datasets
known_cis_sites_regulated <- regression[regression$pair %in% k_s_cis$pair & regression$regulated,]
write.table(x = known_cis_sites_regulated, file = paste0(makeOutDir(resultD = resultD), "known_cis_sites_regulated.txt"), quote = F, sep = "\t", row.names = F)
length(unique(known_cis_sites_regulated$pair))

k_s_cis_regulated <- k_s_cis[k_s_cis$pair %in% regression$pair[regression$pair %in% k_s_cis$pair & regression$regulated],]


# examing cis-regulated kinases -------------------------------------------
## get common cis-regulated kinases
BRCA_cis_kinases <- unique(regression$GENE[regression$regulated & regression$SELF == "cis" & regression$Cancer == "BRCA"])
OV_cis_kinases <- unique(regression$GENE[regression$regulated & regression$SELF == "cis" & regression$Cancer == "OV"])
CO_cis_kinases <- unique(regression$GENE[regression$regulated & regression$SELF == "cis" & regression$Cancer == "CO"])

shared_cis_kinases <- intersect(CO_cis_kinases, intersect(BRCA_cis_kinases, OV_cis_kinases))
sort(unlist(map2TCGApathwaways(gene_list = shared_cis_kinases, pathway_list = tcga_pathways)))
sort(unlist(map2TCGApathwaways(gene_list = shared_cis_kinases, pathway_list = tcga_pathways_pluskegg)))
table(sort(unlist(map2TCGApathwaways(gene_list = shared_cis_kinases, pathway_list = tcga_pathways_pluskegg))))
unlist(map2TCGApathwaways(gene_list = shared_cis_kinases, pathway_list = tcga_pathways_pluskegg_and_pathway))

## get common cis-regulated phosphosites
BRCA_cis_sites <- unique(regression$site[regression$regulated & regression$SELF == "cis" & regression$Cancer == "BRCA"])
OV_cis_sites <- unique(regression$site[regression$regulated & regression$SELF == "cis" & regression$Cancer == "OV"])
CO_cis_sites <- unique(regression$site[regression$regulated & regression$SELF == "cis" & regression$Cancer == "CO"])

shared_cis_sites <- intersect(CO_cis_sites, intersect(BRCA_cis_sites, OV_cis_sites))

## get cancer type specific kinases

## get cancer type specific cis-regulated phosphosites


