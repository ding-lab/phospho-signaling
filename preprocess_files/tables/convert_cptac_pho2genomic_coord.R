# Yige Wu @ WashU 2018 Jan
# convert CPTAC global proteomids phosphosites to genomic positions
# refseq IDs to ensembl peptide ids to genomic coordinates

# source -------------------------------------------------------------------
## for convert uniprot ID to ensembl peptide ID
source("https://bioconductor.org/biocLite.R")
# biocLite("devtools")
# BiocInstaller::biocLite('grimbough/biomaRt')
library(biomaRt)

## for convert ensembl peptide ID + residue to genomic coord
library(httr)
library(jsonlite)
library(xml2)

## others
library(data.table)
library(readr)
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


# source ensembl -----------------------------------------------------------
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org") ## 38
filters = listFilters(ensembl)
ensembl37 = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="grch37.ensembl.org")

# input bobo’s PTMcosmos mapping files --------------------------------
uniprot_ensembl_map <- fread("./cptac2p/resources/PTMcosmos/Data dump/2019-01-11/all_uniprot_ensembl_id_mappings.tsv", data.table = F)
ptm_sites_map <- fread("./cptac2p/resources/PTMcosmos/Data dump/2019-01-11/all_ptm_sites.all_txs.tsv", data.table = F)


# input already mapped phosphosites ---------------------------------------
pho_head_mapped_previous <- fread(input = "./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/20180724/pho_head_cancers_single_site_wpid_hg38_hg19.txt", data.table = F)
pho_head_mapped_previous %>%
  head

pho_head_mapped_previous <- pho_head_mapped_previous %>%
  mutate(id = paste0(SUBSTRATE, "_", SUB_MOD_RSD))

# set cancers to process --------------------------------------------------
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")

# gather all phosphosites -------------------------------------------------
pho_head_cancers <- NULL
for (cancer in cancers2process) {
  if (cancer %in% c("BRCA", "OV", "CO")) {
    pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
  } else if (cancer == "UCEC") {
    pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
  } else if (cancer == "CCRCC") {
    pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
  } else if (cancer == "LIHC") {
    pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
  }
  
  pho_head <- data.frame(Gene = pho_tab$Gene, Phosphosite = pho_tab$Phosphosite)
  pho_head_cancers <- rbind(pho_head_cancers, pho_head)
}
pho_head_cancers <- unique(pho_head_cancers)

# mark whether it's solo phosphosite or multiple phosphosite --------------
tmp <- str_split_fixed(string = pho_head_cancers$Phosphosite, pattern = "[STY]", 3)
pho_head_cancers$is.solo <- (tmp[,3] == "")
pho_head_cancers$p_coord <- tmp[,2]
pho_head_cancers$residue <- sapply(X = pho_head_cancers$Phosphosite, FUN = function(phosphosite) substr(x = phosphosite, start = 1, stop = 1))


# mark whether it has genomic coordination or not -------------------------
pho_head_cancers <- pho_head_cancers %>%
  mutate(id = paste0(Gene, "_", Phosphosite)) %>%
  mutate(is_previously_mapped = (id %in% pho_head_mapped_previous$id))

table(pho_head_cancers$is_previously_mapped)

pho_head_cancers %>%
  filter(is.solo == T) %>%
  select(id) %>%
  unique() %>%
  nrow()

write.table(x = pho_head_cancers, file = paste0(makeOutDir(resultD = resultD), "pho_head_cancers.txt"), quote = F, row.names = F, sep = "\t")

# fill in the ensembl protein ids using the already filled in ids ---------
pho_head_cancers2map <- pho_head_cancers %>%
  filter(is_previously_mapped == F) %>%
  filter(is.solo == T)

pho_head_cancers2map <- merge(pho_head_cancers2map, ptm_sites_map[, c("symbol", "p_coord", "residue", "ensembl_protein_id")], by.x = c("Gene", "p_coord", "residue"), by.y = c("symbol", "p_coord", "residue"), all.x = T)
pho_head_cancers2map_wpid <- pho_head_cancers2map %>%
  filter(!is.na(ensembl_protein_id) & ensembl_protein_id != "")

pho_head_cancers2map_2mappid <- pho_head_cancers2map %>%
  filter(is.na(ensembl_protein_id))
pho_head_cancers2map_2mappid$ensembl_protein_id <- NULL
pho_head_cancers2map_2mappid <- merge(pho_head_cancers2map_2mappid, unique(pho_head_cancers2map_wpid[, c("Gene", "ensembl_protein_id")]), all.x = T)
pho_head_cancers2map_2mappid_mapped <- pho_head_cancers2map_2mappid %>%
  filter(!is.na(ensembl_protein_id))
pho_head_cancers2map_wpid <- rbind(pho_head_cancers2map_wpid, pho_head_cancers2map_2mappid_mapped)

pho_head_cancers2map_2mappid <- pho_head_cancers2map_2mappid %>%
  filter(is.na(ensembl_protein_id))

pho_head_cancers2map_2mappid$ensembl_protein_id <- NULL
pho_head_cancers2map_2mappid$ensembl_protein_id <- NULL
pho_head_cancers2map_2mappid <- merge(pho_head_cancers2map_2mappid, unique(ptm_sites_map[ptm_sites_map$symbol %in% pho_head_cancers2map_2mappid$Gene, c("symbol", "ensembl_protein_id")]), by.x = c("Gene"), by.y = c("symbol"), all.x = T)
table(is.na(pho_head_cancers2map_2mappid$ensembl_protein_id))
pho_head_cancers2map_2mappid_mapped <- pho_head_cancers2map_2mappid %>%
  filter(!is.na(ensembl_protein_id))
pho_head_cancers2map_wpid <- rbind(pho_head_cancers2map_wpid, pho_head_cancers2map_2mappid_mapped)
pho_head_cancers2map_wpid <- unique(pho_head_cancers2map_wpid)

pho_head_cancers_sim <- pho_head_cancers2map_wpid

# query genomic coordinates based on ensembl id ---------------------------
coords <- list()
pids = as.vector(pho_head_cancers_sim$ensembl_protein_id)
p_coords = as.vector(pho_head_cancers_sim$p_coord)
for (n in (length(coords) + 1):nrow(pho_head_cancers_sim)) {
  if (!is.na(pids[n])) {
    result <- tryCatch({
      r <- GET(paste("http://rest.ensembl.org", paste0("/map/translation/", pids[n], "/", p_coords[n], "..", p_coords[n], "?"), sep = ""), content_type("application/json"))
      stop_for_status(r)
      tmp <- fromJSON(toJSON(content(r)))$mappings
      print(n)
      coords[[n]] <- tmp
    }, error = function(err) {
      print(paste("MY_ERROR:  ",err))
    })
  }
}

# saveRDS(object = coords, file = paste0(makeOutDir(resultD = resultD), "coords.RDS"))
df <- matrix(data = NA, nrow = nrow(pho_head_cancers_sim), ncol = ncol(coords[[1]]))
coords_names <- colnames(coords[[1]])
for (i in 1:nrow(pho_head_cancers_sim)) {
  if (!is.null(coords[[i]])) {
    if (nrow(coords[[i]]) == 1) {
      tmp <- unlist(coords[[i]])[coords_names]
      df[i,] <- matrix(tmp, nrow = 1, byrow = T)
    }
  }
}
df <- data.frame(df,stringsAsFactors=FALSE)
colnames(df) <- coords_names
pho_head_cancers_sim <- cbind(pho_head_cancers_sim, df)
pho_head_cancers_sim$is_previously_mapped <- NULL
pho_head_cancers_sim$is.solo <- NULL
pho_head_cancers_sim %>%
  head()
pho_head_mapped_previous2merge <- pho_head_mapped_previous %>%
  mutate(Gene = SUBSTRATE) %>%
  mutate(p_coord = offset) %>%
  mutate(Phosphosite = SUB_MOD_RSD) %>%
  mutate(ensembl_protein_id = ensembl_peptide_id) %>%
  mutate(residue = substr(SUB_MOD_RSD, start = 1, stop = 1)) %>%
  select(colnames(pho_head_cancers_sim)) %>%
  unique()

pho_head_cancers_sim <- rbind(pho_head_cancers_sim, pho_head_mapped_previous2merge)
write.table(x = pho_head_cancers_sim,
             file = paste0(makeOutDir(resultD = resultD), "pho_head_cancers_single_site_wpid_hg38.txt"),
             row.names = F, col.names = T, quote = F, sep = ",")
rm(coords)

pho_head_cancers_sim %>%
  select(id) %>%
  unique() %>%
  nrow()

# # lift over to hg19 -------------------------------------------------------
# pho_head_cancers_sim$chromosome <- paste0("chr", pho_head_cancers_sim$seq_region_name)
# pho_head_cancers_sim$coordinate_hg38 <- paste0(pho_head_cancers_sim$chromosome, ":", pho_head_cancers_sim$start, "-", pho_head_cancers_sim$end)
# test <- pho_head_cancers_sim[duplicated(pho_head_cancers_sim$coordinate_hg38) & !is.na(pho_head_cancers_sim$seq_region_name) ,]
# ## only keep primary sequences and no duplicated coordinates
# bedfile <- pho_head_cancers_sim[!is.na(pho_head_cancers_sim$seq_region_name) & !duplicated(pho_head_cancers_sim$coordinate_hg38) & !grepl(pattern = "HSCHR", x = pho_head_cancers_sim$seq_region_name), c("chromosome", "start", "end")]
# write.table(x = bedfile,
#             file = paste0(makeOutDir(resultD = resultD), "pho_head_cancers_single_site_wpid_hg38_coordinates2liftover.bed"),
#             row.names = F, col.names = F, quote = F, sep = "\t")
# bedfile_hg19 <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid_hg19.bed"), 
#                       data.table = F, col.names = c("chromosome", "start", "end"))
# bedfile_hg19$coordinate_hg19 <- paste0(bedfile_hg19$chromosome, ":", bedfile_hg19$start, "-", bedfile_hg19$end)
# ## read.delim left out a line
# bedfile_failed <- read_delim(file = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid_hg38_coordinates2liftover.err"), 
#                              comment = c("#"), col_names = c("chromosome", "start", "end"), delim = "\t")
# bedfile_failed$liftover_failed <- TRUE
# bedfile_failed$coordinate_hg38 <- paste0(bedfile_failed$chromosome, ":", bedfile_failed$start, "-", bedfile_failed$end)
# bedfile$coordinate_hg38 <- paste0(bedfile$chromosome, ":", bedfile$start, "-", bedfile$end)
# bedfile$liftover_failed <- (bedfile$coordinate_hg38 %in% bedfile_failed$coordinate_hg38)
# bedfile_hg38_liftovered <- bedfile[!bedfile$liftover_failed,]
# nrow(bedfile_hg38_liftovered) == nrow(bedfile_hg19)
# colnames(bedfile_hg19) <- c("chromosome_hg19", "start_hg19", "end_hg19", "coordinate_hg19")
# bedfile_hg38_liftovered <- cbind(bedfile_hg38_liftovered, bedfile_hg19)
# bedfile_hg38_liftovered$chromosome_hg19 <- NULL
# pho_head_cancers_sim <- merge(pho_head_cancers_sim, bedfile_hg38_liftovered, all.x = T)
# liftover_failed <- as.vector(pho_head_cancers_sim$liftover_failed)
# liftover_failed[pho_head_cancers_sim$coordinate_hg38 %in% bedfile_failed$coordinate_hg38] <- TRUE
# pho_head_cancers_sim$liftover_failed <- liftover_failed
# write.table(x = pho_head_cancers_sim,
#             file = paste0(makeOutDir(resultD = resultD), "pho_head_cancers_single_site_wpid_hg38_hg19.txt"),
#             row.names = F, col.names = T, quote = F, sep = ",")
