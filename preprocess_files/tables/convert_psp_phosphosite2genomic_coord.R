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
uniprot_ensembl_map %>%
  head()

uniprot_ensembl_map2merge <- uniprot_ensembl_map %>%
  select(uniprot_entry_id, ensembl_protein_id) %>%
  unique()

# input Phosphosites from PSP ---------------------------------------------
psp_sites <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/Phospho_databases/PhosphositePlus/Feb_04_2019/Phosphorylation_site_dataset", data.table = F)

# filter to human phosphosites --------------------------------------------
psp_sites <- psp_sites %>%
  filter(ORGANISM == "human") %>%
  mutate(Phosphosite = str_split_fixed(string = MOD_RSD, pattern = "\\-", 2)[,1]) %>%
  mutate(p_coord = str_split_fixed(string = Phosphosite, pattern = '[STY]', 2)[,2])

psp_sites %>%
  head()

# merge with ensemble protein ids -----------------------------------------
psp_sites <- merge(psp_sites, uniprot_ensembl_map2merge, by.x = c("ACC_ID"), by.y = c("uniprot_entry_id"), all.x = T)
table(is.na(psp_sites$ensembl_protein_id))

# input already mapped phosphosites ---------------------------------------
pho_head_mapped_previous <- fread(input = "./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid_hg38.txt", data.table = F)

# mark whether it has genomic coordination or not -------------------------
psp_sites <- psp_sites %>%
  mutate(id = paste0(GENE, "_", Phosphosite)) %>%
  mutate(is_previously_mapped = (id %in% pho_head_mapped_previous$id))

table(psp_sites$is_previously_mapped)

psp_sites %>%
  select(id) %>%
  unique() %>%
  nrow()

write.table(x = psp_sites, file = paste0(makeOutDir(resultD = resultD), "psp_sites.txt"), quote = F, row.names = F, sep = "\t")

# query genomic coordinates based on ensembl id ---------------------------
psp_sites2map <- psp_sites %>%
  filter(is_previously_mapped == F) %>%
  filter(!is.na(ensembl_protein_id))

coords <- list()
pids = as.vector(psp_sites2map$ensembl_protein_id)
p_coords = as.vector(psp_sites2map$p_coord)
for (n in (length(coords) + 1):nrow(psp_sites2map)) {
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
df <- matrix(data = NA, nrow = nrow(psp_sites2map), ncol = ncol(coords[[1]]))
coords_names <- colnames(coords[[1]])
for (i in 1:nrow(psp_sites2map)) {
  if (!is.null(coords[[i]])) {
    if (nrow(coords[[i]]) == 1) {
      tmp <- unlist(coords[[i]])[coords_names]
      df[i,] <- matrix(tmp, nrow = 1, byrow = T)
    }
  }
}
df <- data.frame(df,stringsAsFactors=FALSE)
colnames(df) <- coords_names
psp_sites2map <- cbind(psp_sites2map, df)
psp_sites2map$is_previously_mapped <- NULL
psp_sites2map$is.solo <- NULL
psp_sites2map %>%
  head()

psp_sites2merge <- psp_sites %>%
  filter(id %in% pho_head_mapped_previous$id)
psp_sites2merge$is_previously_mapped <- NULL

psp_sites2merge <- merge(psp_sites2merge, pho_head_mapped_previous[, c("id", colnames(df))], by = c("id"), all.x = T)

psp_sites2map <- rbind(psp_sites2map, psp_sites2merge[, colnames(psp_sites2map)])
psp_sites2map <- unique(psp_sites2map)
write.table(x = psp_sites2map,
             file = paste0(makeOutDir(resultD = resultD), "psp_sites_single_site_wpid_hg38.txt"),
             row.names = F, col.names = T, quote = F, sep = "\t")
rm(coords)

psp_sites2map %>%
  select(id) %>%
  unique() %>%
  nrow()

# # lift over to hg19 -------------------------------------------------------
# psp_sites2map$chromosome <- paste0("chr", psp_sites2map$seq_region_name)
# psp_sites2map$coordinate_hg38 <- paste0(psp_sites2map$chromosome, ":", psp_sites2map$start, "-", psp_sites2map$end)
# test <- psp_sites2map[duplicated(psp_sites2map$coordinate_hg38) & !is.na(psp_sites2map$seq_region_name) ,]
# ## only keep primary sequences and no duplicated coordinates
# bedfile <- psp_sites2map[!is.na(psp_sites2map$seq_region_name) & !duplicated(psp_sites2map$coordinate_hg38) & !grepl(pattern = "HSCHR", x = psp_sites2map$seq_region_name), c("chromosome", "start", "end")]
# write.table(x = bedfile,
#             file = paste0(makeOutDir(resultD = resultD), "psp_sites_single_site_wpid_hg38_coordinates2liftover.bed"),
#             row.names = F, col.names = F, quote = F, sep = "\t")
# bedfile_hg19 <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/psp_sites_single_site_wpid_hg19.bed"), 
#                       data.table = F, col.names = c("chromosome", "start", "end"))
# bedfile_hg19$coordinate_hg19 <- paste0(bedfile_hg19$chromosome, ":", bedfile_hg19$start, "-", bedfile_hg19$end)
# ## read.delim left out a line
# bedfile_failed <- read_delim(file = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/psp_sites_single_site_wpid_hg38_coordinates2liftover.err"), 
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
# psp_sites2map <- merge(psp_sites2map, bedfile_hg38_liftovered, all.x = T)
# liftover_failed <- as.vector(psp_sites2map$liftover_failed)
# liftover_failed[psp_sites2map$coordinate_hg38 %in% bedfile_failed$coordinate_hg38] <- TRUE
# psp_sites2map$liftover_failed <- liftover_failed
# write.table(x = psp_sites2map,
#             file = paste0(makeOutDir(resultD = resultD), "psp_sites_single_site_wpid_hg38_hg19.txt"),
#             row.names = F, col.names = T, quote = F, sep = ",")
