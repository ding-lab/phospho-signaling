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



# shared inputs -----------------------------------------------------------
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org") ## 38
filters = listFilters(ensembl)
ensembl37 = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="grch37.ensembl.org")


# gather all phosphosites -------------------------------------------------
pho_head_cancers <- NULL
for (cancer in c("BRCA", "OV", "CO", "UCEC", "CCRCC")) {
  pho <- fread(paste0(cptac_sharedD, cancer, "/", prefix[cancer], "_", "PHO", "_formatted.txt"),
               data.table = F)
  pho_head <- formatPhosphosite(phosphosite_vector = as.vector(pho$Phosphosite), gene_vector = as.vector(pho$Gene))
  pho_head_cancers <- rbind(pho_head_cancers, pho_head)
  pho_head_cancers <- unique(pho_head_cancers)
}
pho_head_cancers$offset <- str_split_fixed(string = pho_head_cancers$SUB_MOD_RSD, pattern = '[STY]', 3)[,2]
pho_head_cancers$offset.2 <- str_split_fixed(string = pho_head_cancers$SUB_MOD_RSD, pattern = '[STY]', 3)[,3]
pho_head_cancers_sim <- pho_head_cancers[pho_head_cancers$offset.2 == "",]

# convert refseq ids to ensembl ids ---------------------------------------
pho_head_cancers_sim$refseq_peptide <- str_split_fixed(string = pho_head_cancers_sim$transcript, pattern = "\\.", 2)[,1]
pho_head_cancers_sim$with_refseq_peptide_predicted <- grepl(x = as.vector(pho_head_cancers_sim$refseq_peptide), pattern = "XP")
refseq_peptides <- unique(pho_head_cancers_sim$refseq_peptide[!pho_head_cancers_sim$with_refseq_peptide_predicted])
refseqs_mapTab = getBM(attributes = c("ensembl_peptide_id", "refseq_peptide", "hgnc_symbol"),
                       filters = c("refseq_peptide"),
                       values = refseq_peptides, mart = ensembl, uniqueRows=T)
refseqs_mapTab$refseq_peptide_predicted <- refseqs_mapTab$refseq_peptide
refseq_peptides_predicted <- unique(pho_head_cancers_sim$refseq_peptide[pho_head_cancers_sim$with_refseq_peptide_predicted])
refseqs_predicted_mapTab = getBM(attributes = c("ensembl_peptide_id", "refseq_peptide", "hgnc_symbol", "refseq_peptide_predicted"),
                                 filters = c("refseq_peptide_predicted"),
                                 values = refseq_peptides_predicted, mart = ensembl, uniqueRows=T)
refseqs_mapTab2m <- rbind(refseqs_mapTab, refseqs_predicted_mapTab)[,c("refseq_peptide_predicted", "ensembl_peptide_id")]
refseqs_mapTab2m <- refseqs_mapTab2m[!duplicated(refseqs_mapTab2m$refseq_peptide_predicted),]
pho_head_cancers_sim <- merge(pho_head_cancers_sim, refseqs_mapTab2m, by.x = c("refseq_peptide"), by.y = c("refseq_peptide_predicted"), all.x = T)

# query genomic coordinates based on ensembl id ---------------------------
coords <- list()
pids = as.vector(pho_head_cancers_sim$ensembl_peptide_id)
offsets = as.vector(pho_head_cancers_sim$offset)
for (n in 1:nrow(pho_head_cancers_sim)) {
  if (!is.na(pids[n])) {
    r <- GET(paste("http://rest.ensembl.org", paste0("/map/translation/", pids[n], "/", offsets[n], "..", offsets[n], "?"), sep = ""), content_type("application/json"))
    stop_for_status(r)
    tmp <- fromJSON(toJSON(content(r)))$mappings
    print(n)
    coords[[n]] <- tmp
  } 
}
# saveRDS(object = coords, file = paste0(makeOutDir(), "coords.RDS"))
df <- matrix(data = NA, nrow = nrow(pho_head_cancers_sim), ncol = ncol(coords[[1]]))
for (i in 1:nrow(pho_head_cancers_sim)) {
  if (!is.null(coords[[i]])) {
    if (nrow(coords[[i]]) == 1) {
      df[i,] <- matrix(unlist(coords[[i]]), nrow = 1, byrow = T)
    }
  }
}
df <- data.frame(df,stringsAsFactors=FALSE)
colnames(df) <- colnames(coords[[1]])
pho_head_cancers_sim <- cbind(pho_head_cancers_sim, df)
swrite.table(x = pho_head_cancers_sim,
            file = paste0(makeOutDir(resultD = resultD), "pho_head_cancers_single_site_wpid_hg38.txt"),
            row.names = F, col.names = T, quote = F, sep = ",")
rm(coords)


# lift over to hg19 -------------------------------------------------------
pho_head_cancers_sim$chromosome <- paste0("chr", pho_head_cancers_sim$seq_region_name)
pho_head_cancers_sim$coordinate_hg38 <- paste0(pho_head_cancers_sim$chromosome, ":", pho_head_cancers_sim$start, "-", pho_head_cancers_sim$end)
test <- pho_head_cancers_sim[duplicated(pho_head_cancers_sim$coordinate_hg38) & !is.na(pho_head_cancers_sim$seq_region_name) ,]
## only keep primary sequences and no duplicated coordinates
bedfile <- pho_head_cancers_sim[!is.na(pho_head_cancers_sim$seq_region_name) & !duplicated(pho_head_cancers_sim$coordinate_hg38) & !grepl(pattern = "HSCHR", x = pho_head_cancers_sim$seq_region_name), c("chromosome", "start", "end")]
write.table(x = bedfile,
            file = paste0(makeOutDir(resultD = resultD), "pho_head_cancers_single_site_wpid_hg38_coordinates2liftover.bed"),
            row.names = F, col.names = F, quote = F, sep = "\t")
bedfile_hg19 <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid_hg19.bed"), 
                      data.table = F, col.names = c("chromosome", "start", "end"))
bedfile_hg19$coordinate_hg19 <- paste0(bedfile_hg19$chromosome, ":", bedfile_hg19$start, "-", bedfile_hg19$end)
## read.delim left out a line
bedfile_failed <- read_delim(file = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid_hg38_coordinates2liftover.err"), 
                             comment = c("#"), col_names = c("chromosome", "start", "end"), delim = "\t")
bedfile_failed$liftover_failed <- TRUE
bedfile_failed$coordinate_hg38 <- paste0(bedfile_failed$chromosome, ":", bedfile_failed$start, "-", bedfile_failed$end)
bedfile$coordinate_hg38 <- paste0(bedfile$chromosome, ":", bedfile$start, "-", bedfile$end)
bedfile$liftover_failed <- (bedfile$coordinate_hg38 %in% bedfile_failed$coordinate_hg38)
bedfile_hg38_liftovered <- bedfile[!bedfile$liftover_failed,]
nrow(bedfile_hg38_liftovered) == nrow(bedfile_hg19)
colnames(bedfile_hg19) <- c("chromosome_hg19", "start_hg19", "end_hg19", "coordinate_hg19")
bedfile_hg38_liftovered <- cbind(bedfile_hg38_liftovered, bedfile_hg19)
bedfile_hg38_liftovered$chromosome_hg19 <- NULL
pho_head_cancers_sim <- merge(pho_head_cancers_sim, bedfile_hg38_liftovered, all.x = T)
liftover_failed <- as.vector(pho_head_cancers_sim$liftover_failed)
liftover_failed[pho_head_cancers_sim$coordinate_hg38 %in% bedfile_failed$coordinate_hg38] <- TRUE
pho_head_cancers_sim$liftover_failed <- liftover_failed
write.table(x = pho_head_cancers_sim,
            file = paste0(makeOutDir(resultD = resultD), "pho_head_cancers_single_site_wpid_hg38_hg19.txt"),
            row.names = F, col.names = T, quote = F, sep = ",")
