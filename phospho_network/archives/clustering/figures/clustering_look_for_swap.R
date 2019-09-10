# Yige Wu @ WashU 2018 Jan
# clustering proteomics and phosphoproteomics data to look for tumor-normal swaps

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
clinical <- read.delim(paste0(inputD, "Specimen_Data_20161005_Yige_20180127.txt"), sep = "\t")
clinical <- data.frame(clinical)
library(readxl)
library(matrixStats)
copy_number_swap <- c("01BR001", "01BR015", "01BR017", "01BR018", "01BR025", "01BR027", "01CO008", "01OV023", "01OV030", "02OV008", "02OV015", "02OV022", "02OV023", "03BR002", "03BR004", "03BR005", "03BR010", "04OV008", "04OV011", "04OV012", "04OV013")


# inputs ------------------------------------------------------------------
samples3can <- list()
for (cancer in c("BRCA", "OV", "CO")) {
  pro <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PRO", "_formatted_normalized_replicate_averaged.txt"),
                           data.table = F)
  pro.n <- get_normal(expression = pro, clinical.m = clinical)
  pro.t <- get_tumor(expression = pro, clinical.m = clinical)
  samples3can[[cancer]] <- list()
  samples3can[[cancer]][["normal"]] <- colnames(pro.n)
  samples3can[[cancer]][["tumor"]] <- colnames(pro.t)
}

cancer <- "BRCA"; pro_BRCA <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PRO", "_formatted_normalized_replicate_averaged.txt"),
                  data.table = F)
cancer <- "OV"; pro_OV <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PRO", "_formatted_normalized_replicate_averaged.txt"),
                  data.table = F)
cancer <- "CO"; pro_CO <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PRO", "_formatted_normalized_replicate_averaged.txt"),
                data.table = F)
tmp <- merge(pro_BRCA, pro_OV, all = T)
pan3can_pro <- merge(tmp, pro_CO, all = T)

pro_BRCA.n <- get_normal(expression = pro_BRCA, clinical.m = clinical); pro_BRCA.t <- get_tumor(expression = pro_BRCA, clinical.m = clinical)
pro_OV.n <- get_normal(expression = pro_OV, clinical.m = clinical); pro_OV.t <- get_tumor(expression = pro_OV, clinical.m = clinical)
pro_CO.n <- get_normal(expression = pro_CO, clinical.m = clinical); pro_CO.t <- get_tumor(expression = pro_CO, clinical.m = clinical)

samples= c(colnames(pro_BRCA.n), colnames(pro_BRCA.t), colnames(pro_OV.n), colnames(pro_OV.t), colnames(pro_CO.n), colnames(pro_CO.t))
tmt_m = as.matrix(pan3can_pro[,samples])
samples[samples %in% colnames(pro_BRCA.n)]="#B2DF8A"
samples[samples %in% colnames(pro_BRCA.t)]="#33A02C"
samples[samples %in% colnames(pro_OV.n)]="#FDBF6F"
samples[samples %in% colnames(pro_OV.t)]="#FF7F00"
samples[samples %in% colnames(pro_CO.n)]="#CAB2D6"
samples[samples %in% colnames(pro_CO.t)]="#6A3D9A"

q99=quantile(tmt_m, probs=0.99, na.rm=T)
q1=quantile(tmt_m, probs=0.01, na.rm=T)
tmt_m2 = matrix(NA,nrow=dim(tmt_m)[1],ncol=dim(tmt_m)[2])
colnames(tmt_m2)=colnames(tmt_m)
row.names(tmt_m2)=row.names(tmt_m)
for (i in 1:nrow(tmt_m)){
  if ( (sum(tmt_m[i,][!is.na(tmt_m[i,])] > q99) + sum(tmt_m[i,][!is.na(tmt_m[i,])] < q1)) < 1){
    tmt_m2[i,]=tmt_m[i,]
  }
}


# clustering 3 cancer types individually ----------------------------------
for (cancer in c("BRCA", "OV", "CO")) {
  tmt_m3 = tmt_m2[rowSums(is.na(tmt_m2)) <= 10, colnames(tmt_m2) %in% c(samples3can[[cancer]][["normal"]], samples3can[[cancer]][["tumor"]])]
  SD=rowSds(tmt_m3, na.rm=TRUE)
  tmt_m4 = tmt_m3[SD>0.7,]
  partIDs <- sampID2partID(sampleID_vector = colnames(tmt_m4), sample_map = clinical)
  partIDs[partIDs %in% copy_number_swap ] <- paste0("*", partIDs[partIDs %in% copy_number_swap ])
  partIDs[colnames(tmt_m4) %in% c(colnames(pro_BRCA.n), colnames(pro_OV.n), colnames(pro_CO.n))] <- paste0(partIDs[colnames(tmt_m4) %in% c(colnames(pro_BRCA.n), colnames(pro_OV.n), colnames(pro_CO.n))], "_normal")
  partIDs[colnames(tmt_m4) %in% c(colnames(pro_BRCA.t), colnames(pro_OV.t), colnames(pro_CO.t))] <- paste0(partIDs[colnames(tmt_m4) %in% c(colnames(pro_BRCA.t), colnames(pro_OV.t), colnames(pro_CO.t))], "_tumor")
  
  ## make color skeme for tumor and normal
  samples_can <- colnames(tmt_m4)
  samples_can[samples_can %in% colnames(pro_BRCA.n)]="#B2DF8A"
  samples_can[samples_can %in% colnames(pro_BRCA.t)]="#33A02C"
  samples_can[samples_can %in% colnames(pro_OV.n)]="#FDBF6F"
  samples_can[samples_can %in% colnames(pro_OV.t)]="#FF7F00"
  samples_can[samples_can %in% colnames(pro_CO.n)]="#CAB2D6"
  samples_can[samples_can %in% colnames(pro_CO.t)]="#6A3D9A"
  
  resultDnow <- makeOutDir()
  pdf(paste0(resultDnow, cancer, '_PRO_naMax10_SD0.7.pdf'),width=30)
  par(oma=c(3,5,3,5))
  tmt_m4_hm = heatmap.2(tmt_m4, trace="none", na.color="white", notecol="black",
                        cexRow=0.8, cexCol=0.8, adjCol= c(0.5, NA), scale="none",dendrogram='column', ColSideColors = samples_can,
                        labRow=NA,labCol=partIDs,col=getPalette, margins=c(5,5), key=F)
  
  par(lend = 1)  
  legend("bottomleft",    # location of the legend on the heatmap plot
         legend = c("BRCA normal", "BRCA tumor", "OV normal", "OV tumor", "CRC normal", "CRC tumor"), # category labels
         col = c("#B2DF8A", "#33A02C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A"),  # color key
         lty= 1,             # line style
         lwd = 10            # line width
  )
  
  dev.off()
}

# clustering all 3 cancers all together -----------------------------------
tmt_m3 = tmt_m2[rowSums(is.na(tmt_m2)) <= 30 ,]
SD=rowSds(tmt_m3, na.rm=TRUE)
tmt_m4 = tmt_m3[SD>0.7,] # 458 observations
partIDs <- sampID2partID(sampleID_vector = colnames(tmt_m4), sample_map = clinical)
partIDs[partIDs %in% copy_number_swap ] <- paste0("*", partIDs[partIDs %in% copy_number_swap ])
partIDs[colnames(tmt_m4) %in% c(colnames(pro_BRCA.n), colnames(pro_OV.n), colnames(pro_CO.n))] <- paste0(partIDs[colnames(tmt_m4) %in% c(colnames(pro_BRCA.n), colnames(pro_OV.n), colnames(pro_CO.n))], "_normal")
partIDs[colnames(tmt_m4) %in% c(colnames(pro_BRCA.t), colnames(pro_OV.t), colnames(pro_CO.t))] <- paste0(partIDs[colnames(tmt_m4) %in% c(colnames(pro_BRCA.t), colnames(pro_OV.t), colnames(pro_CO.t))], "_tumor")

resultDnow <- makeOutDir()
pdf(paste0(resultDnow,'pan3can_PRO_naMax30_SD0.7.pdf'),width=30)
par(oma=c(3,5,3,5))
tmt_m4_hm = heatmap.2(tmt_m4, trace="none", na.color="white", notecol="black",
                        cexRow=0.8,cexCol=0.5, adjCol= c(0.5, NA), scale="none",dendrogram='column', ColSideColors = samples,
                        labRow=NA,labCol=partIDs,col=getPalette, margins=c(5,5), key=F)

par(lend = 1)  
legend("bottomleft",    # location of the legend on the heatmap plot
       legend = c("BRCA normal", "BRCA tumor", "OV normal", "OV tumor", "CRC normal", "CRC tumor"), # category labels
       col = c("#B2DF8A", "#33A02C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)

dev.off()
