# Yige Wu @ WashU 2018 Feb
# clean up clinical annotation file

# input -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
clin_raw <- read.delim(paste0(inputD, "Specimen_Data_20161005_Yige_20180127.txt"), sep = "\t")


# rid the _1 & _2 ---------------------------------------------------------
# tmp <- str_split_fixed(clin_raw$Specimen.Label, "_", 4)
# clin_clean <- clin_raw
# messed_label <- (tmp[,3] == "1" | tmp[,3] == "2")
# labels <- as.vector(clin_clean$Specimen.Label)
# labels[messed_label] <- paste0(tmp[messed_label,1], "_", tmp[messed_label, 2])
# clin_clean$Specimen.Label <- labels
# clin_clean <- clin_raw
# messed_label <- (tmp[,4] == "1" | tmp[,4] == "2")
# labels <- as.vector(clin_clean$Specimen.Label)
# labels[messed_label] <- paste0(tmp[messed_label,1], "_", tmp[messed_label, 2], "_", tmp[messed_label, 3])
# clin_clean$Specimen.Label <- labels
# clin_clean <- unique(clin_clean[,c("Participant.ID","Specimen.Label",
#                                    "tumor_normal","protein", "phosphoprotein",
#                                    "Parent.Specimen.Label", "Parent.Specimen.Title","Specimen.Comments","Pooled.","Specimen.Class" ,"Specimen.Type")])
# write_delim(clin_clean, delim = "\t", paste0(inputD, "Specimen_Data_20161005_Yige_20180209.txt"))

