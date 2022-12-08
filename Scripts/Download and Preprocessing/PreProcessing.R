#Loading in Necessary Libraries (Already Installed in the "DownloadingData" Script)

library(dplyr)
library(tibble)
library(stringr)
if (!require('Rfast')) install.packages('Rfast')
library(Rfast)

### Step 1: Preprocessing Clinical Data

# Combining Dataframes Downloaded Previously using the "DownloadData" Script
load(file="~/desktop/Multimodal/Clinical/LUSC_Clin_Harm.RData")
load(file="~/desktop/Multimodal/Clinical/LUAD_Clin_Harm.RData")
LUSCclin <- bind_rows(LUSCclin, LUADclin)

#Removing unwanted clinical variables
LUSCclin <- LUSCclin %>% dplyr::select(-c(prior_treatment, race, 
                                          cigarettes_per_day, tissue_or_organ_of_origin))

#Taking Outcome Information and Numeric variables
LUSCclin_half <- LUSCclin %>% 
  dplyr::select(c(barcode, patient, shortLetterCode, treatments, days_to_last_follow_up, days_to_death, vital_status))
LUSCclin_other <- LUSCclin %>%dplyr::select(c(age_at_index, pack_years_smoked, shortest_dimension, intermediate_dimension, longest_dimension))

#Creating a new tumor volume variable by multiplying all three dimensions
LUSCclin_other <- LUSCclin_other %>% mutate(volume=shortest_dimension*intermediate_dimension*longest_dimension)
LUSCclin_other <- LUSCclin_other %>% dplyr::select(-c(shortest_dimension,intermediate_dimension,longest_dimension))

#Taking non-numeric variables (to eventually make them numeric)
LUSCclin_quantitative <- LUSCclin %>% 
  dplyr::select(c(ajcc_pathologic_stage, prior_malignancy, ajcc_pathologic_n, ajcc_pathologic_m,ajcc_pathologic_t, gender, synchronous_malignancy, name))

#Replacing Values where needed (turning non-reports into missing values
# and changing variable categories with only one observation into NAs)
for (x in 1:nrow(LUSCclin_quantitative)) {
  if(LUSCclin_quantitative$prior_malignancy[x]=="not reported") {
    LUSCclin_quantitative$prior_malignancy[x] <- NA
  }
  if(LUSCclin_quantitative$synchronous_malignancy[x]=="Not Reported") {
    LUSCclin_quantitative$synchronous_malignancy[x] <- NA
  }
  if(!is.na(LUSCclin_quantitative$ajcc_pathologic_n[x]) & LUSCclin_quantitative$ajcc_pathologic_n[x]=="NX") {
    LUSCclin_quantitative$ajcc_pathologic_n[x] <- NA
  }
  if(!is.na(LUSCclin_quantitative$ajcc_pathologic_m[x]) &LUSCclin_quantitative$ajcc_pathologic_m[x]=="MX") {
    LUSCclin_quantitative$ajcc_pathologic_m[x] <- NA
  }
  if(!is.na(LUSCclin_quantitative$ajcc_pathologic_n[x])) {
    if(LUSCclin_quantitative$ajcc_pathologic_n[x] == "N3") { #N3 has only One Observation
      LUSCclin_quantitative$ajcc_pathologic_n[x] <- NA }
  }
  if(LUSCclin_quantitative$ajcc_pathologic_t[x]=="TX") {
    LUSCclin_quantitative$ajcc_pathologic_t[x] <- NA
  }
}

#Turning factor levels into numeric
LUSCclin_quantitative <- sapply(LUSCclin_quantitative, function(x) if(class(x)=="character") { 
  as.numeric(as.factor(x)) })
LUSCclin_quantitative <- as.data.frame(LUSCclin_quantitative)

#Recombining the data frame
LUSC_all <- cbind(LUSCclin_quantitative, LUSCclin_other)
LUSC_together <- cbind(LUSC_all, LUSCclin_half)

#Separating Tumor and Non-Tumor Samples (only keeping Tumor samples)
LUSCclin_lung <- LUSC_together %>% filter(shortLetterCode=='TP')

# Limiting patients to only those in the Multimodal Analysis
barcodes_both <- LUSCclin_lung$barcode
load(file="~/desktop/Multimodal/Eligible_Barcodes/barcode_disease_mapping.Rdata") #list of 732 patients that qualify in all data categories
LUSCclin_lung$barcode <- substring(LUSCclin_lung$barcode,1,19)
LUSCclin_lung <- LUSCclin_lung %>% filter(barcode %in% barcode_disease_mapping$barcode)

#Getting censored survival times
LUSCclin_lung$censor_time1 <- LUSCclin_lung$days_to_death
for (x in 1:nrow(LUSCclin_lung)) {
  if (is.na(LUSCclin_lung$censor_time1[x])) {
    LUSCclin_lung$censor_time1[x] <- LUSCclin_lung$days_to_last_follow_up[x] }
}

#Removing those with no data or censor time of zero or one days
removal_na_indices <- which(is.na(LUSCclin_lung$censor_time1))
LUSCclin_lung <- LUSCclin_lung[-c(removal_na_indices),]
LUSCclin_lung <- LUSCclin_lung %>% filter(censor_time1 > 1)

#Taking only one distinct sample from each patient
LUSCclin_lung <- distinct(LUSCclin_lung, barcode, .keep_all = TRUE)

#Ordering the samples
LUSCclin_lung <- LUSCclin_lung[order(LUSCclin_lung$barcode),]

# Saving Clinical Data 
save(LUSCclin_lung, file="~/desktop/Multimodal/Clinical/LUSCclin_lung.Rdata")

### Step 2: Extracting lncs from mRNA data and Saving mRNA data

# Data Taken from the 'DownloadingData' script for Preprocessing
load(file="~/desktop/Multimodal/mRNA/LUSC_RNASeq_Harmonized.RData") #LUSC RNA
load(file="~/desktop/Multimodal/mRNA/LUAD_RNASeq_Harmonized.RData") #LUAD RNA

#Processing Large List of LNC RNAs Extracted from Li et al (2020): 
#Downloaded From Web Tool Associated with Paper: http://bio-bigdata.hrbmu.edu.cn/ImmLnc/
#Can download file and put into the LNCRNAS folder

#Taking Downloaded lncRNA Data and Saving It For Ease of Use (Since the Initial File is So Large)
LNCS <- read.delim("~/desktop/Multimodal/LNCRNAS/Lnc_Pathways_All.txt") # LIST OF LNC RNAS
LNCS <- LNCS %>% dplyr::select(c(lncRNA_id, lncRNA_symbol))
LNCS_unique <- unique(LNCS$lncRNA_id)
save(LNCS_unique, file = "~/desktop/Multimodal/LNCRNAS/all_lncs.Rdata")

# Using Saved lncRNA list Above to Extract LNCs From RNA-Seq Data
load(file="~/desktop/Multimodal/LNCRNAS/all_lncs.Rdata")
remove_all <- intersect(LNCS_unique, colnames(LUAD_RNAseq))
LUSC_RNAseq <- as.data.frame(LUSC_RNAseq)
LUAD_RNAseq <- as.data.frame(LUAD_RNAseq)
LUSC_RNAseq <- LUSC_RNAseq %>% dplyr::select(-c(all_of(remove_all)))
LUAD_RNAseq <- LUAD_RNAseq %>% dplyr::select(-c(all_of(remove_all)))
LUAD_RNAseq <- as.matrix(LUAD_RNAseq)
LUSC_RNAseq <- as.matrix(LUSC_RNAseq)

#Saving the mRNA data with LNC RNAs removed
save(LUSC_RNAseq, file="~/desktop/Multimodal/mRNA/LUSC_RNASeq_Harmonized_LNCSRemoved.RData")
save(LUAD_RNAseq, file="~/desktop/Multimodal/mRNA/LUAD_RNASeq_Harmonized_LNCSRemoved.RData")


### Step 3: Processing mRNA data

#Loading mRNA data
load(file="~/desktop/Multimodal/mRNA/LUAD_RNASeq_Harmonized_LNCSRemoved.RData")
load(file="~/desktop/Multimodal/mRNA/LUSC_RNASeq_Harmonized_LNCSRemoved.RData")

#Loading Clinical Data to get Correct Barcodes and Getting rid of Non-Tumor Samples
load(file="~/desktop/Multimodal/Clinical/LUSC_Clin_Harm.RData")
load(file="~/desktop/Multimodal/Clinical/LUAD_Clin_Harm.RData")
load(file="~/desktop/Multimodal/Clinical/LUSCclin_lung.Rdata")
barcodes_both <- c(LUSCclin$barcode, LUADclin$barcode)
LUSC_RNAseq <- rbind(LUSC_RNAseq, LUAD_RNAseq)
LUSC_RNAseq <- cbind(LUSC_RNAseq, barcodes_both)
tumorornot <- c(LUSCclin$shortLetterCode,LUADclin$shortLetterCode)
rownames(LUSC_RNAseq) <- tumorornot
LUSC_RNAseq  <- as.data.frame(subset(LUSC_RNAseq, rownames(LUSC_RNAseq)=="TP")) # Only tumor samples

#Filtering dataset to only include chosen patients
LUSC_RNAseq$barcodes_both <- substring(LUSC_RNAseq$barcodes_both,1,19)
LUSC_RNA_process <- LUSC_RNAseq %>% filter(barcodes_both %in% LUSCclin_lung$barcode)
LUSC_RNA_process <- distinct(LUSC_RNA_process, barcodes_both, .keep_all = TRUE)

# Removing Genes with More than 20% zeros
barcodes_both <- LUSC_RNA_process$barcodes_both
LUSC_RNA_process <- LUSC_RNA_process %>% dplyr::select(-c(barcodes_both))
limit <- nrow(LUSC_RNA_process)*.2
LUSC_RNA_process <- as.matrix(LUSC_RNA_process)
class(LUSC_RNA_process) <- "numeric"
LUSC_RNA_filt <- LUSC_RNA_process[, which(as.numeric(colSums(LUSC_RNA_process == 0)) < limit)]

#Connecting back to the barcodes
LUSC_RNA_zeros <- as.data.frame(cbind(LUSC_RNA_filt, barcodes_both))
colnames(LUSC_RNA_zeros)[ncol(LUSC_RNA_zeros)] <- "barcode"

#Taking only the multimodal patients
LUSC_RNA <- merge(LUSC_RNA_zeros, LUSCclin_lung, by="barcode")

#Saving the file before using it below
save(LUSC_RNA,file='~/desktop/Multimodal/mRNA/RNASeq_Processed.RData')


### Step 4: Processing LNC Data in the Exact Same Way

#Loading in RNA Seq and Our List of lncRNAs
load(file="~/desktop/Multimodal/LNCRNAS/all_lncs.Rdata") #Taken from Li et al (2020); see above
load(file="~/desktop/Multimodal/mRNA/LUSC_RNASeq_Harmonized.RData") #LUSC RNA
load(file="~/desktop/Multimodal/mRNA/LUAD_RNASeq_Harmonized.RData") #LUAD RNA

#Taking only the LNC RNAs from RNA Seq Data
remove_all <- intersect(LNCS_unique, colnames(LUAD_RNAseq))
LUSC_RNAseq <- as.data.frame(LUSC_RNAseq)
LUAD_RNAseq <- as.data.frame(LUAD_RNAseq)
LUSC_LNCS <- LUSC_RNAseq %>% dplyr::select(c(all_of(remove_all)))
LUAD_LNCS <- LUAD_RNAseq %>% dplyr::select(c(all_of(remove_all)))
LUAD_LNCS <- as.matrix(LUAD_LNCS)
LUSC_LNCS <- as.matrix(LUSC_LNCS)

#Clinical Data to get Barcodes and Getting Rid of Non-Tumor samples
load(file="~/desktop/Multimodal/Clinical/LUSC_Clin_Harm.RData")
load(file="~/desktop/Multimodal/Clinical/LUAD_Clin_Harm.RData")
load(file="~/desktop/Multimodal/Clinical/LUSCclin_lung.Rdata")
barcodes_both <- c(LUSCclin$barcode, LUADclin$barcode)
LUSC_LNCS <- rbind(LUSC_LNCS, LUAD_LNCS)
LUSC_LNCS <- cbind(LUSC_LNCS, barcodes_both)
tumorornot <- c(LUSCclin$shortLetterCode,LUADclin$shortLetterCode)
rownames(LUSC_LNCS) <- tumorornot
LUSC_LNCS  <- as.data.frame(subset(LUSC_LNCS, rownames(LUSC_LNCS)=="TP"))

#Filtering dataset to only include chosen patients
LUSC_LNCS$barcodes_both <- substring(LUSC_LNCS$barcodes_both,1,19)
LUSC_LNCS_process <- LUSC_LNCS %>% filter(barcodes_both %in% LUSCclin_lung$barcode)
LUSC_LNCS_process <- distinct(LUSC_LNCS_process, barcodes_both, .keep_all = TRUE)

#Removing genes with more than 20% zeros
barcodes_both <- LUSC_LNCS_process$barcodes_both
LUSC_LNCS_process <- LUSC_LNCS_process %>% dplyr::select(-c(barcodes_both))
limit <- nrow(LUSC_LNCS_process)*.2
LUSC_LNCS_process <- as.matrix(LUSC_LNCS_process)
class(LUSC_LNCS_process) <- "numeric"
LUSC_LNCS_filt <- LUSC_LNCS_process[, which(as.numeric(colSums(LUSC_LNCS_process == 0)) < limit)]

#Connecting back to the barcodes
LUSC_LNCS_zeros <- as.data.frame(cbind(LUSC_LNCS_filt, barcodes_both))
colnames(LUSC_LNCS_zeros)[ncol(LUSC_LNCS_zeros)] <- "barcode"

#Taking only the multimodal patients
LUSC_lncRNA <- merge(LUSC_LNCS_zeros, LUSCclin_lung, by="barcode")

#Saving the file before using it below
save(LUSC_lncRNA,file='~/desktop/Multimodal/LNCRNAS/LUSC_LNCS_process.RData')


### Step 5: Preprocessing Methylation Data (super slow due to large data size)

#Working with original data
  load(file="~/desktop/Multimodal/Methylation/LUSC_Methylation.RData")
  load(file="~/desktop/Multimodal/Methylation/LUAD_Methylation.RData")
  
  #Clinical Methylation data
  load(file="~/desktop/Multimodal/Methylation/LUSC_Methylation_Clin.RData")
  load(file="~/desktop/Multimodal/Methylation/LUAD_Methylation_Clin.RData")
  
  #LUSC and LUAD Row Names
  rownames(LUSC_meth) <- LUSCclin_meth[,1]
  rownames(LUAD_meth) <- LUADclin_meth[,1]
  
  #Matching up barcodes to only chosen cases
  load(file="~/desktop/Multimodal/Eligible_Barcodes/barcode_disease_mapping.Rdata")
  colnames(barcode_disease_mapping)[1] <- "rowname"
  barcode_disease_mapping$rowname <- substring(barcode_disease_mapping$rowname,1,19)
  LUSC_meth <- as.data.frame(LUSC_meth)
  LUAD_meth <- as.data.frame(LUAD_meth)
  
  #Taking only 732 chosen patients
  LUSC_meth <- rownames_to_column(LUSC_meth)
  LUAD_meth <- rownames_to_column(LUAD_meth)
  LUSC_meth$rowname <- substring(LUSC_meth$rowname,1,19)
  LUAD_meth$rowname <- substring(LUAD_meth$rowname,1,19)
  LUSC_meth <- merge(barcode_disease_mapping, LUSC_meth, by="rowname")
  LUAD_meth <- merge(barcode_disease_mapping, LUAD_meth, by="rowname")
  
  #No duplicate patient ids
  LUSC_meth <- distinct(LUSC_meth, rowname,.keep_all= TRUE)
  LUAD_meth <- distinct(LUAD_meth, rowname,.keep_all= TRUE)
  
  #Combining datasets
  LUSC_meth <- bind_rows(LUSC_meth, LUAD_meth)
  LUSC_meth <- column_to_rownames(LUSC_meth, "rowname")
  
  #Removing unneccessary name variable
  LUSC_meth <- LUSC_meth %>% dplyr::select(-c(name))
  
  #Removing SNPs
  LUSC_meth <- LUSC_meth %>% dplyr::select(!starts_with("rs"))
  
  # Remove probes with more than 20% NA values
  limit_LUSC_probes <- nrow(LUSC_meth)*.2
  LUSC_meth_nas <- LUSC_meth[ , which(colSums(is.na(LUSC_meth)) < limit_LUSC_probes)]
  
  # Remove patients with more than 20% NA values (there are none)
  limit_LUSC_patients <- ncol(LUSC_meth_nas)*.2
  LUSC_meth_nas <- LUSC_meth_nas[which(rowSums(is.na(LUSC_meth_nas)) < limit_LUSC_patients), ]
  
  # Mapping probes to cpg islands
  # CpG Data Taken Downloaded Directly from the Illumina Website, file named: "infinium-methylationepic-v-1-0-b5-manifest-file.csv"
  cpg_map <- read.csv("~/desktop/Multimodal/Methylation/CPG_ProbeData/infinium-methylationepic-v-1-0-b5-manifest-file.csv")
  cpg_map <- cpg_map %>% dplyr::select(IlmnID, UCSC_CpG_Islands_Name, Relation_to_UCSC_CpG_Island)
  probe_index2 <- which(cpg_map$Relation_to_UCSC_CpG_Island=="Island")
  probe_id2 <- cpg_map$IlmnID[probe_index2]
  LUSC_meth_cpg <- LUSC_meth_nas[ ,colnames(LUSC_meth_nas) %in% probe_id2]
  cpg_barcodes <- rownames(LUSC_meth_cpg)
  
  # Saving Results (Since I did each step individually due to slow processing time)
  save(cpg_barcodes,file= "~/desktop/Multimodal/Methylation/CPG_islands_barcodes.Rdata")
  save(LUSC_meth_cpg, file="~/desktop/Multimodal/Methylation/CPG_islands_only.Rdata")

  # Taking top 25,000 probes based on highest variance 
  # Loading same files created above so each step can be run independently
  load(file = '~/desktop/Multimodal/Methylation/CPG_islands_only.Rdata')
  load(file= "~/desktop/Multimodal/Methylation/CPG_islands_barcodes.Rdata")
  LUSC_cpg <- as.matrix(LUSC_meth_cpg)
  ranks <- colVars(LUSC_cpg, na.rm=TRUE)
  ranks_sort <- sort(ranks, index.return=TRUE, decreasing=TRUE)
  top_indexes <- ranks_sort$ix[1:25000] #top 25,000 methylation probes by variance
  LUSC_cpg_variance <- LUSC_cpg[,top_indexes]

# Adding barcode data and merging dataframes with clinical to get outcome data
  load(file="~/desktop/Multimodal/Clinical/LUSCclin_lung.Rdata")
  LUSC_cpg_variance <- cbind(LUSC_cpg_variance, cpg_barcodes)
  colnames(LUSC_cpg_variance)[ncol(LUSC_cpg_variance)] <- "barcode"
  LUSC_meth <- as.data.frame(LUSC_cpg_variance)
  LUSC_meth <- merge(LUSC_meth, LUSCclin_lung, by = "barcode")

  #Saving results
  save(LUSC_meth, file="~/desktop/Multimodal/Methylation/LUSC_meth_processed.Rdata")
  
### Step 6: miRNA Processing

# Getting miRNA Data from the Downloading Data
load(file="~/desktop/Multimodal/miRNA/LUSC_miRNA.RData")
load(file="~/desktop/Multimodal/miRNA/LUAD_miRNA.RData")
load(file="~/desktop/Multimodal/Clinical/LUSCclin_lung.Rdata")
LUSC_miRNA <- rbind(LUSC_miRNA, LUAD_miRNA)

# Fixing the Barcodes and filtering only patients with all types of data
LUSC_miRNA <- as.data.frame(LUSC_miRNA)
LUSC_miRNA <- rownames_to_column(LUSC_miRNA)
LUSC_miRNA$rowname <- str_split_fixed(LUSC_miRNA$rowname, "_", 6)[,6]
LUSC_miRNA$rowname <- substring(LUSC_miRNA$rowname,1,19)
colnames(LUSC_miRNA)[1] <- "barcode"
LUSC_miRNA <- merge(LUSC_miRNA, LUSCclin_lung, by="barcode")
LUSC_miRNA <- distinct(LUSC_miRNA, barcode, .keep_all = TRUE)

# Taking only miRNA expression values
LUSC_miRNA <- column_to_rownames(LUSC_miRNA, "barcode")
LUSC_miRNA <- LUSC_miRNA[,1:1881]

# Getting rid of features with more than 20% zeros
limit_LUSC_miRNA <- nrow(LUSC_miRNA)*.2
LUSC_miRNA <- as.data.frame(LUSC_miRNA)
LUSC_miRNA_process <- LUSC_miRNA[ , which(as.numeric(colSums(LUSC_miRNA==0)) < limit_LUSC_miRNA)]
LUSC_miRNA_process <- as.data.frame(LUSC_miRNA_process)

# Preparing Data
LUSC_miRNA_process <- rownames_to_column(LUSC_miRNA_process)
colnames(LUSC_miRNA_process)[1] <- "barcode"

#Combining with Clinical Data
LUSC_miRNA <- merge(LUSC_miRNA_process, LUSCclin_lung, by="barcode")

#Saving Data
save(LUSC_miRNA, file="~/desktop/Multimodal/miRNA/miRNA_processed.RData")
