# Clustering Data Using K-Means

#Loading and Installing Required Packages
if (!require('factoextra')) install.packages('factoextra')
if (!require('survival')) install.packages('survival')
library(survminer)
library(dplyr)
library(ggplot2)
library(survival)
library(factoextra)

# Taking Results From One Pipeline Run that Recorded PCA Data
load(file="~/desktop/Multimodal/MultimodalPipeline_Results_Post/CV_sigmoid0.3Zeros_AllCombinations_EarlyandLateIntegration_denoising_zeros/clustering_data.Rdata")
load(file="~/desktop/Multimodal/MultimodalPipeline_Results_Post/CV_sigmoid0.3Zeros_AllCombinations_EarlyandLateIntegration_denoising_zeros/elastic_variable.Rdata")

#Choosing a Run and Extracting Barcodes
Clustering_Data <- Clustering_Data[[5]] # feel free to select any run at random
elastic_variables <- elastic_variables[[5]] #All Non-Zero Variables
cluster_barcodes <- Clustering_Data$barcode_vec
Clustering_Data <- Clustering_Data %>% dplyr::select(-c(barcode_vec))

# Taking variables that the ElasticNet Used (All Non-Zero Variables)
Clustering_Data <- Clustering_Data %>% dplyr::select(c(all_of(elastic_variables), censor_time, vital_status))
Y <- as.data.frame(sapply(Clustering_Data, as.numeric))
Y2 <- Y %>%  dplyr::select(-c(censor_time, vital_status))

#Silhoutette Metric - Optimal Number Found at 2 Clusters
clust <- fviz_nbclust(Y2, kmeans, nstart=25, method="silhouette") + labs(subtitle="Silhouette Method")
chosen_clust <- which.max(clust$data$y)
print(paste0("The Optimal Number of Clusters is: ", chosen_clust))

# Carrying out K-Means Clustering
km.res <- kmeans(Y2, chosen_clust, iter.max=50, nstart=30)
table(km.res$cluster) #cluster labels to see how many of each

#Doing Kaplan Meier Curves Between Subtypes and Getting the P-Value
Y <- Y %>% mutate(cluster= km.res$cluster)
km_fit <- survfit(Surv(censor_time, vital_status)~factor(cluster), data=Y)
survplot <- ggsurvplot(km_fit, data=Y, ggtheme = theme_bw(), 
                       title = NULL, legend.labs = c("Group 1", "Group 2"),
                       palette = c("#E64B35B2", "#4DBBD5FF"),conf.int=T, pval.coord = c(0, 0.03),
                       pval = "log-rank test: p = 1e-9") 
survplot$plot <- survplot$plot + 
  theme(legend.text = element_text(size = 12)) + labs(x= "Time (days)", y="Survival Probability")
survplot

survplot # Viewing Plot
surv_diff <- survdiff(Surv(censor_time, vital_status)~factor(cluster), data=Y)
surv_diff # Getting Exact P-Value

# Checking Distribution of Types
load(file="~/desktop/Multimodal/Eligible_Barcodes/barcode_disease_mapping.Rdata")
Y <- Y %>% mutate(barcode = cluster_barcodes)
Checking_Type <- merge(Y, barcode_disease_mapping, by="barcode")
table(Checking_Type$name, Checking_Type$cluster)

#Labeling the Clusters and Getting Barcodes
Good_Cluster <- Checking_Type %>% filter(cluster==2)
BadPrognosis_Cluster <- Checking_Type %>% filter(cluster==1)
Good_Barcodes <- Good_Cluster$barcode
Bad_Barcodes <- BadPrognosis_Cluster$barcode

#Saving Files
save(Good_Barcodes,file= "~/desktop/Multimodal/DEA/DEA_Files/Barcodes/Good_Barcodes.Rdata")
save(Bad_Barcodes, file="~/desktop/Multimodal/DEA/DEA_Files/Barcodes/Bad_Barcodes.Rdata")
