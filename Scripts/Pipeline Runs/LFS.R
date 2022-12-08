### Linear Feature Selection Script (With no Autoencoders)
    
    ## Loading Packages (Installing if Necessary)
    if (!require(tidyverse)) install.packages('tidyverse')
    if (!require(glmnet)) install.packages('glmnet')
    if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    if (!require(SummarizedExperiment)) BiocManager::install("MultiAssayExperiment")
    if (!require(TCGAbiolinks)) BiocManager::install('TCGAbiolinks')
    if (!require(AnnotationDbi)) BiocManager::install('AnnotationDbi')
    if (!require(tibble)) install.packages('tibble')
    if (!require(purrr)) install.packages('purrr')
    if (!require(keras)) install.packages('keras')
    if (!require(survminer)) install.packages('survminer')
    if (!require(compound.Cox)) install.packages('compound.Cox')
    if (!require(MultiAssayExperiment)) BiocManager::install('MultiAssayExperiment')
    if (!require(org.Hs.eg.db)) BiocManager::install('org.Hs.eg.db')
    if (!require(tensorflow)) install.packages('tensorflow')
    if (!require(dplyr)) install.packages('dplyr')
    if (!require(ggplot2)) install.packages('ggplot2')
    if (!require(survival)) install.packages('survival')

    ## Making sure libraries are attached
    library(BiocManager)
    library(SummarizedExperiment)
    library(compound.Cox)
    library(tensorflow)
    library(survminer)
    library(ggplot2)
    library(dplyr)
    library(survival)
    library(keras)
    library(glmnet)
    library(purrr)
    library(survminer)
    library(org.Hs.eg.db)
    library(MultiAssayExperiment)
    library(AnnotationDbi)
    library(TCGAbiolinks)
    library(tidyverse)
    library(tibble)
    
    ### Initializing Values
    
    # The number of features the linear feature selection technique will choose for each modality
    n_lncs <- 300
    n_mirna <- 300
    ngenes <- 500
    nprobes <- 500
    n_clinical <- 11
    
    # The number of features LFS will reduce the feature space to
    lnc_features <- 30
    mirna_features <- 30
    meth_features <- 50
    gene_features <- 50
    
    #Total Number of Features (160 features in all)
    single_LFS_features <- meth_features+gene_features+mirna_features+lnc_features
    
    #Run Description
    runinfo <- "CV_LFS"
    one_LFS_early <- "yes"
    
    #Establishing Lists for Storing Loop Values
    overall_results_list <- list()
    total_PCA_list <- list()
    results.df_singleLFS_list <- list()
    mod_combo_results <- list()
    
    # Loading Data Frames Of Each Modality (Already Preprocessed)
    load(file="~/desktop/Multimodal/Clinical/LUSCclin_lung.Rdata") #Clinical Data
    load(file='~/desktop/Multimodal/LNCRNAS/LUSC_LNCS_process.RData') #Loading LNC dataframe
    load(file='~/desktop/Multimodal/mRNA/RNASeq_Processed.RData') # mRNA data
    load(file="~/desktop/Multimodal/miRNA/miRNA_processed.RData") # miRNA Data
    load(file="~/desktop/Multimodal/Methylation/LUSC_meth_processed.Rdata") # Methylation Data
    
    ### Splitting data into five parts to prepare for 5-fold CV:
    
    #Loading barcode data to use for splitting evenly between LUAD/LUSC below
    load(file="~/desktop/Multimodal/Eligible_Barcodes/barcode_disease_mapping.Rdata")
    
    #585 is number of training examples (80%) and 147 is number of test examples (20%)
    LUAD_indices <- which(barcode_disease_mapping$name=="Lung Adenocarcinoma") #408 luad
    LUSC_indices <- which(barcode_disease_mapping$name=="Lung Squamous Cell Carcinoma") #324 LUSC
    
    #Randomly sampling the LUAD data
    set.seed(3) #reproducible
    LUAD_indices_sample <- sample(LUAD_indices)
    
    #Saving and loading random sample indices for reference to use in other LUAD/LUSC Pipelines
    load(file="~/desktop/Multimodal/Indices/LUAD_indices_sample_5fold.Rdata") #Saved random indices used to make sure seed doesn't change results across systems
    
    #Splitting LUAD data into 5 groups for 5-fold CV:
    LUAD_indexes <- list(1:81, 82:162, 163:244, 245:326, 327:408)
    LUAD_folds <- list()
    for (x in 1:5) {
      LUAD_folds[[x]]  <- LUAD_indices_sample[LUAD_indexes[[x]]] }
    
    #LUSC Splitting 
    set.seed(3) #reproducible
    LUSC_indices_sample <- sample(LUSC_indices)
    
    #Saving indices for reference
    load(file="~/desktop/Multimodal/Indices/LUSC_indices_sample_5fold.Rdata")  #Saved random indices used to make sure seed doesn't change results across systems
    
    #Splitting data into 5 groups for 5-fold CV LUSC:
    LUSC_indexes <- list(1:65, 66:130, 131:195, 196:260, 261:324)
    LUSC_folds <- list()
    for (x in 1:5) {
      LUSC_folds[[x]]  <- LUSC_indices_sample[LUSC_indexes[[x]]]
    }
    
    ### CREATING FUNCTIONS THAT WILL BE REPEATED THROUGHOUT THIS SCRIPT:
    
    # 1. Elastic Net Function
      elastic_net_func <- function(train_data, test_data, multimodal_test_LUSC, multimodal_test_LUAD) {
      multimodal_cv <- as.matrix(train_data[1:nrow(train_data), 1:(ncol(train_data)-2)])
      multimodal_cv_test <- as.matrix(test_data[1:nrow(test_data), 1:(ncol(test_data)-2)])
      cv.fit = cv.glmnet(multimodal_cv,  Surv(multimodal_imputed$censor_time,multimodal_imputed$vital_status),
                         alpha = 0.5, # lasso: alpha = 1; ridge: alpha=0
                         family = "cox", type.measure = "C")
      
      #Predictions
      glm_pred <- predict(cv.fit, multimodal_cv_test)
      glm_c_testingdata <- Cindex(glm_pred,Surv(test_data$censor_time,test_data$vital_status))
      
      #Individual Cancer Types
      multimodal_test_LUSC_cv <-as.matrix(multimodal_test_LUSC[1:nrow(multimodal_test_LUSC), 1:(ncol(multimodal_test_LUSC)-2)])
      multimodal_test_LUAD_cv <-as.matrix(multimodal_test_LUAD[1:nrow(multimodal_test_LUAD), 1:(ncol(multimodal_test_LUAD)-2)])
      
      #Cancer Specific - LUSC
      elastic_overfit_preds_LUSC <- predict(cv.fit, multimodal_test_LUSC_cv)
      elastic_testingdata_LUSC <- Cindex(elastic_overfit_preds_LUSC, Surv(multimodal_test_LUSC$censor_time,multimodal_test_LUSC$vital_status))
      
      #Cancer Specific - LUAD
      elastic_overfit_preds_LUAD <- predict(cv.fit, multimodal_test_LUAD_cv)
      elastic_testingdata_LUAD <- Cindex(elastic_overfit_preds_LUAD, Surv(multimodal_test_LUAD$censor_time,multimodal_test_LUAD$vital_status))
      
      # Saving Results in a List
      elastic_net_list <- list(glm_c_testingdata,  elastic_testingdata_LUSC, elastic_testingdata_LUAD)
      return(elastic_net_list)
    }
    
    # 2. Feature Selection Function
    linear_featureselection_func <- function(train_data, test_data, association_matrix, time, survival, num_features) {
      
      #Linear Feature Selection Metholodology - takes top N features by lowest P-value
      set.seed(5)
      associat <- uni.selection(t.vec = time, d.vec = survival, X.mat=association_matrix,
                                P.value = 0.8, randomize=TRUE, K=5)
      associat$P <- associat$P[1:num_features]
      col_filtered <- rownames(as.data.frame(associat$P))
      
      # Saving Results for each modality
      train_data <- train_data %>% dplyr::select(all_of(col_filtered))
      test_data <- test_data %>% dplyr::select(all_of(col_filtered))
      test_data <- as.matrix(test_data)
      
      #Returning List Result
      lfs_list <- list(train_data, test_data)
      return(lfs_list)
    }
    
    # 3. PCA Visualization Function
    pca_visualization_func <- function(X, modality_name) {
      X <- merge(X, barcode_disease_mapping, by="barcode")
      
      #Taking only patients who either died before 5 years or were censored afterwards
      X_deceased <- X %>% filter(vital_status=="Dead" & censor_time < 1825) #1825 days = 5 years
      X_alive_5years <- X %>% filter(censor_time > 1825)
      X_alive_5years$vital_status <- "Alive"
      X <- rbind(X_deceased, X_alive_5years)
      
      #Creating four categories of combinations of alive/dead and NSCLC type using a loop
      disease_vital <- NULL
      for (i in 1:nrow(X)) {
        if (X$vital_status[i]=="Dead" & X$name.y[i]=="Lung Squamous Cell Carcinoma") {
          disease_vital[i] <- "LUSC - Dead"
        }
        if (X$vital_status[i]=="Alive" & X$name.y[i]=="Lung Squamous Cell Carcinoma") {
          disease_vital[i] <- "LUSC - Alive"
        }
        if (X$vital_status[i]=="Dead" & X$name.y[i]=="Lung Adenocarcinoma") {
          disease_vital[i] <- "LUAD - Dead"
        }
        if (X$vital_status[i]=="Alive" & X$name.y[i]=="Lung Adenocarcinoma") {
          disease_vital[i] <- "LUAD - Alive"
        }
      }
      
      #Recreating Data Frame and Selecting Relevant Variables
      X <- as.data.frame(cbind(X, disease_vital))
      survival_vector <- X$disease_vital
      barcode_vector <- X$barcode
      X <- X %>% dplyr::select(-c(vital_status, censor_time, TestTrain, barcode, disease_vital, name.y, volume))
      X <- X %>% dplyr::select(starts_with("V"))
      X <- scale(X)
      
      #Creating and Displaying PCA Plot
      PCA_Data <- prcomp(X)
      PCA_Plot <- ggplot(as.data.frame(PCA_Data$x), aes(x = PC1, y = PC2, col = factor(disease_vital))) + 
        geom_point() + scale_colour_manual(name = "Survival Status", values = c("#EE3377", "dark red", "cadetblue1", "blue4"))+ 
        labs(title=modality_name) + theme_classic()
      
      #Saving PCA Data
      X <- as.data.frame(X)
      PCA_func_list <- list(PCA_Data, barcode_vector, survival_vector, PCA_Plot)
      return(PCA_func_list)
    }
    
    # 4. PreProcessing Data Function
    preprocessing_func <- function(starting_data, starts_with_letter) {
      testing_data <- starting_data %>% filter(barcode %in% test_barcodes)
      training_data <- starting_data %>% filter(barcode %in% train_barcodes)
      
      #Saving Vital Status (Binary variable life/death), censored time and barcodes
      vital_vector <- training_data$vital_status
      censor_time1 <- training_data$censor_time1
      training_barcodes <- training_data$barcode
      testing_barcodes <- testing_data$barcode
      
      #Only taking LNC RNA Expression Variables (which start with E)
      training_data <- training_data %>% dplyr::select(starts_with(starts_with_letter))
      testing_data <- testing_data %>% dplyr::select(starts_with(starts_with_letter))
      training_data <- as.data.frame(sapply(training_data, as.numeric))
      testing_data  <- as.data.frame(sapply(testing_data, as.numeric))
      
      #Scaling Testing and Training Separately 
      training_data <- as.data.frame(scale(training_data))
      testing_data <- as.data.frame(scale(testing_data))
      
      # Saving Objects in a List
      preprocessing_list <- list(training_data, testing_data,
                                 vital_vector, censor_time1, training_barcodes,testing_barcodes)
      
      # Returning Variables
      return(preprocessing_list)
    }
    
    # CROSS VALIDATION LOOP: IMPLEMENTING FIVE FOLD CROSS VALIDATION:
    for (cv in 1:5) {
      
      LUAD_train_sample <- setdiff(LUAD_indices_sample, LUAD_folds[[cv]])
      LUAD_test_sample <- intersect(LUAD_indices_sample, LUAD_folds[[cv]])
      
      LUSC_train_sample <- setdiff(LUSC_indices_sample, LUSC_folds[[cv]])
      LUSC_test_sample <- intersect(LUSC_indices_sample, LUSC_folds[[cv]])
      
      #Combining Train and Test Indexes
      train_indices <- c(LUSC_train_sample, LUAD_train_sample)
      test_indices <- c(LUSC_test_sample, LUAD_test_sample)
      
      #Getting the actual barcodes in each group
      LUAD_train_barcode <- barcode_disease_mapping$barcode[LUAD_train_sample]
      LUAD_test_barcode <- barcode_disease_mapping$barcode[LUAD_test_sample]
      LUSC_train_barcode <- barcode_disease_mapping$barcode[LUSC_train_sample]
      LUSC_test_barcode <- barcode_disease_mapping$barcode[LUSC_test_sample]
      train_barcodes <- c(LUAD_train_barcode, LUSC_train_barcode)
      test_barcodes <- c(LUAD_test_barcode, LUSC_test_barcode)
      
      #Early Stopping for Training to Prevent Overfitting
      early_stop <- keras::callback_early_stopping(monitor = "val_loss", 
                                                   min_delta = 0.001, patience = 5, restore_best_weights = TRUE,
                                                   verbose = 1)
      
      # Notifying Run Number
      paste0("This is fold number: ", cv, " (out of 5)")
      
      ### Data Processing for Each Modality
      
      # LNC RNA Modality
     
      # Preprocessing
      lnc_list <- preprocessing_func(LUSC_lncRNA, "E")
      LUSC_lncRNA_process <- lnc_list[[1]]
      LUSC_lncRNA_process_test <- lnc_list[[2]]
      vital_vector <- lnc_list[[3]]
      censor_time1 <- lnc_list[[4]]
      lncrna_barcode_train <- lnc_list[[5]]
      lncrna_barcode_test <- lnc_list[[6]]
    
      #Linear Feature Selection (First Round - needed for early integration LFS; will further reduced features afterward)
      y <- (as.numeric(factor(vital_vector))-1)
      association_mat <- as.matrix(LUSC_lncRNA_process)
      lfs_lnc_list <- linear_featureselection_func(LUSC_lncRNA_process, LUSC_lncRNA_process_test, association_mat, censor_time1, y, n_lncs)
      LUSC_lncRNA_process_first <- lfs_lnc_list[[1]]
      LUSC_lncRNA_process_test_first <- lfs_lnc_list[[2]]
      
      #Extracting for altogether LFS (rather than all modalities separately)
      lncrna_single_LFS <- LUSC_lncRNA_process_first
      lncrna_single_LFS_test <- LUSC_lncRNA_process_test_first
      
      #Connecting to Barcodes
      lncrna_single_LFS$barcode <- lncrna_barcode_train
      lncrna_single_LFS_test <- as.data.frame(cbind(lncrna_single_LFS_test, lncrna_barcode_test))
      colnames(lncrna_single_LFS_test)[ncol(lncrna_single_LFS_test)] <- "barcode"

      #Linear Feature Selection:(reducing LNC directly to 30 features)
      lfs_lnc_list <- linear_featureselection_func(LUSC_lncRNA_process, LUSC_lncRNA_process_test, association_mat, censor_time1, y, lnc_features)
      LUSC_lncRNA_process <- lfs_lnc_list[[1]]
      LUSC_lncRNA_process_test <- lfs_lnc_list[[2]]
      
      # Changing Column Names for Training and Testing Data
      LUSC_lncRNA_auto <-  as.data.frame(LUSC_lncRNA_process)
      for (i in 1:ncol(LUSC_lncRNA_auto)) {
        colnames(LUSC_lncRNA_auto)[i] <- paste0("V",i) }
      LUSC_lncRNA_auto_output<-  as.data.frame(LUSC_lncRNA_process_test)
      for (i in 1:ncol(LUSC_lncRNA_auto_output)) {
        colnames(LUSC_lncRNA_auto_output)[i] <- paste0("V",i) }
      
      ### Putting Together the Dataframe
      
      #Creating patient ID columns
      LUSCclin_lung$barcode <- substring(LUSCclin_lung$barcode,1,19)
      LUSC_lncRNA_auto$barcode <- lncrna_barcode_train
      LUSC_lncRNA_auto$barcode <- substring(LUSC_lncRNA_auto$barcode,1,19)
      LUSC_lncRNA_auto <- LUSC_lncRNA_auto %>% mutate(TestTrain="Train")
      
      #Creating dataframes with miRNA, mRNA and both (can then test all of them)
      multimodal_both <- merge(LUSCclin_lung, LUSC_lncRNA_auto, by = "barcode") 
      
      # Testing Data
      LUSC_lncRNA_test <- LUSC_lncRNA_auto_output %>% mutate(TestTrain="Test")
      LUSC_lncRNA_test$barcode <- lncrna_barcode_test
      LUSC_lncRNA_test$barcode <- substring(LUSC_lncRNA_test$barcode,1,19)
      multimodal_both2 <- merge(LUSCclin_lung, LUSC_lncRNA_test, by = "barcode")
      multimodal_both <- rbind(multimodal_both, multimodal_both2)
      
      ### Setting up data for the survival models
      multimodal_prediction <- multimodal_both %>%
        dplyr::select(-c(patient, treatments, days_to_last_follow_up, days_to_death,  shortLetterCode))
      testtrain_vector <- multimodal_prediction$TestTrain
      vital_vector <- multimodal_prediction$vital_status
      censor_time <- multimodal_prediction$censor_time1
      barcode_lnc <- multimodal_prediction$barcode
      
      # Removing Variables Unnecessary for Predictions
      multimodal_prediction <- multimodal_prediction %>%
        dplyr::select(-c(vital_status, censor_time1, TestTrain, barcode))
      multimodal_prediction <- as.data.frame(lapply(multimodal_prediction, as.numeric))
      
      #Factor Variables
      multimodal_prediction$vital_status <- vital_vector
      multimodal_prediction$censor_time <- censor_time
      multimodal_prediction$vital_status <- as.factor(multimodal_prediction$vital_status)
      multimodal_prediction$TestTrain <- testtrain_vector
      multimodal_prediction$barcode <- barcode_lnc
      
      #Resplitting the data - same testing and training
      multimodal_imputed <- multimodal_prediction
      multimodal_imputed_train <- multimodal_imputed %>% filter(TestTrain=="Train")
      multimodal_imputed_test <- multimodal_imputed %>% filter(TestTrain=="Test")
      censor_time_filtered_train <- multimodal_imputed_train$censor_time
      censor_time_filtered_test <- multimodal_imputed_test$censor_time
      
      #Taking only the Desired Variables (All start with 'V')
      multimodal_imputed_train <- multimodal_imputed_train %>% dplyr::select(starts_with("v")|starts_with("V"), barcode)
      multimodal_imputed_test <- multimodal_imputed_test %>% dplyr::select(c(starts_with("v")|starts_with("V"), barcode))
      
      #Removing Volume since it is a clinical variable (but happens to start with 'V')
      multimodal_imputed_train <- multimodal_imputed_train %>% dplyr::select(-c(volume))
      multimodal_imputed_test <- multimodal_imputed_test %>% dplyr::select(-c(volume))
      
      #Readding the target variables to both test and train datasets
      multimodal_imputed_train$censor_time <- censor_time_filtered_train
      multimodal_imputed_train$vital_status <- as.numeric(multimodal_imputed_train$vital_status)
      multimodal_imputed <- multimodal_imputed_train
      multimodal_imputed_test$censor_time <- censor_time_filtered_test
      multimodal_imputed_test$vital_status <- as.numeric(multimodal_imputed_test$vital_status)
      
      # Renaming the Variables 'LNC' with numbers, as they will be added to one large
      #multimodal dataframe with all biological variables later on
      for (i in 1:((ncol(multimodal_imputed)-3))) {
        colnames(multimodal_imputed)[i] <- paste0("LNC", i) }
      for (i in 1:((ncol(multimodal_imputed_test)-3))) {
        colnames(multimodal_imputed_test)[i] <- paste0("LNC", i) }
      
      #Saving Test and Train Data for the multimodal dataset
      lnc_test <- multimodal_imputed_test
      lnc_train <- multimodal_imputed
      
      #Separating by Train/Test and by Cancer Type
      lnc_LUSC_indices <- merge(multimodal_imputed, barcode_disease_mapping, by="barcode")
      lnc_LUAD_indices_train <- which(lnc_LUSC_indices$name=="Lung Adenocarcinoma") 
      lnc_LUSC_indices_train <- which(lnc_LUSC_indices$name=="Lung Squamous Cell Carcinoma")
      lnc_LUSC_indices_test <- merge(multimodal_imputed_test, barcode_disease_mapping, by="barcode")
      lnc_LUAD_indices_test <- which(lnc_LUSC_indices_test$name=="Lung Adenocarcinoma") 
      lnc_LUSC_indices_test <- which(lnc_LUSC_indices_test$name=="Lung Squamous Cell Carcinoma")
      
      #Removing barcode variable and making sure variables are all numeric
      multimodal_imputed_test <- multimodal_imputed_test %>% dplyr::select(c(-barcode))
      multimodal_imputed <- multimodal_imputed %>% dplyr::select(c(-barcode))
      multimodal_imputed_test <- as.data.frame(sapply(multimodal_imputed_test, as.numeric))
      multimodal_imputed <- as.data.frame(sapply(multimodal_imputed, as.numeric))
      
      
      ### Elastic Net
      multimodal_test_LUSC <- multimodal_imputed_test[c(lnc_LUSC_indices_test), ]
      multimodal_test_LUAD <- multimodal_imputed_test[c(lnc_LUAD_indices_test), ]
      elasticnet_cindex_list <- elastic_net_func(multimodal_imputed, multimodal_imputed_test, multimodal_test_LUSC, multimodal_test_LUAD)
      Cindex_Test <- elasticnet_cindex_list[[1]]
      Cindex_LUSC <- elasticnet_cindex_list[[2]]
      Cindex_LUAD <- elasticnet_cindex_list[[3]]
      
      ### Results and Saving Data
      Models <- "ElasticNet"
      Data_Types <- "LNC RNAs"
      Dimension_Reduction <- "LFS"
      results.df <- as.data.frame(cbind(Models, Cindex_Test, Dimension_Reduction, Data_Types,
                                        Cindex_LUSC, Cindex_LUAD))
      results.df$Cindex_Test <- as.numeric(results.df$Cindex_Test)
      results.df$Cindex_Test <- round(results.df$Cindex_Test, 3)
      results.df_lncs <- results.df
      
      #PCA VISUALIZATION
      PCA_data_list <- pca_visualization_func(multimodal_prediction, "LNCRNA")
      LUSC_PCA_LNC <- PCA_data_list[[1]]
      PCA_barcode_LNC <- PCA_data_list[[2]]
      disease_vitalstatus_LNC <- PCA_data_list[[3]]
      PCA_data_list[[4]] # Visualizing PCA Plot
      
      
      # GENE EXPRESSION MODALITY
      
      # Preprocessing
      rna_list <- preprocessing_func(LUSC_RNA,"E")
      LUSC_RNA_process <- rna_list[[1]]
      LUSC_RNA_process_test <- rna_list[[2]]
      vital_vector <- rna_list[[3]]
      censor_time1 <- rna_list[[4]]
      mrna_barcode_train <- rna_list[[5]]
      mrna_barcode_test <- rna_list[[6]]
      
      #Linear Feature Selection
      y <- (as.numeric(factor(vital_vector))-1)
      association_mat <- as.matrix(LUSC_RNA_process)
      lfs_rna_list <- linear_featureselection_func(LUSC_RNA_process, LUSC_RNA_process_test, association_mat, censor_time1, y, ngenes)
      LUSC_RNA_process_first <- lfs_rna_list[[1]]
      LUSC_RNA_process_tests_first <- lfs_rna_list[[2]]
      
      #Extracting for altogether LFS
      mrna_single_LFS <- LUSC_RNA_process_first
      mrna_single_LFS_test <- LUSC_RNA_process_tests_first
      mrna_single_LFS$barcode <- mrna_barcode_train
      mrna_single_LFS_test <- as.data.frame(cbind(mrna_single_LFS_test, mrna_barcode_test))
      colnames(mrna_single_LFS_test)[ncol(mrna_single_LFS_test)] <- "barcode"
      
      # Reducing mRNA data directly to 50 features
      lfs_rna_list <- linear_featureselection_func(LUSC_RNA_process, as.data.frame(LUSC_RNA_process_test), association_mat, censor_time1, y, gene_features)
      LUSC_RNA_process <- lfs_rna_list[[1]]
      LUSC_RNA_process_test <- lfs_rna_list[[2]]
      
      # Changing Column Names for Training and Testing Data
      LUSC_RNA_auto <-  as.data.frame(LUSC_RNA_process)
      for (i in 1:ncol(LUSC_RNA_auto)) {
        colnames(LUSC_RNA_auto)[i] <- paste0("V",i) }
      LUSC_RNA_auto_output<-  as.data.frame(LUSC_RNA_process_test)
      for (i in 1:ncol(LUSC_RNA_auto_output)) {
        colnames(LUSC_RNA_auto_output)[i] <- paste0("V",i) }
      
      #Creating patient ID columns
      LUSC_RNA_auto$barcode <- mrna_barcode_train
      LUSC_RNA_auto$barcode <- substring(LUSC_RNA_auto$barcode,1,19)
      LUSC_RNA_auto <- LUSC_RNA_auto %>% mutate(TestTrain="Train")
      multimodal_both <- merge(LUSCclin_lung, LUSC_RNA_auto, by = "barcode") 
      
      # Testing Data
      LUSC_RNA_test <- LUSC_RNA_auto_output %>% mutate(TestTrain="Test")
      LUSC_RNA_test$barcode <- mrna_barcode_test
      LUSC_RNA_test$barcode <- substring(LUSC_RNA_test$barcode,1,19)
      multimodal_both2 <- merge(LUSCclin_lung, LUSC_RNA_test, by = "barcode")
      multimodal_both <- rbind(multimodal_both, multimodal_both2)
      
      # Putting Together the Dataframe for Testing Survival Models
      multimodal_prediction <- multimodal_both %>%
        dplyr::select(-c(patient, treatments, days_to_last_follow_up, days_to_death, shortLetterCode))
      testtrain_vector <- multimodal_prediction$TestTrain
      vital_vector <- multimodal_prediction$vital_status
      censor_time <- multimodal_prediction$censor_time1
      barcode_gene <- multimodal_prediction$barcode
      multimodal_prediction <- multimodal_prediction %>%
        dplyr::select(-c(vital_status, censor_time1, TestTrain, barcode))
      multimodal_prediction <- as.data.frame(sapply(multimodal_prediction, as.numeric))
      
      #Factor Variables
      multimodal_prediction$vital_status <- vital_vector
      multimodal_prediction$censor_time <- censor_time
      multimodal_prediction$vital_status <- as.factor(multimodal_prediction$vital_status)
      multimodal_prediction$TestTrain <- testtrain_vector
      multimodal_prediction$barcode <- barcode_gene
      
      #Resplitting the data - same testing and training
      multimodal_imputed <- multimodal_prediction
      multimodal_imputed_train <- multimodal_imputed %>% filter(TestTrain=="Train")
      multimodal_imputed_test <- multimodal_imputed %>% filter(TestTrain=="Test")
      censor_time_filtered_train <- multimodal_imputed_train$censor_time
      censor_time_filtered_test <- multimodal_imputed_test$censor_time
      
      #Filtering out Autoencoded Features
      multimodal_imputed_train <- multimodal_imputed_train %>% dplyr::select(c(starts_with("v")|starts_with("V"), barcode))
      multimodal_imputed_test <- multimodal_imputed_test %>% dplyr::select(c(starts_with("v")|starts_with("V"), barcode))
      
      #Removing volume (also starts with a 'V')
      multimodal_imputed_train <- multimodal_imputed_train %>% dplyr::select(-c(volume))
      multimodal_imputed_test <- multimodal_imputed_test %>% dplyr::select(-c(volume))
      multimodal_imputed_train$censor_time <- censor_time_filtered_train
      multimodal_imputed_train$vital_status <- as.numeric(multimodal_imputed_train$vital_status)
      multimodal_imputed <- multimodal_imputed_train
      
      multimodal_imputed_test$censor_time <- censor_time_filtered_test
      multimodal_imputed_test$vital_status <- as.numeric(multimodal_imputed_test$vital_status)
      
      # Taking dataframe for multimodal section
      for (i in 1:((ncol(multimodal_imputed)-3))) {
        colnames(multimodal_imputed)[i] <- paste0("Gene Expression", i) }
      for (i in 1:((ncol(multimodal_imputed_test)-3))) {
        colnames(multimodal_imputed_test)[i] <- paste0("Gene Expression", i) }
      gene_test <- multimodal_imputed_test
      gene_train <- multimodal_imputed
      
      gene_LUSC_indices <- merge(multimodal_imputed, barcode_disease_mapping, by="barcode")
      gene_LUAD_indices_train <- which(gene_LUSC_indices$name=="Lung Adenocarcinoma") 
      gene_LUSC_indices_train <- which(gene_LUSC_indices$name=="Lung Squamous Cell Carcinoma")
      gene_LUSC_indices_test <- merge(multimodal_imputed_test, barcode_disease_mapping, by="barcode")
      gene_LUAD_indices_test <- which(gene_LUSC_indices_test$name=="Lung Adenocarcinoma") 
      gene_LUSC_indices_test <- which(gene_LUSC_indices_test$name=="Lung Squamous Cell Carcinoma")
      multimodal_imputed_test <- multimodal_imputed_test %>% dplyr::select(c(-barcode))
      multimodal_imputed <- multimodal_imputed %>% dplyr::select(c(-barcode))
      multimodal_imputed_test <- as.data.frame(sapply(multimodal_imputed_test, as.numeric))
      multimodal_imputed <- as.data.frame(sapply(multimodal_imputed, as.numeric))
      
      # Elastic Net Model
      multimodal_test_LUSC <- multimodal_imputed_test[c(gene_LUSC_indices_test), ]
      multimodal_test_LUAD <- multimodal_imputed_test[c(gene_LUAD_indices_test), ]
      elasticnet_cindex_list <- elastic_net_func(multimodal_imputed, multimodal_imputed_test, multimodal_test_LUSC, multimodal_test_LUAD)
      Cindex_Test <- elasticnet_cindex_list[[1]]
      Cindex_LUSC <- elasticnet_cindex_list[[2]]
      Cindex_LUAD <- elasticnet_cindex_list[[3]]
      
      ### Results and Saving Data
      Data_Types <- "Gene Expression"
      results.df <- as.data.frame(cbind(Models, Cindex_Test, Dimension_Reduction, Data_Types,
                                        Cindex_LUSC, Cindex_LUAD))
      results.df$Cindex_Test <- as.numeric(results.df$Cindex_Test)
      results.df$Cindex_Test <- round(results.df$Cindex_Test, 3)
      results.df_gene <- results.df
      
      #PCA VISUALIZATION
      PCA_data_list <- pca_visualization_func(multimodal_prediction, "Gene")
      LUSC_PCA_GENE <- PCA_data_list[[1]]
      PCA_barcode_GENE <- PCA_data_list[[2]]
      disease_vitalstatus_GENE <- PCA_data_list[[3]]
      PCA_data_list[[4]] # Visualizing PCA Plot
      
      # MIRNA MODALITY
      
      # Preprocessing
      mirna_list <- preprocessing_func(LUSC_miRNA, "h")
      LUSC_miRNA_process <- mirna_list[[1]]
      LUSC_miRNA_process_test <- mirna_list[[2]]
      vital_vector <- mirna_list[[3]]
      censor_time1 <- mirna_list[[4]]
      miRNA_barcode_train <- mirna_list[[5]]
      miRNA_barcode_test <- mirna_list[[6]]
      
      #Linear Feature Selection
      y <- (as.numeric(factor(vital_vector))-1)
      association_mat <- as.matrix(LUSC_miRNA_process)
      lfs_mirna_list <- linear_featureselection_func(LUSC_miRNA_process, LUSC_miRNA_process_test, association_mat, censor_time1, y, n_mirna)
      LUSC_miRNA_process_first <- lfs_mirna_list[[1]]
      LUSC_miRNA_process_test_first <- lfs_mirna_list[[2]]

      #Extracting for altogether LFS
      mirna_single_LFS <- LUSC_miRNA_process_first
      mirna_single_LFS_test <- LUSC_miRNA_process_test_first
      mirna_single_LFS$barcode <- miRNA_barcode_train
      mirna_single_LFS_test <- as.data.frame(cbind(mirna_single_LFS_test, miRNA_barcode_test))
      colnames(mirna_single_LFS_test)[ncol(mirna_single_LFS_test)] <- "barcode"
      
      #Linear Feature Selection
      lfs_mirna_list <- linear_featureselection_func(LUSC_miRNA_process, LUSC_miRNA_process_test, association_mat, censor_time1, y, mirna_features)
      LUSC_miRNA_process <- lfs_mirna_list[[1]]
      LUSC_miRNA_process_test <- lfs_mirna_list[[2]]
      
      # Changing Column Names for Training and Testing Data
      LUSC_miRNA_auto <-  as.data.frame(LUSC_miRNA_process)
      for (i in 1:ncol(LUSC_miRNA_auto)) {
        colnames(LUSC_miRNA_auto)[i] <- paste0("V",i) }
      LUSC_miRNA_auto_output <-  as.data.frame(LUSC_miRNA_process_test)
      for (i in 1:ncol(LUSC_miRNA_auto_output)) {
        colnames(LUSC_miRNA_auto_output)[i] <- paste0("V",i) }
      
      #Creating patient ID columns
      LUSC_miRNA_auto$barcode <- miRNA_barcode_train
      LUSC_miRNA_auto$barcode <- substring(LUSC_miRNA_auto$barcode,1,19)
      LUSC_miRNA_auto <- LUSC_miRNA_auto %>% mutate(TestTrain="Train")
      multimodal_both <- merge(LUSCclin_lung, LUSC_miRNA_auto, by = "barcode") 
      
      # Testing Data
      LUSC_miRNA_test <- LUSC_miRNA_auto_output %>% mutate(TestTrain="Test")
      LUSC_miRNA_test$barcode <- miRNA_barcode_test
      multimodal_both2 <- merge(LUSCclin_lung, LUSC_miRNA_test, by = "barcode")
      multimodal_both <- rbind(multimodal_both, multimodal_both2)
      
      ### Setting up the Model Data
      multimodal_prediction <- multimodal_both %>%
        dplyr::select(-c(patient, treatments, days_to_last_follow_up, days_to_death, shortLetterCode))
      testtrain_vector <- multimodal_prediction$TestTrain
      vital_vector <- multimodal_prediction$vital_status
      censor_time <- multimodal_prediction$censor_time1
      barcode_of_mirna <- multimodal_prediction$barcode
      multimodal_prediction <- multimodal_prediction %>%
        dplyr::select(-c(vital_status, censor_time1, TestTrain, barcode))
      multimodal_prediction <- as.data.frame(sapply(multimodal_prediction, as.numeric))
      
      #Factor Variable
      multimodal_prediction$vital_status <- vital_vector
      multimodal_prediction$censor_time <- censor_time
      multimodal_prediction$vital_status <- as.factor(multimodal_prediction$vital_status)
      multimodal_prediction$TestTrain <- testtrain_vector
      multimodal_prediction$barcode <- barcode_of_mirna
      
      #Resplitting the data - same testing and training
      multimodal_imputed <- multimodal_prediction
      multimodal_imputed_train <- multimodal_imputed %>% filter(TestTrain=="Train")
      multimodal_imputed_test <- multimodal_imputed %>% filter(TestTrain=="Test")
      censor_time_filtered_train <- multimodal_imputed_train$censor_time
      censor_time_filtered_test <- multimodal_imputed_test$censor_time
      
      #Taking only autoencoded features
      multimodal_imputed_train <- multimodal_imputed_train %>% dplyr::select(c(starts_with("v")|starts_with("V"), barcode))
      multimodal_imputed_test <- multimodal_imputed_test %>% dplyr::select(c(starts_with("v")|starts_with("V"), barcode))
      multimodal_imputed_train <- multimodal_imputed_train %>% dplyr::select(-c(volume))
      multimodal_imputed_test <- multimodal_imputed_test %>% dplyr::select(-c(volume))
      
      #Adding back in target variables to train and test sets
      multimodal_imputed_train$censor_time <- censor_time_filtered_train
      multimodal_imputed_train$vital_status <- as.numeric(multimodal_imputed_train$vital_status)
      multimodal_imputed <- multimodal_imputed_train
      multimodal_imputed_test$censor_time <- censor_time_filtered_test
      multimodal_imputed_test$vital_status <- as.numeric(multimodal_imputed_test$vital_status)
      
      # Taking dataframe for multimodal section
      for (i in 1:((ncol(multimodal_imputed)-3))) {
        colnames(multimodal_imputed)[i] <- paste0("miRNA", i)
      }
      for (i in 1:((ncol(multimodal_imputed_test)-3))) {
        colnames(multimodal_imputed_test)[i] <- paste0("miRNA", i)
      }
      
      mirna_test <- multimodal_imputed_test
      mirna_train <- multimodal_imputed
      
      mirna_LUSC_indices <- merge(multimodal_imputed, barcode_disease_mapping, by="barcode")
      mirna_LUAD_indices_train <- which(mirna_LUSC_indices$name=="Lung Adenocarcinoma") 
      mirna_LUSC_indices_train <- which(mirna_LUSC_indices$name=="Lung Squamous Cell Carcinoma")
      mirna_LUSC_indices_test <- merge(multimodal_imputed_test, barcode_disease_mapping, by="barcode")
      mirna_LUAD_indices_test <- which(mirna_LUSC_indices_test$name=="Lung Adenocarcinoma") 
      mirna_LUSC_indices_test <- which(mirna_LUSC_indices_test$name=="Lung Squamous Cell Carcinoma")
      
      #Finalizing Data for Model
      multimodal_imputed_test <- multimodal_imputed_test %>% dplyr::select(c(-barcode))
      multimodal_imputed <- multimodal_imputed %>% dplyr::select(c(-barcode))
      multimodal_imputed_test <- as.data.frame(sapply(multimodal_imputed_test, as.numeric))
      multimodal_imputed <- as.data.frame(sapply(multimodal_imputed, as.numeric))
      
      ### Elastic Net
      multimodal_test_LUSC <- multimodal_imputed_test[c(mirna_LUSC_indices_test), ]
      multimodal_test_LUAD <- multimodal_imputed_test[c(mirna_LUAD_indices_test), ]
      elasticnet_cindex_list <- elastic_net_func(multimodal_imputed, multimodal_imputed_test, multimodal_test_LUSC, multimodal_test_LUAD)
      Cindex_Test <- elasticnet_cindex_list[[1]]
      Cindex_LUSC <- elasticnet_cindex_list[[2]]
      Cindex_LUAD <- elasticnet_cindex_list[[3]]
      
      ### Results and Saving Data
      Data_Types <- "miRNA"
      results.df <- as.data.frame(cbind(Models, Cindex_Test, Dimension_Reduction, Data_Types,
                                        Cindex_LUSC, Cindex_LUAD))
      results.df$Cindex_Test <- as.numeric(results.df$Cindex_Test)
      results.df$Cindex_Test <- round(results.df$Cindex_Test, 3)
      results.df_mirna <- results.df
      
      #PCA VISUALIZATION
      PCA_data_list <- pca_visualization_func(multimodal_prediction, "miRNA")
      LUSC_PCA_MIRNA <- PCA_data_list[[1]]
      PCA_barcode_MIRNA <- PCA_data_list[[2]]
      disease_vitalstatus_MIRNA <- PCA_data_list[[3]]
      PCA_data_list[[4]] # Visualizing PCA Plot
      
      
      # DNA Methylation Modality
      
      # Preprocessing
      meth_list <- preprocessing_func(LUSC_meth, "c")
      LUSC_meth_process <- meth_list[[1]]
      LUSC_meth_process_test <- meth_list[[2]]
      vital_vector <- meth_list[[3]]
      censor_time1 <- meth_list[[4]]
      meth_barcode_train <- meth_list[[5]]
      meth_barcode_test <- meth_list[[6]]
      
      #Taking only Methylation features
      LUSC_meth_process <- LUSC_meth_process %>% dplyr::select(-c(censor_time1))
      LUSC_meth_process_test <- LUSC_meth_process_test %>% dplyr::select(-c(censor_time1))
      
      #Scaling Data 
      LUSC_meth_process <- as.data.frame(scale(LUSC_meth_process))
      LUSC_meth_process_test <- as.data.frame(scale(LUSC_meth_process_test))
      
      #Separate Imputation - Training and Testing
      #Imputation loop - training data
      for (i in 1:ncol(LUSC_meth_process)) {
        LUSC_meth_process[ , i][is.na(LUSC_meth_process[ , i])] <- median(LUSC_meth_process[ , i], 
                                                                          na.rm=TRUE) }
      
      #Imputation loop - testing data
      for (i in 1:ncol(LUSC_meth_process_test)) {
        LUSC_meth_process_test[ , i][is.na(LUSC_meth_process_test[ , i])] <- median(LUSC_meth_process_test[ , i], na.rm=TRUE)
      }
      
      # Linear Feature Selection
      y <- (as.numeric(factor(vital_vector))-1)
      association_mat <- as.matrix(LUSC_meth_process)
      lfs_meth_list <- linear_featureselection_func(LUSC_meth_process, LUSC_meth_process_test, association_mat, censor_time1, y, nprobes)
      LUSC_meth_process_first <- lfs_meth_list[[1]]
      LUSC_meth_process_test_first <- lfs_meth_list[[2]]
      
      #Extracting for altogether LFS
      meth_single_LFS <- LUSC_meth_process_first
      meth_single_LFS_test <- LUSC_meth_process_test_first
      meth_single_LFS <- as.data.frame(cbind(meth_single_LFS, meth_barcode_train))
      meth_single_LFS_test <- as.data.frame(cbind(meth_single_LFS_test, meth_barcode_test))
      colnames(meth_single_LFS)[ncol(meth_single_LFS)] <- "barcode"
      colnames(meth_single_LFS_test)[ncol(meth_single_LFS_test)] <- "barcode"
      
      # Linear Feature Selection - Directly Reducing Data to 50 Features
      lfs_meth_list <- linear_featureselection_func(LUSC_meth_process, LUSC_meth_process_test, association_mat, censor_time1, y, meth_features)
      LUSC_meth_process <- lfs_meth_list[[1]]
      LUSC_meth_process_test <- lfs_meth_list[[2]]
      
      # Changing Column Names for Training and Testing Data
      LUSC_meth_auto <-  as.data.frame(LUSC_meth_process)
      for (i in 1:ncol(LUSC_meth_auto)) {
        colnames(LUSC_meth_auto)[i] <- paste0("V",i) }
      LUSC_meth_auto_output<-  as.data.frame(LUSC_meth_process_test)
      for (i in 1:ncol(LUSC_meth_auto_output)) {
        colnames(LUSC_meth_auto_output)[i] <- paste0("V",i) }
      
      ### Putting Together the Dataframe
      LUSC_meth_auto$barcode <- meth_barcode_train
      LUSC_meth_auto <- LUSC_meth_auto %>% mutate(TestTrain="Train")
      multimodal_both <- merge(LUSCclin_lung, LUSC_meth_auto, by = "barcode")
      
      # Testing Data
      LUSC_meth_test <- LUSC_meth_auto_output %>% mutate(TestTrain="Test")
      LUSC_meth_test$barcode <- meth_barcode_test
      multimodal_both2 <- merge(LUSCclin_lung, LUSC_meth_test, by = "barcode")
      multimodal_both <- rbind(multimodal_both, multimodal_both2)
      
      # Removing unused variables
      multimodal_prediction <- multimodal_both %>% dplyr::select(-c(patient, treatments, days_to_last_follow_up, days_to_death, shortLetterCode))
      
      # Setting up data for survival models
      testtrain_vector <- multimodal_prediction$TestTrain
      vital_vector <- multimodal_prediction$vital_status
      censor_time <- multimodal_prediction$censor_time1
      barcodes_of_meth <- multimodal_prediction$barcode
      multimodal_prediction <- multimodal_prediction %>%
        dplyr::select(-c(vital_status, censor_time1, TestTrain, barcode))
      multimodal_prediction <- as.data.frame(sapply(multimodal_prediction, as.numeric))
      
      #Factor Variables
      multimodal_prediction$vital_status <- vital_vector
      multimodal_prediction$censor_time <- censor_time
      multimodal_prediction$vital_status <- as.factor(multimodal_prediction$vital_status)
      multimodal_prediction$TestTrain <- testtrain_vector
      multimodal_prediction$barcode <- barcodes_of_meth
      
      # Resplitting the data - same testing and training
      multimodal_imputed <- multimodal_prediction
      multimodal_imputed_train <- multimodal_imputed %>% filter(TestTrain=="Train")
      multimodal_imputed_test <- multimodal_imputed %>% filter(TestTrain=="Test")
      censor_time_filtered_train <- multimodal_imputed_train$censor_time
      censor_time_filtered_test <- multimodal_imputed_test$censor_time
      
      #Taking only autoencoded features
      multimodal_imputed_train <- multimodal_imputed_train %>% dplyr::select(c(starts_with("v")|starts_with("V"), barcode))
      multimodal_imputed_test <- multimodal_imputed_test %>% dplyr::select(c(starts_with("v")|starts_with("V"), barcode))
      multimodal_imputed_train <- multimodal_imputed_train %>% dplyr::select(-c(volume))
      multimodal_imputed_test <- multimodal_imputed_test %>% dplyr::select(-c(volume))
      
      #Setting up target variables
      multimodal_imputed_train$censor_time <- censor_time_filtered_train
      multimodal_imputed_train$vital_status <- as.numeric(multimodal_imputed_train$vital_status)
      multimodal_imputed <- multimodal_imputed_train
      multimodal_imputed_test$censor_time <- censor_time_filtered_test
      multimodal_imputed_test$vital_status <- as.numeric(multimodal_imputed_test$vital_status)
      
      # Saving dataframe for multimodal section
      for (i in 1:((ncol(multimodal_imputed)-3))) {
        colnames(multimodal_imputed)[i] <- paste0("Methylation", i) }
      for (i in 1:((ncol(multimodal_imputed_test)-3))) {
        colnames(multimodal_imputed_test)[i] <- paste0("Methylation", i) }
      meth_test <- multimodal_imputed_test
      meth_train <- multimodal_imputed
      
      #Train
      meth_LUSC_indices <- merge(multimodal_imputed, barcode_disease_mapping, by="barcode")
      meth_LUAD_indices_train <- which(meth_LUSC_indices$name=="Lung Adenocarcinoma") 
      meth_LUSC_indices_train <- which(meth_LUSC_indices$name=="Lung Squamous Cell Carcinoma")
      #Test
      meth_LUSC_indices_test <- merge(multimodal_imputed_test, barcode_disease_mapping, by="barcode")
      meth_LUAD_indices_test <- which(meth_LUSC_indices_test$name=="Lung Adenocarcinoma") 
      meth_LUSC_indices_test <- which(meth_LUSC_indices_test$name=="Lung Squamous Cell Carcinoma")
      
      #Removing barcodes
      multimodal_imputed_test <- multimodal_imputed_test %>% dplyr::select(c(-barcode))
      multimodal_imputed <- multimodal_imputed %>% dplyr::select(c(-barcode))
      multimodal_imputed_test <- as.data.frame(sapply(multimodal_imputed_test, as.numeric))
      multimodal_imputed <- as.data.frame(sapply(multimodal_imputed, as.numeric))
      
      ### Elastic Net
      multimodal_test_LUSC <- multimodal_imputed_test[c(meth_LUSC_indices_test), ]
      multimodal_test_LUAD <- multimodal_imputed_test[c(meth_LUAD_indices_test), ]
      elasticnet_cindex_list <- elastic_net_func(multimodal_imputed, multimodal_imputed_test, multimodal_test_LUSC, multimodal_test_LUAD)
      Cindex_Test <- elasticnet_cindex_list[[1]]
      Cindex_LUSC <- elasticnet_cindex_list[[2]]
      Cindex_LUAD <- elasticnet_cindex_list[[3]]
      
      # Saving Results
      Data_Types <- "Methylation"
      results.df <- as.data.frame(cbind(Models, Cindex_Test, Dimension_Reduction, Data_Types,
                                        Cindex_LUSC, Cindex_LUAD))
      results.df$Cindex_Test <- as.numeric(results.df$Cindex_Test)
      results.df$Cindex_Test <- round(results.df$Cindex_Test, 3)
      results.df_meth <- results.df
      
      #PCA VISUALIZATION
      PCA_data_list <- pca_visualization_func(multimodal_prediction, "Methylation")
      LUSC_PCA_METH <- PCA_data_list[[1]]
      PCA_barcode_METH <- PCA_data_list[[2]]
      disease_vitalstatus_METH <- PCA_data_list[[3]]
      PCA_data_list[[4]] # Visualizing PCA Plot
      
      
      # Clinical Modality
      
      # Splitting into training and testing
      LUSC_process_test <- LUSCclin_lung %>% filter(barcode %in% test_barcodes)
      LUSC_process <- LUSCclin_lung %>% filter(barcode %in% train_barcodes)
      LUSC_process <- LUSC_process %>% mutate(TestTrain="Train")
      LUSC_process_test <- LUSC_process_test %>% mutate(TestTrain="Test")
      multimodal_both <- as.data.frame(rbind(LUSC_process, LUSC_process_test))
      
      ### Setting up Data for Survival Models
      multimodal_prediction <- multimodal_both %>%
        dplyr::select(-c(patient, treatments, days_to_last_follow_up, days_to_death,  shortLetterCode))
      testtrain_vector <- multimodal_prediction$TestTrain
      vital_vector <- multimodal_prediction$vital_status
      clinbarcodes <- multimodal_prediction$barcode
      censor_clin <- multimodal_prediction$censor_time1
      multimodal_prediction <- multimodal_prediction %>%
        dplyr::select(-c(censor_time1, vital_status,TestTrain, barcode))
      multimodal_prediction <- as.data.frame(sapply(multimodal_prediction, as.numeric))
      multimodal_prediction$vital_status <- vital_vector
      multimodal_prediction$TestTrain <- testtrain_vector
      multimodal_prediction$barcode <- clinbarcodes
      multimodal_prediction$censor_time1 <- censor_clin
      
      #Train 
      multimodal_imputed <- multimodal_prediction
      multimodal_imputed_train <- multimodal_imputed %>% filter(TestTrain=="Train")
      multimodal_imputed_test <- multimodal_imputed %>% filter(TestTrain=="Test")
      
      #Test
      multimodal_imputed <- multimodal_imputed_train
      multimodal_imputed_test <- multimodal_imputed_test %>% dplyr::select(-c(TestTrain))
      multimodal_imputed <- multimodal_imputed %>% dplyr::select(-c(TestTrain))
      
      #Putting into numeric form
      multimodal_imputed$vital_status <- as.numeric(as.factor(multimodal_imputed$vital_status))
      multimodal_imputed_test$vital_status <- as.numeric(as.factor(multimodal_imputed_test$vital_status))
      
      #Scaling 11 clinical variables separately for training and testing
      multimodal_prediction[c(1:11)] <- scale(multimodal_prediction[c(1:11)])
      multimodal_prediction <- as.data.frame(multimodal_prediction)
      multimodal_imputed[c(1:11)] <- scale(multimodal_imputed[c(1:11)])
      multimodal_imputed <- as.data.frame(multimodal_imputed)
      multimodal_imputed_test[c(1:11)] <- scale(multimodal_imputed_test[c(1:11)])
      multimodal_imputed_test <- as.data.frame(multimodal_imputed_test)
      
      # Not looping to replace barcodes below
      stop_iteration <- which(colnames(multimodal_imputed)=="barcode")
      
      #SEPARATE IMPUTATION FOR TRAINING AND TESTING
      #Imputation loop - training data
      for (i in 1:(stop_iteration-1)) {
        multimodal_imputed[ , i][is.na(multimodal_imputed[ , i])] <- median(multimodal_imputed[ , i], na.rm=TRUE) }
      
      #Imputation loop - testing data
      for (i in 1:(stop_iteration-1)) {
        multimodal_imputed_test[ , i][is.na(multimodal_imputed_test[ , i])] <- median(multimodal_imputed_test[ , i], na.rm=TRUE) }
      
      #Imputation loop - overall PCA data
      for (i in 1:(ncol(multimodal_prediction)-4)) {
        multimodal_prediction[ , i][is.na(multimodal_prediction[ , i])] <- median(multimodal_prediction[ , i], na.rm=TRUE) }
      
      # Taking dataframe for multimodal section
      clin_test <- multimodal_imputed_test
      clin_train <- multimodal_imputed
      
      # Setting Indices
      clin_LUSC_indices <- merge(multimodal_imputed, barcode_disease_mapping, by="barcode")
      clin_LUAD_indices_train <- which(clin_LUSC_indices$name.y=="Lung Adenocarcinoma") 
      clin_LUSC_indices_train <- which(clin_LUSC_indices$name.y=="Lung Squamous Cell Carcinoma")
      clin_LUSC_indices_test <- merge(multimodal_imputed_test, barcode_disease_mapping, by="barcode")
      clin_LUAD_indices_test <- which(clin_LUSC_indices_test$name.y=="Lung Adenocarcinoma") 
      clin_LUSC_indices_test <- which(clin_LUSC_indices_test$name.y=="Lung Squamous Cell Carcinoma")
      
      # Finalizing Data Frame
      multimodal_imputed_test <- multimodal_imputed_test %>% dplyr::select(c(-barcode))
      multimodal_imputed <- multimodal_imputed %>% dplyr::select(c(-barcode))
      multimodal_imputed_test <- as.data.frame(sapply(multimodal_imputed_test, as.numeric))
      multimodal_imputed <- as.data.frame(sapply(multimodal_imputed, as.numeric))
      
      ### Elastic Net
      multimodal_test_LUSC <- multimodal_imputed_test[c(clin_LUSC_indices_test), ]
      multimodal_test_LUAD <- multimodal_imputed_test[c(clin_LUAD_indices_test), ]
      elasticnet_cindex_list <- elastic_net_func(multimodal_imputed, multimodal_imputed_test, multimodal_test_LUSC, multimodal_test_LUAD)
      Cindex_Test <- elasticnet_cindex_list[[1]]
      Cindex_LUSC <- elasticnet_cindex_list[[2]]
      Cindex_LUAD <- elasticnet_cindex_list[[3]]
      
      # Saving Results
      Data_Types <- "Clinical"
      results.df <- as.data.frame(cbind(Models, Cindex_Test, Dimension_Reduction, Data_Types,
                                        Cindex_LUSC, Cindex_LUAD))
      results.df$Cindex_Test <- as.numeric(results.df$Cindex_Test)
      results.df$Cindex_Test <- round(results.df$Cindex_Test, 3)
      results.df_clin <- results.df
      
      # PCA VISUALIZATION
      X <- multimodal_prediction
      X <- merge(X, barcode_disease_mapping, by="barcode")
      X_deceased <- X %>% filter(vital_status=="Dead" & censor_time1 < 1825)
      X_alive_5years <- X %>% filter(censor_time1 > 1825)
      X_alive_5years$vital_status <- "Alive"
      X <- rbind(X_deceased, X_alive_5years)
      
      #Creating Four Groups
      disease_vital <- NULL
      for (i in 1:nrow(X)) {
        if (X$vital_status[i]=="Dead" & X$name.y[i]=="Lung Squamous Cell Carcinoma") { disease_vital[i] <- "LUSC - Dead" }
        if (X$vital_status[i]=="Alive" & X$name.y[i]=="Lung Squamous Cell Carcinoma") {disease_vital[i] <- "LUSC - Alive"}
        if (X$vital_status[i]=="Dead" & X$name.y[i]=="Lung Adenocarcinoma") {disease_vital[i] <- "LUAD - Dead"}
        if (X$vital_status[i]=="Alive" & X$name.y[i]=="Lung Adenocarcinoma") {disease_vital[i] <- "LUAD - Alive"}
      }
      
      # Saving Data and Setting Up Datafrane
      X <- as.data.frame(cbind(X, disease_vital))
      PCA_barcode_CLIN <- X$barcode
      disease_vitalstatus_CLIN <- X$disease_vital
      X <- X %>% dplyr::select(-c(vital_status, censor_time1, TestTrain, barcode, disease_vital, name.x))
      X$name.y <- as.numeric(as.factor(X$name.y))
      X <- scale(X)
      
      # Visualization
      LUSC_PCA_CLIN <- prcomp(X)
      clin_plot <- ggplot(as.data.frame(LUSC_PCA_CLIN$x), aes(x = PC1, y = PC2, col = factor(disease_vitalstatus_CLIN))) + geom_point() + 
        scale_colour_manual(name = "Survival Status", values = c("#EE3377", "dark red", "cadetblue1", "blue4"))+ 
        labs(title="Clinical") +theme_classic()
      clin_plot
      
      
      # Multimodal Analysis
      
      #Setting up the dataframe for training
      mirna_train <- mirna_train %>% dplyr::select(-c(vital_status, censor_time))
      gene_train <- gene_train %>% dplyr::select(-c(vital_status, censor_time))
      lnc_train <- lnc_train %>% dplyr::select(-c(vital_status, censor_time))
      meth_train <- meth_train %>% dplyr::select(-c(vital_status, censor_time))
      
      #Testing
      mirna_test <- mirna_test %>% dplyr::select(-c(vital_status, censor_time))
      gene_test <- gene_test %>% dplyr::select(-c(vital_status, censor_time))
      lnc_test <- lnc_test %>% dplyr::select(-c(vital_status, censor_time))
      meth_test <- meth_test %>% dplyr::select(-c(vital_status, censor_time))
      
      #Training Data Combination
      index <- which(colnames(clin_train)=="censor_time1")
      colnames(clin_train)[index] <- "censor_time"
      
      df2 <- merge(lnc_train, mirna_train, by="barcode")
      df2 <- merge(df2, gene_train, by="barcode")
      df2 <- merge(df2, meth_train, by="barcode")
      multimodal_train <- merge(df2,clin_train, by="barcode")
      
      #Testing Data Combination
      index <- which(colnames(clin_test)=="censor_time1")
      colnames(clin_test)[index] <- "censor_time"
      
      df3 <- merge(lnc_test, mirna_test, by ="barcode")
      df3 <- merge(df3, gene_test, by="barcode")
      df3 <- merge(df3, meth_test, by="barcode")
      multimodal_test <- merge(df3, clin_test,by="barcode")
      
      #Mapping Barcodes
      multimodal_LUSC_indices <- merge(multimodal_train, barcode_disease_mapping, by="barcode")
      multimodal_LUAD_indices_train <- which(multimodal_LUSC_indices$name.y=="Lung Adenocarcinoma") 
      multimodal_LUSC_indices_train <- which(multimodal_LUSC_indices$name.y=="Lung Squamous Cell Carcinoma")
      
      #Testing barcodes
      multimodal_LUSC_indices_test_multi <- merge(multimodal_test, barcode_disease_mapping, by="barcode")
      multimodal_LUAD_indices_test <- which(multimodal_LUSC_indices_test_multi$name.y=="Lung Adenocarcinoma") 
      multimodal_LUSC_indices_test <- which(multimodal_LUSC_indices_test_multi$name.y=="Lung Squamous Cell Carcinoma")
      
      #Taking barcodes for clustering
      train_code <- multimodal_train$barcode
      test_code <- multimodal_test$barcode
      
      #Removing barcodes
      multimodal_test <- multimodal_test %>% dplyr::select(-c(barcode))
      multimodal_train <- multimodal_train %>% dplyr::select(-c(barcode))
      multimodal_test <- as.data.frame(sapply(multimodal_test, as.numeric))
      multimodal_train <- as.data.frame(sapply(multimodal_train, as.numeric))
      
      ### Elastic Net
      multimodal_test_LUSC <- multimodal_test[c(multimodal_LUSC_indices_test), ]
      multimodal_test_LUAD <- multimodal_test[c(multimodal_LUAD_indices_test), ]
      elasticnet_cindex_list <- elastic_net_func(multimodal_train, multimodal_test, multimodal_test_LUSC, multimodal_test_LUAD)
      Cindex_Test <- elasticnet_cindex_list[[1]]
      Cindex_LUSC <- elasticnet_cindex_list[[2]]
      Cindex_LUAD <- elasticnet_cindex_list[[3]]
      
      ### Results 
      Data_Types <- "Multimodal"
      results.df <- as.data.frame(cbind(Models, Cindex_Test, Dimension_Reduction, Data_Types,
                                        Cindex_LUSC, Cindex_LUAD))
      results.df$Cindex_Test <- as.numeric(results.df$Cindex_Test)
      results.df$Cindex_Test <- round(results.df$Cindex_Test, 3)
      results.df_multimodal <- results.df
      
      
      # PCA Visualization
      multimodal_prediction <- rbind(multimodal_LUSC_indices, multimodal_LUSC_indices_test_multi)
      X <- multimodal_prediction
      X <- merge(X, barcode_disease_mapping, by="barcode")
      X$vital_status <- as.factor(X$vital_status)
      levels(X$vital_status) <- c("Alive", "Dead")
      X_deceased <- X %>% filter(vital_status=="Dead" & censor_time < 1825)
      X_alive_5years <- X %>% filter(censor_time > 1825)
      X_alive_5years$vital_status <- "Alive"
      X_both <- rbind(X_deceased, X_alive_5years)
      X <- X_both
      
      #Four groups
      disease_vital <- NULL
      for (i in 1:nrow(X)) {
        if (X$vital_status[i]=="Dead" & X$name.y[i]=="Lung Squamous Cell Carcinoma") { disease_vital[i] <- "LUSC - Dead" }
        if (X$vital_status[i]=="Alive" & X$name.y[i]=="Lung Squamous Cell Carcinoma") { disease_vital[i] <- "LUSC - Alive" }
        if (X$vital_status[i]=="Dead" & X$name.y[i]=="Lung Adenocarcinoma") { disease_vital[i] <- "LUAD - Dead" }
        if (X$vital_status[i]=="Alive" & X$name.y[i]=="Lung Adenocarcinoma") { disease_vital[i] <- "LUAD - Alive" } }
      X <- as.data.frame(cbind(X, disease_vital))
      PCA_barcode_MULTI <- X$barcode
      disease_vitalstatus_MULTI <- X$disease_vital
      X <- X %>% dplyr::select(-c(vital_status, censor_time, barcode, disease_vital, name.y, name.x))
      X$name <- as.numeric(as.factor(X$name))
      X <- scale(X)
      
      # Visualization
      LUSC_PCA_MULTI <- prcomp(X)
      multi_plot <- ggplot(as.data.frame(LUSC_PCA_MULTI$x), aes(x = PC1, y = PC2, col = factor(disease_vitalstatus_MULTI))) + 
        geom_point() + scale_colour_manual(name = "Survival Status",                                                                                                                                                 values = c("#EE3377", "dark red", "cadetblue1", "blue4"))  + labs(title="Multimodal") +theme_classic()
      multi_plot
      
      ### Early Integration using LFS Technique
      
      if (one_LFS_early=="yes") {
        #Training Data
        df3 <- merge(lncrna_single_LFS, mirna_single_LFS, by="barcode")
        df3 <- merge(df3, mrna_single_LFS, by="barcode")
        single_LFS <- merge(df3, meth_single_LFS, by="barcode")
        
        #Testing Data
        df4 <- merge(lncrna_single_LFS_test, mirna_single_LFS_test, by="barcode")
        df4 <- merge(df4, mrna_single_LFS_test, by="barcode")
        single_LFS_test <- merge(df4, meth_single_LFS_test, by="barcode")
        
        #Combining with Outcome Data
        outcome_data_LFS <- LUSCclin_lung %>% dplyr::select(c(barcode, vital_status, censor_time1))
        single_LFS <- merge(single_LFS, outcome_data_LFS, by="barcode")
        single_LFS_test <- merge(single_LFS_test, outcome_data_LFS, by="barcode")
        
        #Save vital status and outcome into a dataframe
        single_LFS_outcome <- single_LFS %>% dplyr::select(c(barcode, vital_status, censor_time1))
        single_LFS_outcome_test <- single_LFS_test %>% dplyr::select(c(barcode, vital_status, censor_time1))
        
        #Correct Indices Train
        full_LUSC_indices <- merge(single_LFS, barcode_disease_mapping, by="barcode")
        full_LUAD_indices_train <- which(full_LUSC_indices$name=="Lung Adenocarcinoma") 
        full_LUSC_indices_train <- which(full_LUSC_indices$name=="Lung Squamous Cell Carcinoma")
        
        #Correct Indices Test
        full_LUSC_indices_test <- merge(single_LFS_test, barcode_disease_mapping, by="barcode")
        full_LUAD_indices_test <- which(full_LUSC_indices_test$name=="Lung Adenocarcinoma") 
        full_LUSC_indices_test <- which(full_LUSC_indices_test$name=="Lung Squamous Cell Carcinoma")
        
        #Removing Barcodes
        single_LFS <- single_LFS %>% dplyr::select(-c(barcode, vital_status, censor_time1))
        single_LFS_test <- single_LFS_test %>% dplyr::select(-c(barcode, vital_status, censor_time1))
        single_LFS <- as.data.frame(sapply(single_LFS, as.numeric))
        single_LFS_test <- as.data.frame(sapply(single_LFS_test, as.numeric))

        #Linear Feature Selection
        y <- (as.numeric(factor(single_LFS_outcome$vital_status))-1)
        association_mat <- as.matrix(single_LFS)
        lfs_list <- linear_featureselection_func(single_LFS, single_LFS_test, association_mat, single_LFS_outcome$censor_time1, y, single_LFS_features)
        single_LFS <- lfs_list[[1]]
        single_LFStest <- lfs_list[[2]]
            
        #incorporating outcome data
        LUSC_LFS <- as.data.frame(cbind(single_LFS, single_LFS_outcome))
        LUSC_LFS_test <- as.data.frame(cbind(single_LFStest,single_LFS_outcome_test))
            
        #incorporating clinical data
        LUSC_LFS <- LUSC_LFS %>% dplyr::select(-c(vital_status, censor_time1))
        LUSC_LFS_test <- LUSC_LFS_test %>% dplyr::select(-c(vital_status, censor_time1))
        LUSC_single_LFS <- merge(LUSC_LFS, clin_train, by="barcode")
        LUSC_single_LFS_test <- merge(LUSC_LFS_test, clin_test, by="barcode")
            
        #Making Numeric
        LUSC_single_LFS <- LUSC_single_LFS %>% dplyr::select(-c(barcode))
        LUSC_single_LFS_test <- LUSC_single_LFS_test %>% dplyr::select(-c(barcode))
        LUSC_single_LFS <- as.data.frame(sapply(LUSC_single_LFS, as.numeric))
        LUSC_single_LFS_test <- as.data.frame(sapply(LUSC_single_LFS_test, as.numeric))
            
        #Elastic Net
        multimodal_test_LUSC <- LUSC_single_LFS_test[c(full_LUSC_indices_test), ]
        multimodal_test_LUAD <- LUSC_single_LFS_test[c(full_LUAD_indices_test), ]
        elasticnet_cindex_list <- elastic_net_func(LUSC_single_LFS, LUSC_single_LFS_test, multimodal_test_LUSC, multimodal_test_LUAD)
        Cindex_Test <- elasticnet_cindex_list[[1]]
        Cindex_LUSC <- elasticnet_cindex_list[[2]]
        Cindex_LUAD <- elasticnet_cindex_list[[3]]
            
        # Results
        Data_Types <- "Full_LFS"
        results.df <- as.data.frame(cbind(Models, Cindex_Test, Dimension_Reduction, Data_Types,
                                              Cindex_LUSC, Cindex_LUAD))
        results.df$Cindex_Test <- as.numeric(results.df$Cindex_Test)
        results.df$Cindex_Test <- round(results.df$Cindex_Test, 3)
        results.df_fullLFS <- results.df
      }
      
      
      # Joining Overall Results of Different Modalities 
      df1 <- rbind(results.df_gene, results.df_meth)
      df1 <- rbind(df1, results.df_lncs)
      df1 <- rbind(df1, results.df_mirna)
      df1 <- rbind(df1, results.df_clin)
      
      #Adding in all Combination Data and Saving Data to a List
      if (one_LFS_early=="yes") { df1 <- rbind(df1,results.df_fullLFS) }
      overall_results <- rbind(df1, results.df_multimodal)
      overall_results_list[[cv]] <- overall_results
      
      #Saving PCA Data in a list
      total_PCA <- rbind(as.data.frame(LUSC_PCA_GENE$x[,1:2]), as.data.frame(LUSC_PCA_MIRNA$x[,1:2]), 
                         as.data.frame(LUSC_PCA_METH$x[,1:2]), as.data.frame(LUSC_PCA_LNC$x[,1:2]), 
                         as.data.frame(LUSC_PCA_CLIN$x[,1:2]),as.data.frame(LUSC_PCA_MULTI$x[,1:2]))
      Genes <- rep("Gene Expression", nrow(LUSC_PCA_GENE$x))
      mirna <- rep("miRNA", nrow(LUSC_PCA_MIRNA$x))
      meth <- rep("Methylation", nrow(LUSC_PCA_METH$x))
      lnc <- rep("LNC RNA", nrow(LUSC_PCA_LNC$x))
      multi <- rep("Multimodal", nrow(LUSC_PCA_MULTI$x))
      clinical <- rep("Clinical", nrow(LUSC_PCA_CLIN$x))
      Modality <- c(Genes,mirna,meth,lnc,clinical,multi)
      total_PCA <- cbind(total_PCA, Modality)
      barcode_total <- c(PCA_barcode_GENE, PCA_barcode_MIRNA, PCA_barcode_METH,
                         PCA_barcode_LNC, PCA_barcode_CLIN, PCA_barcode_MULTI)
      totalvitalstatus <- c(disease_vitalstatus_GENE, disease_vitalstatus_MIRNA,disease_vitalstatus_METH, 
                            disease_vitalstatus_LNC,disease_vitalstatus_CLIN,disease_vitalstatus_MULTI)
      total_PCA <- cbind(total_PCA, barcode_total)
      total_PCA <- cbind(total_PCA, totalvitalstatus)
      Run_Number <- rep(paste0("Run", cv), nrow(total_PCA))
      total_PCA <- as.data.frame(cbind(total_PCA, Run_Number))
      total_PCA_list[[cv]] <- total_PCA
    }
    
    # Saving Results from Each Full Pipeline Run
    
    #Saving Results in a New Directory
    dir.create(paste0("~/desktop/Multimodal/MultimodalPipeline_Results_Post/", runinfo, "/"))
    
    #Turning Results List into Dataframe
    overall_results <- bind_rows(overall_results_list)
    
    # Saving Overall Table
    write.csv(overall_results, file=paste0("~/desktop/Multimodal/MultimodalPipeline_Results_Post/", runinfo, "/overall_results_run.csv"))
    
    #Saving PCA Data in List Form for Future Visualizations
    save(total_PCA_list, file=paste0("~/desktop/Multimodal/MultimodalPipeline_Results_Post/", 
                                     runinfo, "/PCA_Visualization_List.Rdata"))
