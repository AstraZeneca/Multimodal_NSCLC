### Parameters to Set: (sigmoid, 0.3 proportion of zeros and denoising zeros autoencoder
# were chosen for this study)

# Choose your autoencoder type 
autoencoder_type <- "denoising_zeros"
# Other Autoencoder Option: "denoising_gaussian",

#Choose your activation function (could also choose "relu" function):
activation_function_set <- c("sigmoid" #,"tanh","relu"
)

#Percentage of matrix to replace with zeros (fights overfitting by the autoencoder)
zeros_percentage_set <- c(0.3 #, 0.2, 0.4
                          )
gaussian_set <- c(0.1, 0.5, 1)


# Looping Through The Chosen Activation Functions and Chosen Hyperparameters (If you want to test multiple combinations)

for (u in 1:length(activation_function_set)) {
  for (z in 1:length(zeros_percentage_set)) {
    
    ## Loading Packages (Installing if Necessary)
    if (!require(tidyverse)) install.packages('tidyverse')
    if (!require(glmnet)) install.packages('glmnet')
    if (!require(BiocManager, quietly = TRUE)) install.packages("BiocManager")
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
    library(ggplot2)
    library(BiocManager)
    library(SummarizedExperiment)
    library(compound.Cox)
    library(tensorflow)
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
    
    # The number of features the denoising autoencoder will reduce the feature space to
    lnc_features <- 30
    mirna_features <- 30
    meth_features <- 50
    gene_features <- 50
    
    #Total Number of Features (150 features in all)
    fullAE_features <- meth_features+gene_features+mirna_features+lnc_features
    
    # Whether You Want All Combinations of Features (Or Just Uni- and Multimodal)
    all_modalities_combinations_late <- "yes"
    one_AE_early <- "yes"
    one_AE_early_allcombos <- "yes"
    
    # Varying the Training Datasets If Desired (Still Tests on Both Datasets)
    only_train_on_LUSC <- "no"
    only_train_on_LUAD <- "no"
    
    # Setting Chosen Parameters for Each Run
    activation_function <- activation_function_set[u]
    if (autoencoder_type == "denoising_gaussian") { gaussian_sd <- gaussian_set[z] }
    if (autoencoder_type == "denoising_zeros") { zeros_percentage <- zeros_percentage_set[z] }
    
    #Run Description (Will Become Title of Folder Data is Saved in)
    label <- paste0(activation_function, zeros_percentage, "Zeros_AllCombinations_EarlyandLateIntegration_")
    runinfo <- paste0("CV_", label, autoencoder_type)
    
    #Establishing Lists for Storing Loop Values
    overall_results_list <- list()
    feat_importance_list <- list()
    autoencoder_list <- list()
    total_PCA_list <- list()
    elastic_variables <- list()
    Clustering_Data <- list()
    results.df_fullAE_list <- list()
    
    # Loading Data Frames Of Each Modality (Already Preprocessed in the Preprocessing Script)
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
    #save(LUAD_indices_sample, file="~/desktop/Multimodal/Indices/LUAD_indices_sample_5fold.Rdata")
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
    #save(LUAD_indices_sample, file="~/desktop/Multimodal/Indices/LUAD_indices_sample_5fold.Rdata")
    load(file="~/desktop/Multimodal/Indices/LUSC_indices_sample_5fold.Rdata")  #Saved random indices used to make sure seed doesn't change results across systems
    
    #Splitting data into 5 groups for 5-fold CV LUSC:
    LUSC_indexes <- list(1:65, 66:130, 131:195, 196:260, 261:324)
    LUSC_folds <- list()
    for (x in 1:5) {
      LUSC_folds[[x]]  <- LUSC_indices_sample[LUSC_indexes[[x]]]
    }

    ### CREATING FUNCTIONS THAT WILL BE REPEATED THROUGHOUT THIS SCRIPT:
    
    # 1. Autoencoder Functions
      # A. Denoising Zeros Autoencoder
    
    # Function to Add Random Zeros to Input Matrix Using the Chosen 'Zero Percentage'
    corrupt_with_ones <- function(x) {
      n_to_sample <- floor(length(x) * zeros_percentage)
      elements_to_corrupt <- sample(seq_along(x), n_to_sample, replace = FALSE)
      x[elements_to_corrupt] <- 0
      return(x) }
     
    # Autoencoder Function - takes in training dataset, testing dataset and number of features to reduce dataset to
     denoising_zeros_func <- function(train_data, test_data, num_features) {
        #Reproducible
        set.seed(5) 
       
        #Creating a corrupted matrix with zeros in random places
        inputs_currupted_ones <- train_data %>%
          as.data.frame() %>%
          purrr::map_df(corrupt_with_ones)
        features <- as.matrix(train_data)
        inputs_currupted_ones <- as.matrix(inputs_currupted_ones)
        test_data <- as.matrix(test_data)
        
        #Setting up the model with only one layer (bottleneck)
        model1 <- keras_model_sequential()
        model1 %>%
          layer_dense(units = num_features, activation = activation_function, input_shape = ncol(inputs_currupted_ones), name = "BottleNeck") %>%
          layer_dense(units = ncol(inputs_currupted_ones))
        
        # Compiling Model
        model1 %>% keras::compile(
          loss = "mean_squared_error",
          optimizer = optimizer_adam(lr=0.001))
        
        #Training the Model with Early Stopping and 10% Validation Data
        history <- model1 %>% keras::fit(
          x = inputs_currupted_ones, y = features,
          epochs = 100, validation_split=0.1,
          callbacks = list(early_stop))
        
        #Taking the Bottleneck Layer and Evaluating RMSE
        train_RMSE <- evaluate(model1, features,inputs_currupted_ones)
        intermediate_layer_model1 <- keras_model(inputs = model1$input, outputs = get_layer(model1, "BottleNeck")$output)
        test_RMSE <- evaluate(model1,test_data, test_data)
        
        #Saving Train and Test Results
        denoising_zeros_list <- list(predict(intermediate_layer_model1, inputs_currupted_ones), 
                                     predict(intermediate_layer_model1, test_data), train_RMSE, test_RMSE)
        return(denoising_zeros_list)
    }
    
     # B. Denoising Gaussian Autoencoder
     denoising_gaussian_func <- function(train_data, test_data, num_features) {
       
       #Reproducible
       set.seed(5) 
       
       #Creating a corrupted matrix with zeros in random places
       inputs_currupted_ones <- train_data %>%
         as.data.frame() %>%
         purrr::map_df(corrupt_with_ones)
       features <- as.matrix(train_data)
       inputs_currupted_ones <- as.matrix(inputs_currupted_ones)
       test_data <- as.matrix(test_data)
       
       #Setting up the model with only one layer (bottleneck)
       model1 <- keras_model_sequential()
       model1 %>%
         layer_gaussian_noise(stddev = gaussian_sd, input_shape = ncol(features)) %>% #adding noise
         layer_dense(units = num_features, activation = activation_function, input_shape = ncol(inputs_currupted_ones), name = "BottleNeck") %>%
         layer_dense(units = ncol(inputs_currupted_ones))
       
       # Compiling Model
       model1 %>% keras::compile(
         loss = "mean_squared_error",
         optimizer = optimizer_adam(lr=0.001))
       
       #Training the Model with Early Stopping and 10% Validation Data
       history <- model1 %>% keras::fit(
         x = inputs_currupted_ones, y = features,
         epochs = 100, validation_split=0.1,
         callbacks = list(early_stop))
       
       #Taking the Bottleneck Layer and Evaluating RMSE
       train_RMSE <- evaluate(model1, features,inputs_currupted_ones)
       intermediate_layer_model1 <- keras_model(inputs = model1$input, outputs = get_layer(model1, "BottleNeck")$output)
       test_RMSE <- evaluate(model1,test_data, test_data)
       
       #Saving Train and Test Results
       denoising_gaussian_list <- list(predict(intermediate_layer_model1, inputs_currupted_ones), 
                                    predict(intermediate_layer_model1, test_data), train_RMSE, test_RMSE)
       return(denoising_gaussian_list)
     }
    
    # 2. Elastic Net Function
     elastic_net_func <- function(train_data, test_data, multimodal_test_LUSC, multimodal_test_LUAD) {
     multimodal_cv <- as.matrix(train_data[1:nrow(train_data), 1:(ncol(train_data)-2)])
     multimodal_cv_test <- as.matrix(test_data[1:nrow(test_data), 1:(ncol(test_data)-2)])
     cv.fit = cv.glmnet(multimodal_cv,  Surv(multimodal_imputed$censor_time,multimodal_imputed$vital_status),
                        alpha = 0.5, # lasso: alpha = 1; ridge: alpha=0
                        family = "cox", type.measure = "C")
     est.coef <- coef(cv.fit, s = cv.fit$lambda.min)
     
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
     elastic_net_list <- list(glm_c_testingdata,  elastic_testingdata_LUSC, elastic_testingdata_LUAD, est.coef)
     return(elastic_net_list)
     }
     
    # 3. Feature Selection Function
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
     
    # 4. PCA Visualization Function
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
     
     # 5. PreProcessing Data Function
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
     
    ###  CROSS VALIDATION LOOP: IMPLEMENTING FIVE FOLD CROSS VALIDATION:
    for (cv in 1:5) {
      
      # Setting A Unique Train/Test Split for Each Cross Validation Fold
      LUAD_train_sample <- setdiff(LUAD_indices_sample, LUAD_folds[[cv]])
      LUAD_test_sample <- intersect(LUAD_indices_sample, LUAD_folds[[cv]])
      LUSC_train_sample <- setdiff(LUSC_indices_sample, LUSC_folds[[cv]])
      LUSC_test_sample <- intersect(LUSC_indices_sample, LUSC_folds[[cv]])
      
      #Combining Train and Test Indexes for LUSC and LUAD
      train_indices <- c(LUSC_train_sample, LUAD_train_sample)
      test_indices <- c(LUSC_test_sample, LUAD_test_sample)
      
      #Extracting the actual barcodes in each group
      LUAD_train_barcode <- barcode_disease_mapping$barcode[LUAD_train_sample]
      LUAD_test_barcode <- barcode_disease_mapping$barcode[LUAD_test_sample]
      LUSC_train_barcode <- barcode_disease_mapping$barcode[LUSC_train_sample]
      LUSC_test_barcode <- barcode_disease_mapping$barcode[LUSC_test_sample]
      train_barcodes <- c(LUAD_train_barcode, LUSC_train_barcode)
      test_barcodes <- c(LUAD_test_barcode, LUSC_test_barcode)
      
      # Varying Training Data If Desired
      # Only LUAD
      if (only_train_on_LUAD=="yes") { 
      train_indices <- c(LUAD_train_sample) 
      train_barcodes <- c(LUAD_train_barcode) }
      
      # Only LUSC
      if (only_train_on_LUSC=="yes") { 
      train_indices <- c(LUSC_train_sample) 
      train_barcodes <- c(LUSC_train_barcode) }
      
      #Early Stopping Mechanism for Training to Prevent Overfitting
      early_stop <- keras::callback_early_stopping(monitor = "val_loss", min_delta = 0.001, 
                                                   patience = 5, restore_best_weights = TRUE, verbose = 1)
      
      # Notifying Run Number
      print(paste0("This is fold number: ", cv, " (out of 5)"))
      
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
      
      #Linear Feature Selection
      y <- (as.numeric(factor(vital_vector))-1)
      association_mat <- as.matrix(LUSC_lncRNA_process)
      lfs_lnc_list <- linear_featureselection_func(LUSC_lncRNA_process, LUSC_lncRNA_process_test, association_mat, censor_time1, y, n_lncs)
      LUSC_lncRNA_process <- lfs_lnc_list[[1]]
      LUSC_lncRNA_process_test <- lfs_lnc_list[[2]]
      
      #Extracting for altogether AutoEncoder (rather than all modalities separately)
      lncrna_fullAE <- LUSC_lncRNA_process
      lncrna_fullAE_test <- LUSC_lncRNA_process_test
      
      #Connecting to Barcodes
      lncrna_fullAE$barcode <- lncrna_barcode_train
      lncrna_fullAE_test <- as.data.frame(cbind(lncrna_fullAE_test, lncrna_barcode_test))
      colnames(lncrna_fullAE_test)[ncol(lncrna_fullAE_test)] <- "barcode"
      
      #AUTOENCODER IMPLEMENTATIONL
      
      ### Denoising Zeros Autoencoder
      if (autoencoder_type == "denoising_zeros") { 
        denoising_zeros_list <- denoising_zeros_func(LUSC_lncRNA_process, LUSC_lncRNA_process_test, lnc_features)
        LUSC_lncRNA_auto <- denoising_zeros_list[[1]]
        LUSC_lncRNA_auto_output <- denoising_zeros_list[[2]]
        train_RMSE <- denoising_zeros_list[[3]]
        test_RMSE <- denoising_zeros_list[[4]]
        }
      
      ### Other Autoencoder Type: Denoising (Gaussian), 
      if (autoencoder_type=="denoising_gaussian") {
        denoising_gaussian_list <- denoising_gaussian_func(LUSC_lncRNA_process, LUSC_lncRNA_process_test, lnc_features) 
        LUSC_lncRNA_auto <- denoising_gaussian_list[[1]]
        LUSC_lncRNA_auto_output <- denoising_gaussian_list[[2]] 
        train_RMSE <- denoising_gaussian_list[[3]]
        test_RMSE <- denoising_gaussian_list[[4]] 
        }
      
      ### Putting Together the Dataframe
      
      # Turning into dataframes
      LUSC_lncRNA_auto <- as.data.frame(LUSC_lncRNA_auto)
      LUSCclin_lung$barcode <- substring(LUSCclin_lung$barcode,1,19)
      
      # Creating patient ID columns
      LUSC_lncRNA_auto$barcode <- lncrna_barcode_train
      LUSC_lncRNA_auto$barcode <- substring(LUSC_lncRNA_auto$barcode,1,19)
      LUSC_lncRNA_auto <- LUSC_lncRNA_auto %>% mutate(TestTrain="Train")
      
      #C reating dataframes with miRNA, mRNA and both (can then test all of them)
      multimodal_both <- merge(LUSCclin_lung, LUSC_lncRNA_auto, by = "barcode") 
      
      # Testing Data
      LUSC_lncRNA_auto_output <- as.data.frame(LUSC_lncRNA_auto_output)
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
      
      #Taking only the Autoencoder Variables (All start with 'V')
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
      multimodal_imputed_test <- multimodal_imputed_test %>% dplyr::select(-c(barcode))
      multimodal_imputed <- multimodal_imputed %>% dplyr::select(-c(barcode))
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
      Dimension_Reduction <- autoencoder_type
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
      
      ### Autoencoder Performance Saved through AE_RMSE Vector
        AE_RMSE <- NULL
        AE_RMSE_train <- NULL
        AE_RMSE[1] <- test_RMSE
        AE_RMSE_train[1] <- train_RMSE

        
      # GENE EXPRESSION MODALITY
      
      # Preprocessing
      rna_list <- preprocessing_func(LUSC_RNA,"E")
      LUSC_RNA_process <- rna_list[[1]]
      LUSC_RNA_process_test <- rna_list[[2]]
      vital_vector <- rna_list[[3]]
      censor_time1 <- rna_list[[4]]
      rna_barcode_train <- rna_list[[5]]
      rna_barcode_test <- rna_list[[6]]
      
      #Linear Feature Selection
      y <- (as.numeric(factor(vital_vector))-1)
      association_mat <- as.matrix(LUSC_RNA_process)
      lfs_rna_list <- linear_featureselection_func(LUSC_RNA_process, LUSC_RNA_process_test, association_mat, censor_time1, y, ngenes)
      LUSC_RNA_process <- lfs_rna_list[[1]]
      LUSC_RNA_process_test <- lfs_rna_list[[2]]
      
      #Extracting for altogether AE (early integration)
      mrna_fullAE <- LUSC_RNA_process
      mrna_fullAE_test <- LUSC_RNA_process_test
      mrna_fullAE$barcode <- rna_barcode_train
      mrna_fullAE_test <- as.data.frame(cbind(mrna_fullAE_test, rna_barcode_test))
      colnames(mrna_fullAE_test)[ncol(mrna_fullAE_test)] <- "barcode"
      
      #AUTOENCODERS
      
      ### Denoising Zeros Autoencoder
      if (autoencoder_type == "denoising_zeros") { 
        denoising_zeros_list <- denoising_zeros_func(LUSC_RNA_process, LUSC_RNA_process_test, gene_features)
        LUSC_RNA_auto <- denoising_zeros_list[[1]]
        LUSC_RNA_auto_output <- denoising_zeros_list[[2]] 
        train_RMSE <- denoising_zeros_list[[3]]
        test_RMSE <- denoising_zeros_list[[4]] }
      
      ### Other Autoencoder Type: Denoising Gaussian
      if (autoencoder_type=="denoising_gaussian") {
        denoising_gaussian_list <- denoising_gaussian_func(LUSC_RNA_process, LUSC_RNA_process_test, gene_features)
        LUSC_RNA_auto <- denoising_gaussian_list[[1]]
        LUSC_RNA_auto_output <- denoising_gaussian_list[[2]] 
        train_RMSE <- denoising_gaussian_list[[3]]
        test_RMSE <- denoising_gaussian_list[[4]]
      }
      
      #Making sure variables are numeric
      if (class(LUSC_RNA_process_test[1,1])=="character") {
      LUSC_RNA_process_test <- as.data.frame(sapply(LUSC_RNA_process_test, as.numeric)) }
      LUSC_RNA_auto <- as.data.frame(LUSC_RNA_auto)
      
      #Creating patient ID columns
      LUSC_RNA_auto$barcode <- rna_barcode_train
      LUSC_RNA_auto$barcode <- substring(LUSC_RNA_auto$barcode,1,19)
      LUSC_RNA_auto <- LUSC_RNA_auto %>% mutate(TestTrain="Train")
      multimodal_both <- merge(LUSCclin_lung, LUSC_RNA_auto, by = "barcode") 
      
      # Setting Up Testing Data
      LUSC_RNA_auto_output <- as.data.frame(LUSC_RNA_auto_output)
      LUSC_RNA_test <- LUSC_RNA_auto_output %>% mutate(TestTrain="Test")
      LUSC_RNA_test$barcode <- rna_barcode_test
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
        colnames(multimodal_imputed)[i] <- paste0("Gene Expression", i)
      }
      for (i in 1:((ncol(multimodal_imputed_test)-3))) {
        colnames(multimodal_imputed_test)[i] <- paste0("Gene Expression", i)
      }
      gene_test <- multimodal_imputed_test
      gene_train <- multimodal_imputed
      
      # Saving Train and Test Indices of LUSC/LUAD for Individual Evaluations
      gene_LUSC_indices <- merge(multimodal_imputed, barcode_disease_mapping, by="barcode")
      gene_LUAD_indices_train <- which(gene_LUSC_indices$name=="Lung Adenocarcinoma") 
      gene_LUSC_indices_train <- which(gene_LUSC_indices$name=="Lung Squamous Cell Carcinoma")
      gene_LUSC_indices_test <- merge(multimodal_imputed_test, barcode_disease_mapping, by="barcode")
      gene_LUAD_indices_test <- which(gene_LUSC_indices_test$name=="Lung Adenocarcinoma") 
      gene_LUSC_indices_test <- which(gene_LUSC_indices_test$name=="Lung Squamous Cell Carcinoma")
      
      # Getting rid of barcode variable and making sure data is numeric
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
      
      ### Autoencoder Performance
      AE_RMSE[2] <- test_RMSE
      AE_RMSE_train[2] <- train_RMSE
      
      
      # MIRNA MODALITY
      
      # Preprocessing
      mirna_list <- preprocessing_func(LUSC_miRNA, "h")
      LUSC_miRNA_process <- mirna_list[[1]]
      LUSC_miRNA_process_test <- mirna_list[[2]]
      vital_vector <- mirna_list[[3]]
      censor_time1 <- mirna_list[[4]]
      mirna_barcode_train <- mirna_list[[5]]
      mirna_barcode_test <- mirna_list[[6]]
      
      #Linear Feature Selection
      y <- (as.numeric(factor(vital_vector))-1)
      association_mat <- as.matrix(LUSC_miRNA_process)
      lfs_mirna_list <- linear_featureselection_func(LUSC_miRNA_process, LUSC_miRNA_process_test, association_mat, censor_time1, y, n_mirna)
      LUSC_miRNA_process <- lfs_mirna_list[[1]]
      LUSC_miRNA_process_test <- lfs_mirna_list[[2]]
      
      #Extracting for altogether AE
      mirna_fullAE <- LUSC_miRNA_process
      mirna_fullAE_test <- LUSC_miRNA_process_test
      mirna_fullAE$barcode <- mirna_barcode_train
      mirna_fullAE_test <- as.data.frame(cbind(mirna_fullAE_test, mirna_barcode_test))
      colnames(mirna_fullAE_test)[ncol(mirna_fullAE_test)] <- "barcode"
      
      ### Denoising Autoencoder FOR miRNA
      
      ### Denoising Zeros Autoencoder
      if (autoencoder_type == "denoising_zeros") { 
        denoising_zeros_list <- denoising_zeros_func(LUSC_miRNA_process, LUSC_miRNA_process_test, mirna_features) 
        LUSC_miRNA_auto <- denoising_zeros_list[[1]]
        LUSC_miRNA_auto_output <- denoising_zeros_list[[2]] 
        train_RMSE <- denoising_zeros_list[[3]]
        test_RMSE <- denoising_zeros_list[[4]]
        }
      
      ### Other Autoencoder Type: Denoising (Gaussian), 
      if (autoencoder_type=="denoising_gaussian") {
        denoising_gaussian_list <- denoising_gaussian_func(LUSC_miRNA_process, LUSC_miRNA_process_test, mirna_features)
        LUSC_miRNA_auto <- denoising_gaussian_list[[1]]
        LUSC_miRNA_auto_output <- denoising_gaussian_list[[2]] 
        train_RMSE <- denoising_gaussian_list[[3]]
        test_RMSE <- denoising_gaussian_list[[4]]
        }
      
      #Turning into dataframes and creating patient ID columns
      LUSC_miRNA_auto <- as.data.frame(LUSC_miRNA_auto)
      LUSC_miRNA_auto$barcode <- mirna_barcode_train
      LUSC_miRNA_auto$barcode <- substring(LUSC_miRNA_auto$barcode,1,19)
      LUSC_miRNA_auto <- LUSC_miRNA_auto %>% mutate(TestTrain="Train")
      multimodal_both <- merge(LUSCclin_lung, LUSC_miRNA_auto, by = "barcode") 
      
      # Testing Data
      LUSC_miRNA_auto_output <- as.data.frame(LUSC_miRNA_auto_output)
      LUSC_miRNA_test <- LUSC_miRNA_auto_output %>% mutate(TestTrain="Test")
      LUSC_miRNA_test$barcode <- mirna_barcode_test
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
        colnames(multimodal_imputed)[i] <- paste0("miRNA", i) }
      for (i in 1:((ncol(multimodal_imputed_test)-3))) {
        colnames(multimodal_imputed_test)[i] <- paste0("miRNA", i) }
      
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
      
      # Autoencoder Performance
      AE_RMSE[3] <- test_RMSE
      AE_RMSE_train[3] <- train_RMSE

      
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
      
      #Linear Feature Selection
      y <- (as.numeric(factor(vital_vector))-1)
      association_mat <- as.matrix(LUSC_meth_process)
      lfs_meth_list <- linear_featureselection_func(LUSC_meth_process, LUSC_meth_process_test, association_mat, censor_time1, y, nprobes)
      LUSC_meth_process <- lfs_meth_list[[1]]
      LUSC_meth_process_test <- lfs_meth_list[[2]]
      
      #Extracting for altogether AE
      meth_fullAE <- LUSC_meth_process
      meth_fullAE_test <- LUSC_meth_process_test
      meth_fullAE <- as.data.frame(cbind(meth_fullAE, meth_barcode_train))
      meth_fullAE_test <- as.data.frame(cbind(meth_fullAE_test, meth_barcode_test))
      colnames(meth_fullAE)[ncol(meth_fullAE)] <- "barcode"
      colnames(meth_fullAE_test)[ncol(meth_fullAE_test)] <- "barcode"
      
      ### Denoising Zeros Autoencoder
      if (autoencoder_type == "denoising_zeros") { 
        denoising_zeros_list <- denoising_zeros_func(LUSC_meth_process, LUSC_meth_process_test, meth_features) 
        LUSC_meth_auto <- denoising_zeros_list[[1]]
        LUSC_meth_auto_output <- denoising_zeros_list[[2]]
        train_RMSE <- denoising_zeros_list[[3]]
        test_RMSE <- denoising_zeros_list[[4]]
        }
      
      ### Other Autoencoder Type: Denoising Gaussian
      if (autoencoder_type=="denoising_gaussian") {
        denoising_gaussian_list <- denoising_gaussian_func(LUSC_meth_process, LUSC_meth_process_test, meth_features) 
        LUSC_meth_auto <- denoising_gaussian_list[[1]]
        LUSC_meth_auto_output <- denoising_gaussian_list[[2]] 
        train_RMSE <- denoising_gaussian_list[[3]]
        test_RMSE <- denoising_gaussian_list[[4]] }
      
      ### Putting Together the Dataframe
      LUSC_meth_auto <- as.data.frame(LUSC_meth_auto)
      LUSC_meth_auto$barcode <- meth_barcode_train
      LUSC_meth_auto <- LUSC_meth_auto %>% mutate(TestTrain="Train")
      multimodal_both <- merge(LUSCclin_lung, LUSC_meth_auto, by = "barcode")
      
      # Testing Data
      LUSC_meth_auto_output <- as.data.frame(LUSC_meth_auto_output)
      LUSC_meth_test <- LUSC_meth_auto_output %>% mutate(TestTrain="Test")
      LUSC_meth_test$barcode <- meth_barcode_test
      multimodal_both2 <- merge(LUSCclin_lung, LUSC_meth_test, by = "barcode")
      multimodal_both <- rbind(multimodal_both, multimodal_both2)
      
      # Removing unused variables
      multimodal_prediction <- multimodal_both %>%
        dplyr::select(-c(patient, treatments, days_to_last_follow_up, days_to_death, shortLetterCode))
      
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
      
      # Autoencoder Performance
      AE_RMSE[4] <- test_RMSE
      AE_RMSE_train[4] <- train_RMSE
      
      
      # Clinical Modality
      
      ### Splitting into training and testing
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
      multimodal_prediction <- multimodal_prediction %>% dplyr::select(-c(censor_time1, vital_status,TestTrain, barcode))
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
        multimodal_imputed[ , i][is.na(multimodal_imputed[ , i])] <- median(multimodal_imputed[ , i], na.rm=TRUE)
      }
      
      #Imputation loop - testing data
      for (i in 1:(stop_iteration-1)) {
        multimodal_imputed_test[ , i][is.na(multimodal_imputed_test[ , i])] <- median(multimodal_imputed_test[ , i], na.rm=TRUE)
      }
      
      #Imputation loop - overall PCA data
      for (i in 1:(ncol(multimodal_prediction)-4)) {
        multimodal_prediction[ , i][is.na(multimodal_prediction[ , i])] <- median(multimodal_prediction[ , i], na.rm=TRUE)
      }
      
      # Taking dataframe for multimodal section
      clin_test <- multimodal_imputed_test
      clin_train <- multimodal_imputed
      clin_LUSC_indices <- merge(multimodal_imputed, barcode_disease_mapping, by="barcode")
      clin_LUAD_indices_train <- which(clin_LUSC_indices$name.y=="Lung Adenocarcinoma") 
      clin_LUSC_indices_train <- which(clin_LUSC_indices$name.y=="Lung Squamous Cell Carcinoma")
      clin_LUSC_indices_test <- merge(multimodal_imputed_test, barcode_disease_mapping, by="barcode")
      clin_LUAD_indices_test <- which(clin_LUSC_indices_test$name.y=="Lung Adenocarcinoma") 
      clin_LUSC_indices_test <- which(clin_LUSC_indices_test$name.y=="Lung Squamous Cell Carcinoma")
      multimodal_imputed_test <- multimodal_imputed_test %>% dplyr::select(c(-barcode))
      multimodal_imputed <- multimodal_imputed %>% dplyr::select(c(-barcode))
      multimodal_imputed_test <- as.data.frame(sapply(multimodal_imputed_test, as.numeric))
      multimodal_imputed <- as.data.frame(sapply(multimodal_imputed, as.numeric))
      
      # Name Variable Means Nothing if Training on Only One DataType
      if (only_train_on_LUAD=="yes" | only_train_on_LUSC=="yes") {
        multimodal_imputed <- multimodal_imputed %>% dplyr::select(-c(name))
        multimodal_imputed_test <- multimodal_imputed_test %>% dplyr::select(-c(name))
        }

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
      
      
      ### Multimodal Analysis
      
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
      
      # Name Variable Means Nothing if Training on Only One DataType
      if (only_train_on_LUAD=="yes" | only_train_on_LUSC=="yes") {
        multimodal_train <- multimodal_train %>% dplyr::select(-c(name))
        multimodal_test <- multimodal_test %>% dplyr::select(-c(name))
      }
      
      # Making Numeric
      multimodal_test <- as.data.frame(sapply(multimodal_test, as.numeric))
      multimodal_train <- as.data.frame(sapply(multimodal_train, as.numeric))

      ### Elastic Net
      multimodal_test_LUSC_multi <- multimodal_test[c(multimodal_LUSC_indices_test), ]
      multimodal_test_LUAD_multi <- multimodal_test[c(multimodal_LUAD_indices_test), ]
      elasticnet_cindex_list <- elastic_net_func(multimodal_train, multimodal_test, multimodal_test_LUSC_multi, multimodal_test_LUAD_multi)
      Cindex_Test <- elasticnet_cindex_list[[1]]
      Cindex_LUSC <- elasticnet_cindex_list[[2]]
      Cindex_LUAD <- elasticnet_cindex_list[[3]]
      est.coef <- elasticnet_cindex_list[[4]]
      
      # Elastic Net Multimodal Feature Importance Saving
      coef_mat <- as.matrix(est.coef)
      coef_mat <- as.data.frame(coef_mat)
      feat_importance <- cbind(coef_mat, colnames(multimodal_train)[1:(ncol(multimodal_train)-2)])
      Run_Number <- rep(paste0("Run", cv), nrow(feat_importance))
      feat_importance <- cbind(feat_importance, Run_Number)
      colnames(feat_importance) <- c("Coefficients", "Variable", "Run_Number")
      feat_importance_list[[cv]] <- feat_importance
      
      ### Results 
      Data_Types <- "Multimodal (Late)"
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
      
      # Saving feaure importance data and data for clustering analysis
      elastic_variables[[cv]] <- feat_importance$Variable[which(feat_importance$Coefficients!=0)]
      Y <- as.data.frame(rbind(multimodal_train, multimodal_test))
      barcode_vec <- c(train_code, test_code)
      Clustering_Data[[cv]] <- as.data.frame(cbind(Y, barcode_vec))
      
      ### Looping through all combinations of modalities (LATE INTEGRATION)
      
      if (all_modalities_combinations_late=="yes") {
      
      # Setting Indices
      censor_index <- which(colnames(multimodal_train)=="censor_time")
      vital_index <- which(colnames(multimodal_train)=="vital_status")
      outcomes_index <- c(censor_index, vital_index)
      meth_ind <- c(111:160)
      gene_ind <- c(61:110)
      mirna_ind <- c(1:30)
      lnc_ind <- c(31:60)
      clin_ind <- c(161:171)
      modality_list <- list(meth_ind, gene_ind, mirna_ind, lnc_ind, clin_ind, outcomes_index)
      
      #Modality Combinations 
      #(Blocked Out Clinical Data Combinations Since There is No Early/Late Integration Difference for One Modality + Clinical)
      #Feel Free to Uncomment if Desired
      data_combos <- c("Meth-Gene", "Meth-LNC", "Meth-miRNA", #"Meth-Clin", 
                       "Gene-LNC", "Gene-miRNA", 
                       #"Gene-Clin", 
                       "LNC-miRNA", 
                       #"LNC-Clin", "miRNA-Clin",
                       "Meth-Gene-LNC", "Meth-Gene-Clin", "Meth-Gene-miRNA", "Meth-Clin-LNC", "Meth-Clin-miRNA", 
                       "Meth-LNC-miRNA", "Gene-miRNA-LNC", "Gene-miRNA-Clin", "Gene-LNC-Clin",
                       "miRNA-LNC-Clin", "Meth-Gene-LNC-Clin","Meth-Gene-miRNA-Clin","Meth-Gene-LNC-miRNA",
                       "Meth-LNC-Clin-miRNA","Gene-LNC-Clin-miRNA")
      combo_list <- list(c(meth_ind, gene_ind, outcomes_index), c(meth_ind, lnc_ind, outcomes_index),c(meth_ind, mirna_ind, outcomes_index),
                         #c(meth_ind, clin_ind, outcomes_index),
                         c(lnc_ind, gene_ind, outcomes_index),c(mirna_ind, gene_ind, outcomes_index),#c(clin_ind, gene_ind, outcomes_index),
                         c(lnc_ind, mirna_ind, outcomes_index),
                         #c(lnc_ind, clin_ind, outcomes_index),c(mirna_ind, clin_ind, outcomes_index),
                         #Three-Omics
                         c(meth_ind, gene_ind, lnc_ind, outcomes_index),c(meth_ind, gene_ind, clin_ind, outcomes_index),c(meth_ind, gene_ind, mirna_ind, outcomes_index),
                         c(meth_ind, clin_ind, lnc_ind, outcomes_index),c(meth_ind, mirna_ind, clin_ind, outcomes_index),c(meth_ind, mirna_ind, lnc_ind, outcomes_index),c(gene_ind, mirna_ind, lnc_ind, outcomes_index),
                         c(gene_ind, mirna_ind, clin_ind, outcomes_index),c(gene_ind, clin_ind, lnc_ind, outcomes_index),c(mirna_ind, clin_ind, lnc_ind, outcomes_index),
                         #Four-Omics
                         c(meth_ind, gene_ind, lnc_ind, clin_ind,outcomes_index),c(meth_ind, gene_ind, mirna_ind, clin_ind,outcomes_index),
                         c(meth_ind, gene_ind, lnc_ind, mirna_ind,outcomes_index),
                         c(meth_ind, mirna_ind, lnc_ind, clin_ind,outcomes_index),
                         c(gene_ind, mirna_ind, lnc_ind, clin_ind,outcomes_index))
      mod_combo_results <- list()
      
      for (b in 1:length(data_combos)) {
        combinations_ind <- combo_list[[b]]
        data_label <- data_combos[b]
        
        #Getting rid of clinical features
        multimodal_combo_train <- multimodal_train[, combinations_ind]
        multimodal_combo_test <- multimodal_test[, combinations_ind]
      
        ### Elastic Net
        multimodal_combo_test_LUAD <- multimodal_combo_test[c(multimodal_LUAD_indices_test), ]
        multimodal_combo_test_LUSC <- multimodal_combo_test[c(multimodal_LUSC_indices_test), ]
        elasticnet_cindex_list <- elastic_net_func(multimodal_combo_train, multimodal_combo_test, multimodal_combo_test_LUSC, multimodal_combo_test_LUAD)
        Cindex_Test <- elasticnet_cindex_list[[1]]
        Cindex_LUSC <- elasticnet_cindex_list[[2]]
        Cindex_LUAD <- elasticnet_cindex_list[[3]]
        
        # Saving Results 
        Data_Types <- data_label
        results.df <- as.data.frame(cbind(Models, Cindex_Test, Dimension_Reduction, Data_Types,
                                          Cindex_LUSC, Cindex_LUAD))
        results.df$Cindex_Test <- as.numeric(results.df$Cindex_Test)
        results.df$Cindex_Test <- round(results.df$Cindex_Test, 3)
        mod_combo_results[[b]] <- results.df
      }
      all_results_df <- bind_rows(mod_combo_results)
    }
      
      ### Testing All Combinations of Data (Early Integration)
      
      if (one_AE_early=="yes") {
        #Training Data
        df3 <- merge(lncrna_fullAE, mirna_fullAE, by="barcode")
        df3 <- merge(df3, mrna_fullAE, by="barcode")
        fullAE <- merge(df3, meth_fullAE, by="barcode")
        
        #Testing Data
        df4 <- merge(lncrna_fullAE_test, mirna_fullAE_test, by="barcode")
        df4 <- merge(df4, mrna_fullAE_test, by="barcode")
        fullAE_test <- merge(df4, meth_fullAE_test, by="barcode")
        
        #Save vital status and outcome into a dataframe
        fullAE_outcome <- fullAE %>% dplyr::select(c(barcode))
        fullAE_outcome_test <- fullAE_test %>% dplyr::select(c(barcode))
        
        #Correct Indices Train
        full_LUSC_indices <- merge(fullAE, barcode_disease_mapping, by="barcode")
        full_LUAD_indices_train <- which(full_LUSC_indices$name=="Lung Adenocarcinoma") 
        full_LUSC_indices_train <- which(full_LUSC_indices$name=="Lung Squamous Cell Carcinoma")
        
        #Correct Indices Test
        full_LUSC_indices_test <- merge(fullAE_test, barcode_disease_mapping, by="barcode")
        full_LUAD_indices_test <- which(full_LUSC_indices_test$name=="Lung Adenocarcinoma") 
        full_LUSC_indices_test <- which(full_LUSC_indices_test$name=="Lung Squamous Cell Carcinoma")
        
        #Removing Barcodes
        fullAE <- fullAE %>% dplyr::select(-c(barcode))
        fullAE_test <- fullAE_test %>% dplyr::select(-c(barcode))
        fullAE <- as.data.frame(sapply(fullAE, as.numeric))
        fullAE_test <- as.data.frame(sapply(fullAE_test, as.numeric))
        
        #Scaling Data
        fullAE_test <- as.matrix(fullAE_test)
        fullAE_test <- scale(fullAE_test)
        fullAE <- scale(fullAE)
        
        #Getting all combinations
        if (one_AE_early_allcombos == "yes") {
        lnc_ind_full <- c(1:300)
        mirna_ind_full  <- c(301:600)
        gene_ind_full  <- c(601:1100)
        meth_ind_full  <- c(1101:1600)
        modality_list <- list(lnc_ind_full, mirna_ind_full, meth_ind_full, gene_ind_full)
        #All Combinations of Modalities and Feature Numbers
        data_combos_full <- c("Meth-Clin", "Gene-Clin", "LNC-Clin", "miRNA-Clin",
          "Meth-Gene","Meth-Gene-Clin", "Meth-LNC", "Meth-LNC-Clin","Meth-miRNA", "Meth-miRNA-Clin", "Gene-LNC", "Gene-LNC-Clin",
                              "Gene-miRNA", "Gene-miRNA-Clin", "LNC-miRNA", "LNC-miRNA-Clin",
                              "Meth-Gene-LNC", "Meth-Gene-LNC-Clin", "Meth-Gene-miRNA","Meth-Gene-miRNA-Clin", 
                              "Meth-LNC-miRNA", "Meth-LNC-miRNA-Clin", "Gene-miRNA-LNC",  "Gene-miRNA-LNC-Clin",
                              "Meth-Gene-LNC-miRNA", "Multimodal (Early)")
        feature_length <- c(50, 50, 30, 30, 100, 100, 80, 80, 80, 80, 80, 80, 80, 80, 60,60,
                            130, 130, 130, 130, 110, 110, 110, 110, 160)
        combo_list_full <- list(c(meth_ind_full), c(gene_ind_full), c(lnc_ind_full), c(mirna_ind_full),
                                c(meth_ind_full, gene_ind_full), c(meth_ind_full, gene_ind_full),c(meth_ind_full, lnc_ind_full),c(meth_ind_full, lnc_ind_full),
                                c(meth_ind_full, mirna_ind_full), c(meth_ind_full, mirna_ind_full),
                                c(lnc_ind_full, gene_ind_full),c(lnc_ind_full, gene_ind_full),c(mirna_ind_full, gene_ind_full),c(mirna_ind_full, gene_ind_full),
                                c(lnc_ind_full, mirna_ind_full),c(lnc_ind_full, mirna_ind_full),
                                #Three-Omics
                                c(meth_ind_full, gene_ind_full, lnc_ind_full),c(meth_ind_full, gene_ind_full, lnc_ind_full),c(meth_ind_full, gene_ind_full, mirna_ind_full),c(meth_ind_full, gene_ind_full, mirna_ind_full),
                                c(meth_ind_full, mirna_ind_full, lnc_ind_full),c(meth_ind_full, mirna_ind_full, lnc_ind_full),c(gene_ind_full, mirna_ind_full, lnc_ind_full),c(gene_ind_full, mirna_ind_full, lnc_ind_full),
                                #Four-Omics
                                c(meth_ind_full, gene_ind_full, lnc_ind_full, mirna_ind_full),
                                # Multimodal
                                c(meth_ind_full, gene_ind_full, lnc_ind_full, mirna_ind_full))
        
        #Looping through combinations
        for (c in 1:length(data_combos_full)) {
          combinations_ind_full <- combo_list_full[[c]]
          data_label_full <- data_combos_full[c]
          
          fullAE_combo <- fullAE[ ,combinations_ind_full]
          fullAEtest_combo <- fullAE_test[ ,combinations_ind_full]
        
          ### Denoising Zeros Autoencoder
          if (autoencoder_type == "denoising_zeros") { 
            denoising_zeros_func(fullAE_combo, fullAEtest_combo, fullAE_features) 
            LUSC_multimodal_auto <- denoising_zeros_list[[1]]
            LUSC_multimodal_auto_output <- denoising_zeros_list[[2]] 
            train_RMSE <- denoising_zeros_list[[3]]
            test_RMSE <- denoising_zeros_list[[4]] }
          
          ### Other Autoencoder Type: Denoising Gaussian
          if (autoencoder_type=="denoising_gaussian") {
            denoising_gaussian_list <- denoising_gaussian_func(fullAE_combo, fullAEtest_combo, fullAE_features) 
            LUSC_multimodal_auto <- denoising_gaussian_list[[1]]
            LUSC_multimodal_auto_output <- denoising_gaussian_list[[2]] 
            train_RMSE <- denoising_gaussian_list[[3]]
            test_RMSE <- denoising_gaussian_list[[4]] }  
        
        ### Saving Autoencoder Performance
        AE_RMSE[5] <- test_RMSE
        AE_RMSE_train[5] <- train_RMSE
        
        # Incorporating outcome data
        LUSC_AE <- as.data.frame(cbind(LUSC_multimodal_auto, fullAE_outcome))
        LUSC_AE_test <- as.data.frame(cbind(LUSC_multimodal_auto_output,fullAE_outcome_test))
        colnames(LUSC_AE)[ncol(LUSC_AE)] <- "barcode"
        colnames(LUSC_AE_test)[ncol(LUSC_AE_test)] <- "barcode"
        
        # Incorporating clinical data
        if (grepl("Clin", data_label_full) | data_label_full == "Multimodal") {
        LUSC_fullAE <- merge(LUSC_AE, clin_train, by="barcode")
        LUSC_fullAE_test <- merge(LUSC_AE_test, clin_test, by="barcode") }
        if ("barcode" %in% colnames(LUSC_fullAE)) { LUSC_fullAE_barcode <- LUSC_fullAE$barcode }
        if ("barcode" %in% colnames(LUSC_fullAE_test)) { LUSC_fullAE_barcode_test <- LUSC_fullAE_test$barcode }
        
        #Removing Barcodes
        if ("barcode" %in% colnames(LUSC_fullAE)) { LUSC_fullAE <- LUSC_fullAE %>% dplyr::select(-c(barcode)) }
        if ("barcode" %in% colnames(LUSC_fullAE_test)) { LUSC_fullAE_test <- LUSC_fullAE_test %>% dplyr::select(-c(barcode)) }
        LUSC_fullAE <- as.data.frame(sapply(LUSC_fullAE, as.numeric))
        LUSC_fullAE_test <- as.data.frame(sapply(LUSC_fullAE_test, as.numeric))
        
        #Elastic Net
        multimodal_test_LUSC <- LUSC_fullAE_test[c(full_LUSC_indices_test), ]
        multimodal_test_LUAD <- LUSC_fullAE_test[c(full_LUAD_indices_test), ]
        elasticnet_cindex_list <- elastic_net_func(LUSC_fullAE, LUSC_fullAE_test, multimodal_test_LUSC, multimodal_test_LUAD)
        Cindex_Test <- elasticnet_cindex_list[[1]]
        Cindex_LUSC <- elasticnet_cindex_list[[2]]
        Cindex_LUAD <- elasticnet_cindex_list[[3]]
        
        # Results
        ### Results 
        Dimension_Reduction <- "One AE"
        Data_Types <- data_label_full
        results.df <- as.data.frame(cbind(Models, Cindex_Test, Dimension_Reduction, Data_Types,
                                          Cindex_LUSC, Cindex_LUAD))
        results.df$Cindex_Test <- as.numeric(results.df$Cindex_Test)
        results.df$Cindex_Test <- round(results.df$Cindex_Test, 3)
        results.df_fullAE <- results.df
        results.df_fullAE_list[[c]] <- results.df_fullAE
        }
        results.df_fullAE_allcombos <- bind_rows(results.df_fullAE_list)
        }
      }
      
      #Saving Autoencoder Results (RMSE for Training and Testing)
      Modalities <- c("LNC_RNA", "mRNA", "miRNA", "Methylation")
        AE_Info <- as.data.frame(cbind(AE_RMSE, AE_RMSE_train)) 
        colnames(AE_Info) <- c("Test_RMSE", "Train_RMSE")
      if (one_AE_early=="yes") { Modalities <- c(Modalities, "one_AE_early") }
      AE_Info <- as.data.frame(cbind(AE_Info, Modalities))
      autoencoder_list[[cv]] <- AE_Info
      
      # Joining Overall Results of Different Modalities 
      df1 <- rbind(results.df_gene, results.df_meth)
      df1 <- rbind(df1, results.df_lncs)
      df1 <- rbind(df1, results.df_mirna)
      df1 <- rbind(df1, results.df_clin)
      
      #Adding in all Combination Data and Saving Data to a List
      if (all_modalities_combinations_late=="yes") { df1 <- rbind(df1, all_results_df) }
      if (one_AE_early_allcombos=="yes") { df1 <- rbind(df1,results.df_fullAE_allcombos) }
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
    
    #Saving Autoencoder Results
    autoencoder_results <- bind_rows(autoencoder_list)
    write.csv(autoencoder_results, file=paste0("~/desktop/Multimodal/MultimodalPipeline_Results_Post/", runinfo, "/AE_Information.csv"))
    
    #Saving Feature Importance Results
    feat_importance_elastic <- bind_rows(feat_importance_list)
    write.csv(feat_importance_elastic, file=paste0("~/desktop/Multimodal/MultimodalPipeline_Results_Post/", runinfo, "/Elastic_Feat_Importance.csv"))
    
    #Saving PCA Data in List Form for Future Visualizations
    save(total_PCA_list, file=paste0("~/desktop/Multimodal/MultimodalPipeline_Results_Post/", 
                                     runinfo, "/PCA_Visualization_List.Rdata"))
    
    # Saving Multimodal Data for Clustering and Elastic Feature Importance Data
    save(Clustering_Data, file=paste0("~/desktop/Multimodal/MultimodalPipeline_Results_Post/", 
                                      runinfo, "/clustering_data.Rdata"))
    save(elastic_variables, file=paste0("~/desktop/Multimodal/MultimodalPipeline_Results_Post/", 
                                        runinfo, "/elastic_variable.Rdata"))
  }
}
