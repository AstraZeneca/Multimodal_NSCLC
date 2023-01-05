# Autoencoder-Based Multimodal Prediction of Non-Small Cell Lung Cancer Survival

This study presents a novel, autoencoder-based pipeline for both adenocarcinoma (LUAD) and squamous cell carcinoma (LUSC) patients that aims to predict survival more effectively and identify NSCLC survival subtypes. We implemented a combined linear feature selection and denoising autoencoder approach for the compression of multiple modalities of data, including microRNA (miRNA), mRNA, DNA methylation, long non-coding RNAs (lncRNAs) and clinical data. Survival performance, as measured by the concordance index (C-index), was compared across data dimensionality reduction techniques, modality combinations, time of data integration and training data types. A regularized Cox proportional hazards model, Elastic Net, was chosen to make survival predictions using a 5-fold cross validation scheme on 732 NSCLC (408 LUAD; 324 LUSC) patients from the Cancer Genome Atlas (TCGA) dataset.

# Implementation Instructions

1. Download the zipped 'Multimodal' folder here: [Multimodal.zip](https://github.com/AstraZeneca/Multimodal_NSCLC/blob/master/Multimodal.zip?raw=true)
2. Unzip the Multimodal folder and add it to your desktop.
2. If you have not downloaded R, go to https://www.r-project.org/ and click the “Download R” link under the “Getting Started” header. Next, go to https://rstudio.com/products/rstudio/download/, choose the “RStudio Desktop” Option and click “Download.”
3. Open and run the 'DownloadingData.Rmd' script above in R or R Studio, which will download mRNA, miRNA and DNA Methylation data directly from TCGA and save it into the Multimodal folder. 
4. Run the 'PreProcessing.R' script to process all of the data for all five modalities, including extracting lncRNA data from the mRNA file, removing biological transcripts with too many missing values/zeros (see paper) and taking only the 732 eligible patients. Again, the processed files will be saved in the Multimodal folder.
5. Run the files below (see descriptions) to get results and then run the 'PostProcessing.Rmd' script to analyze those results and create summary figures.

NOTE: Updates to the TCGA database and how the NSCLC data are stored over time can create errors in the data downloading script so feel free to reach out to jgellen4@gmail.com with any questions about TCGA updates.

# Pipeline Script Descriptions

**1. AllModality_Comparison.R**

- This is the main script that implements the entire pipeline described in the paper. Starting with the preprocessed files, the data is split into five folds for 5-fold cross validation. The entire pipeline is run for each fold and all of the results are saved in the Multimodal folder. For each train/test split, the training data for each modality is scaled, run through linear feature selection (chosen features are used for the testing data), run through a denoising autoencoder and then used to predict survival using the Elastic Net model. The denoising autoencoder type and hyperparameters can be adjusted at the top of the script. Predictions are made for each single modality and for every modality together. Later, the script will also test all combinations of modalities using both late integration and early integration (see paper for description). Since this script tests all combinations, it will likely take a few hours to run on a CPU. Lastly, this script gives the option to only train on LUSC data or only train on LUAD data to see how the results might change.

**2. LFS.R**

- The 'LFS.R' script works similarly to the 'AllModality_Comparison.R' script, but only implements linear feature selection (with no autoencoder) in order to compare our combined feature selection and denoising autoencoder approach to a pure feature selection approach.

**3. PCA.R**

- This script uses PCA as its dimensionality reduction technique instead of linear feature selection and/or denoising autoencoders.

**4. OnlyAutoencoder.R**

- This script uses only autoencoders to compare our combined feature selection and denoising autoencoder approach to a pure denoising autoencoder approach to see if linear feature selection is even necessary.

**5. PostProcessing.Rmd**

- This post-processing script takes all of the results from above and analyzes the data and produces figures. Figures will be created comparing multimodal and 2-omics approaches to single modalities by survival performance, comparing data dimensionality reduction techniques performance, examining the differences between late and early integration of multimodal data and visualizing how training on one or both data types (LUAD and LUSC) can affect performance. There is the option to save these figures directly into the Multimodal folder if desired. This script will only work fully if all of the above scripts are run, but most figures can still be made with just the 'AllModalities_Comparison.R' output.

# Survival Subtype Script Descriptions

**1. KMeans_Clustering.R**

- This script takes the multimodal feature space from a chosen pipeline run, uses the Silhouette Method to find the optimal number of clusters and then performs K-means clustering with the determined cluster number. A Kaplan-Meier curve comparing the survival subtypes is generated and the barcodes of each group are saved in order to perform differential expression analyses below. Thus, it is important to note that you must run the 'AllModality_Comparison.R' script above in order to get the multimodal feature space used for K Means Clustering. You also must run this script before running any of the DEA scripts below for each modality.

**2. DEA Scripts: DEA_mRNA.Rmd, DEA_miRNA.Rmd, DEA_Methylation.Rmd and DEA_LNCRNAS.Rmd**

- Each of these scripts performs differential expression analysis using the groups generated from the K-means clustering above. For the miRNA, mRNA and lncRNA scripts, count data is downloaded, differential expression analysis is carried out and volcano plots are generated showing transcripts that are differentially expressed. In the multimodal folder, .csv files are saved with the list of 

# Important File Descriptions

**1. LUAD_indices_sample_5fold.Rdata and LUSC_indices_sample_5fold.Rdata**

- NOTE: These files are contained Within the Multimodal folder downloaded above.
- Using the set.seed function in R, these are the five randomly generated cross validation fold groups used in this study for both LUAD and LUSC. The script will run without them, but this is just to make sure that differences in R versions don't lead to differences in these fold groups.

**2. barcode_disease_mapping.Rdata**

- List of all 732 barcodes of patients utilized in this study.
