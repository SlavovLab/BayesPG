# BayesPG

[Website](https://scp.slavovlab.net/Khan_Elcheikhali_et_al_2024) &nbsp; | &nbsp; [Preprint](https://www.biorxiv.org/content/...) &nbsp; | &nbsp; [Raw data](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000096034) &nbsp; | &nbsp; [Processed data](https://drive.google.com/drive/folders/1oMaW74bjDPPPjyPTqocPD_SjMnHy7UeM?usp=sharing)

*Currently a work in progress. Please bear with us while we upload all code to repositories in the coming days.*

This repository hosts code for our workflow to systematically infer post transcriptional regulation as reflected in differences between relative protein and mRNA levels. To account for and distinguish between sources of variation in our data, we developed BayesPG, a hierarchical Bayesian model. BayesPG jointly models cell type level mRNA and protein estimates obtained from single-cell data. We employed BayesPG to quantitatively explore post transcriptional regulation in primary human testis cells. 

### Workflow Overview:

The workflow is divided into three main stages:

1. **Data Generation and Preprocessing** (`001-MSData`)  
   This stage includes scripts used for acquiring and processing pSCoPE and plexDIA data. This includes various helper scripts for data acquisition or database searches, such as inclusion list generation and spectral library IM prediction.

2. **Single-Cell Matrix Processing** (`002-singleCell_data`)  
   Scripts for generating, processing and working with single-cell matrices (protein and mRNA). This stage includes:
   - Normalization and batch correction.
   - Feature selection for dataset alignment using correlation vector analysis.
   - Alignment using **LIGER** and validation of alignment results.
   - Analysis of post-translationally modified proteins within clusters and across single cells (**Figure 6** in manuscript).

3. **Running BayesPG on Cluster-Level Data** (`003-BayesPG`)  
   This section contains code forked from the main BayesPG repository, as well as a small script to reproduce **Figure 4** from the manuscript.
