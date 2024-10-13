# singleCell_data


This sub folder contains scripts relating to processing all single cell matrices (mRNA and protein) in different forms, selecting the feature space for LIGER alignment, evaluating alignment metrics, making figures etc. 

All file paths (read and write) listed in these script correspond to the processed data folder directory structure that can be found [here](https://drive.google.com/drive/folders/1oMaW74bjDPPPjyPTqocPD_SjMnHy7UeM?usp=sharing).

Below are brief descriptions for each script, enumerated in (broadly) the order with which they are to be used, following in the steps from processing MS data into single cell matrices to LIGER alignments to the PTM analysis. 

The scripts themselves have annotations in the headers that provide more details.   

1. **bayesPG_SCoPE2_singleCellProcessing.R**: A script for processing DDA/SCoPE2 style acquisitions (datasets/npops 1- 4). This is a faster implementation of the standard [SCoPE2 pipeline] (https://github.com/SlavovLab/SCoPE2) for processing Maxquant evidence files (updated by [DART-ID] (https://github.com/SlavovLab/DART-ID)) into single cell matrices. 
	1. **functions_parameters_speedUpdate**: source functions that accompany this pipeline. 
2. **bayesPG_plexDIA_singleCellProcessing.R**: A script to process [plexDIA] (https://github.com/SlavovLab/plexDIA) data acquired on the TIMS TOF SCP. This script is implemented using the [QuantQC] (https://github.com/SlavovLab/QuantQC) package.
3. **bayesPG_preProcessMatrices_featureSelect.R**: A rather large script that combines multiple tasks which rely on similarly large objects that are loaded into the workspace. Processes single cell mRNA-Seq datasets (normalization, batch correction) and protein datasets (Batch correction). The feature space for dataset alignment is then selected via correlation vector analysis. The eventual outputs of the alignment are then evaluated (extended data figure 2B). Finally, code for the ANOVA analysis that constitutes extended data figure 1 is added at the end of the script.
4. **bayesPG_LigerAlignment_extendFig2.R**: The script implements the mRNA and protein dataset alignments using [LIGER] (https://github.com/welch-lab/liger). Extended data figure 2A is produced. Post alignment, the transfer of joint cluster identities (for all datasets, both mRNA and protein) is carried out.
5. **bayesPG_fig6_PTM.preprocess.R**: Script to process the variable modification Maxquant evidence files (updated by [DART-ID] (https://github.com/SlavovLab/DART-ID)). Filtration steps to compute a local PTM specific FDR (and filter using) as well as filtering for PTMs based on observations across cells is carried out. Finally, figure 6A is made.
6. **bayesPG_fig6_PTM.Analysis.R**: Script to carry out the PTM Analysis across single cells. Computes and tests significance of celltype level differences in the relative abundances of phosphopeptides and kinases,  finally displaying results as a clustered dot plot (figure 6b) across all 6 cell types. An illustrative scatter plot (figure 6c) is produced for an example phosphopeptide<>kinase pair across the 3 spermatogonial stem cell types (Spermatogonia -> Spermatocytes -> Spermatids. The differences in correlations between phosphopeptide<>kinase pairs across these three cell types is computed and the significance of the difference in correlations is computed. The results are plotted in a dot plot (figure 6D)


