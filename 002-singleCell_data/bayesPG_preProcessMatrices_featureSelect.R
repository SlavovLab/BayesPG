##### Combined script for working up single cell data for BayesPG
### This involves
## Step 1 = Preprocessing and batch correcting mRNA data ; batch correcting protein data
## Step 2 = Carry out correlation vector analysis to select variable features for cross modality alignments
## Step 3 = Evaluate alignments/compute aligment metrics and make extended data figure 2B
## Step Xtra = since we have mRNA data loaded and ready, also added in code to carry out the ANOVA analysis and make extended data figure 1

### Load libraries
library(Seurat)
library(tidyverse)
library(seqinr)
library(Hmisc)
library(matrixStats)
library(HiClimR)
library(RANN)
library(sva)
library(bluster)
library(data.table)
library(ggpubr)

##### Data preprocessing
#### Reading in Protein and mRNA matrices
### Single cell protein matrices
## Imputed
testesOne <- read.table(file = "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/002-Protein/npop1/npop1_proteinMatrix_Imputed.NoBCByMSRun.txt", sep = "\t", header = TRUE, stringsAsFactors = F)
testesTwo <- read.table(file = "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/002-Protein/npop2/npop2_proteinMatrix_Imputed.NoBC.txt", sep = "\t", header = TRUE, stringsAsFactors = F)
testesThree <- read.table(file = "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/002-Protein/npop3/npop3_proteinMatrix_Imputed.NoBC.txt", sep = "\t", header = TRUE, stringsAsFactors = F)
testesFour <- read.table(file = "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/002-Protein/npop4/npop4_proteinMatrix_Imputed.NoBC.txt", sep = "\t", header = TRUE, stringsAsFactors = F)
testesFive <- read.table(file = "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/002-Protein/npop5/npop5_proteinMatrix_Imputed.BCmTRAQ.txt", sep = "\t", header = TRUE, stringsAsFactors = F)

## Unimputed
testesOne_unimp <- read.table(file = "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/002-Protein/npop1/npop1_proteinMatrix_NoImp.NoBCByMsRun.txt", sep = "\t", header = TRUE, stringsAsFactors = F)
testesTwo_unimp <- read.table(file = "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/002-Protein/npop2/npop2_proteinMatrix_NoImp.NoBC.txt", sep = "\t", header = TRUE, stringsAsFactors = F)
testesThree_unimp <- read.table(file = "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/002-Protein/npop3/npop3_proteinMatrix_NoImp.NoBC.txt", sep = "\t", header = TRUE, stringsAsFactors = F)
testesFour_unimp <- read.table(file = "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/002-Protein/npop4/npop4_proteinMatrix_NoImp.NoBC.txt", sep = "\t", header = TRUE, stringsAsFactors = F)
testesFive_unimp <- read.table(file = "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/002-Protein/npop5/npop5_proteinMatrix_NoImp.BCmTRAQ.txt", sep = "\t", header = TRUE, stringsAsFactors = F)

## paste prep number to get unique cell IDs

colnames(testesOne) <- paste("One", colnames(testesOne), sep = "_")
colnames(testesTwo) <- paste("Two", colnames(testesTwo), sep = "_")
colnames(testesThree) <- paste("Three", colnames(testesThree), sep = "_")
colnames(testesFour) <- paste("Four", colnames(testesFour), sep = "_")
colnames(testesFive) <- paste("Five", colnames(testesFive), sep = "_")

colnames(testesOne_unimp) <- paste("One", colnames(testesOne_unimp), sep = "_")
colnames(testesTwo_unimp) <- paste("Two", colnames(testesTwo_unimp), sep = "_")
colnames(testesThree_unimp) <- paste("Three", colnames(testesThree_unimp), sep = "_")
colnames(testesFour_unimp) <- paste("Four", colnames(testesFour_unimp), sep = "_")
colnames(testesFive_unimp) <- paste("Five", colnames(testesFive_unimp), sep = "_")


### Single cell mRNA matrices
## dataset 1 = 2019 paper, Adults 1 and 2
# Adult 1
testRNA_one <- Read10X(data.dir = "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/001-mRNA/datasetOne_sohni.et.al_2019/adultOne/")
testRNA_one_ready <- CreateSeuratObject(counts = testRNA_one)
testRNA_one_ready[["percent.mt"]] <- PercentageFeatureSet(testRNA_one_ready, pattern = "^MT-")
testRNA_one_ready[["percent.rb"]] <- PercentageFeatureSet(testRNA_one_ready, pattern = "^RP[SL]")
testesOne_RNA <- subset(testRNA_one_ready, subset = nFeature_RNA > 900 & nFeature_RNA < 6000 & percent.mt < 12 & percent.mt > 1)

# Adult 2
testRNA_two <- Read10X(data.dir = "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/001-mRNA/datasetOne_sohni.et.al_2019/adultTwo/")
testRNA_two_ready <- CreateSeuratObject(counts = testRNA_two)
testRNA_two_ready[["percent.mt"]] <- PercentageFeatureSet(testRNA_two_ready, pattern = "^MT-")
testRNA_two_ready[["percent.rb"]] <- PercentageFeatureSet(testRNA_two_ready, pattern = "^RP[SL]")
testesTwo_RNA <- subset(testRNA_two_ready, subset = nFeature_RNA > 900 & nFeature_RNA < 6000 & percent.mt < 12 & percent.mt > 1)

## dataset 2 = 2020 paper
# did not filter as matrices were merged and looked to be pre-filtered (distributions were similar to what we filtered for)
evo <- as.matrix(fread("2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/001-mRNA/datasetTwo_shami.et.al_2020/GSE142585_MergedHumanTestis4_DGE.txt"), rownames = 1)
evo_ready <- CreateSeuratObject(counts = evo)
evo_ready[["percent.mt"]] <- PercentageFeatureSet(evo_ready, pattern = "^MT-")
evo_ready[["percent.rb"]] <- PercentageFeatureSet(evo_ready, pattern = "^RP[SL]")


### Short hand for using a Uniprot ID to Gene Name conversion map/dictionary
# uniprot DB for conversion
uniProtDB <-  read.fasta(file = "2024_Khan.Elcheikhali_testes_rPTR/001-MSData/002-auxiliaryFiles/000-FASTA_etc/uniprot-proteome_UP000005640.fasta")

### Convert uniprot IDs to Entrez Gene Names
## Clean up and get UniprotDB in order
grabHeads <- getAnnot(uniProtDB)
unlistedHeads <- unlist(lapply(grabHeads, `[[`, 1))
splitHeads <- strsplit(unlistedHeads, " +" )
sp_noGN <- splitHeads[-which(grepl(pattern="GN=.*", splitHeads))]
sp_noGN_upIds <- unlist(lapply(strsplit(unlist(lapply(sp_noGN, `[[`, 1)), "\\|" ), `[[`, 2))
splitHeads <- splitHeads[which(grepl(pattern="GN=.*", splitHeads))]

Prot_to_Gene_map <- data.frame(Protein = unlist(lapply(strsplit(unlist(lapply(splitHeads, `[[`, 1)), "\\|" ), `[[`, 2)), Gene = unlist(lapply(splitHeads, grep, pattern="GN=.*", value=TRUE)))
Prot_to_Gene_map$Gene <-gsub("GN=","", Prot_to_Gene_map$Gene)


## Getting different subsets of intersects 
# RNA
intRNA <- Reduce(intersect, list(rownames(testes_one_altWay), rownames(testes_two_altWay), rownames(evo_altWay)))
intProt <- Reduce(intersect, list(rownames(testesOne),rownames(testesTwo),rownames(testesThree), rownames(testesFour), rownames(testesFive)))
intProt_geneName <- unique(Prot_to_Gene_map[which(Prot_to_Gene_map$Protein %in% intProt),]$Gene)
intProtRna <- intersect(intRNA, intProt_geneName)
intProtRNA_uniprot <- unique(Prot_to_Gene_map[which(Prot_to_Gene_map$Gene %in% intProtRna),]$Protein)

### Process proteomic data = batch correction
testesOne_int <- testesOne[which(rownames(testesOne) %in% intProtRNA_uniprot),]
testesTwo_int <- testesTwo[which(rownames(testesTwo) %in% intProtRNA_uniprot),]
testesThree_int <- testesThree[which(rownames(testesThree) %in% intProtRNA_uniprot),]
testesFour_int <- testesFour[which(rownames(testesFour) %in% intProtRNA_uniprot),]
testesFive_int <- testesFive[which(rownames(testesFive) %in% intProtRNA_uniprot),]

mergeProt <- as.matrix(cbind(testesOne_int, testesTwo_int,testesThree_int,testesFour_int, testesFive_int))
rownames(mergeProt) <- Prot_to_Gene_map$Gene[match(rownames(mergeProt), Prot_to_Gene_map$Protein)]

## Batch correct using ComBat
mergeProt_ref <- mergeProt
mergeProt[is.na(mergeProt)] <- 0
mergeProt[mergeProt == -Inf] <- 0

mergeProt_combat <- reshape2::melt(mergeProt, varnames = c("Protein", "id"))
mergeProt_combat <- mergeProt_combat %>% dplyr::mutate("prep" = ifelse(grepl("One", id), "Npop1", ifelse(grepl("Two", id), "Npop2", ifelse(grepl("Three", id), "Npop3", ifelse(grepl("Four", id), "Npop4", "Npop5")))))


batch.covs_mergeProt <- mergeProt_combat$prep[match(colnames(mergeProt), mergeProt_combat$id)]
mergeProt_batch <- ComBat(mergeProt, batch=batch.covs_mergeProt)
mergeProt_batch_ref <- mergeProt_batch

### Ugh, quickly sort out the whole unimp to imp situation
testesOne_unimp_int <- testesOne_unimp[which(rownames(testesOne_unimp) %in% intProtRNA_uniprot),]
testesTwo_unimp_int <- testesTwo_unimp[which(rownames(testesTwo_unimp) %in% intProtRNA_uniprot),]
testesThree_unimp_int <- testesThree_unimp[which(rownames(testesThree_unimp) %in% intProtRNA_uniprot),]
testesFour_unimp_int <- testesFour_unimp[which(rownames(testesFour_unimp) %in% intProtRNA_uniprot),]
testesFive_unimp_int <- testesFive_unimp[which(rownames(testesFive_unimp) %in% intProtRNA_uniprot),]

mergeProt_unimp_int <- as.matrix(cbind(testesOne_unimp_int, testesTwo_unimp_int,testesThree_unimp_int,testesFour_unimp_int, testesFive_unimp_int))
rownames(mergeProt_unimp_int) <- Prot_to_Gene_map$Gene[match(rownames(mergeProt_unimp_int), Prot_to_Gene_map$Protein)]

# This is the batch corrected matrix we have available and read in to our LIGER based dataset alignment
#write.table(mergeProt_batch, "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/002-Protein/Protein_AllFivePreps_batchCorrected.txt", sep = "\t", row.names = TRUE)

mergeProt_batch[is.na(mergeProt_unimp_int)] <- NA

rm(list = c('testesOne_unimp_int','testesTwo_unimp_int', 'testesThree_unimp_int', 'testesFour_unimp_int', 'testesFive_unimp_int', 'testesOne_unimp','testesTwo_unimp', 'testesThree_unimp', 'testesFour_unimp', 'testesFive_unimp'))


### Process mRNA data = normalization and batch correction
## SCTransform (normalization etc)
testes_one_altWay <- SCTransform(testesOne_RNA, vars.to.regress = "percent.mt", variable.features.n = 5000, return.only.var.genes = FALSE)
testes_two_altWay <- SCTransform(testesTwo_RNA, vars.to.regress = "percent.mt", variable.features.n = 5000, return.only.var.genes = FALSE)
evo_altWay <- SCTransform(evo_ready, vars.to.regress = "percent.mt", variable.features.n = 5000, return.only.var.genes = FALSE)


## Merge the objects into one and subset to intersect between mRNA and protein
allRNA_altWay <- merge(evo_altWay, y = c(testes_one_altWay, testes_two_altWay), add.cell.ids = c("evo", "dOne", "dTwo"), project = "Testicles")

allRNA_altWay$human <- unlist(lapply(strsplit(rownames(allRNA_altWay@meta.data),"_"), "[[" , 1))
allRNA_altWay@meta.data$orig.ident[allRNA_altWay@meta.data$human == "dOne"] <- "dOne"
allRNA_altWay@meta.data$orig.ident[allRNA_altWay@meta.data$human == "dTwo"] <- "dTwo"

allRNA_altWay_intProtRNA <- allRNA_altWay
VariableFeatures(allRNA_altWay_intProtRNA) <- intProtRna

## Carry out Anchor based integration
testicles.list_intProtRNA <- SplitObject(allRNA_altWay_intProtRNA, split.by = "human")

# select integration features and prep step
features_intProtRNA <- SelectIntegrationFeatures(testicles.list_intProtRNA)

testicles.list_intProtRNA <- PrepSCTIntegration(
  testicles.list_intProtRNA,
  anchor.features = features_intProtRNA
)

# downstream integration steps
anchors_intProtRNA <- FindIntegrationAnchors(
  testicles.list_intProtRNA,
  normalization.method = "SCT",
  anchor.features = features_intProtRNA
)
testicles.integrated_intProtRNA <- IntegrateData(anchors_intProtRNA, normalization.method = "SCT")




##### Carry out correlation vector analysis to select features to be use for alignment
all_SCTrans_integrate <- as.matrix(testicles.integrated_intProtRNA@assays$integrated@scale.data)

mergeRNACorrs <- Hmisc::rcorr(t(all_SCTrans_integrate_init[order(rownames(all_SCTrans_integrate_init)),]))$r
mergeProt_corrs <- Hmisc::rcorr(t(mergeProt_batch[order(rownames(mergeProt_batch)),]))$r

mrgeRNACorVect <- mergeRNACorrs
mrgeProtCorVect <- mergeProt_corrs

mrgeRNACorVect[mrgeRNACorVect == 1] <- NA
mrgeProtCorVect[mrgeProtCorVect == 1] <- NA

## Hmmm, not convinced by our correlation matrix stuffs, gonna try to actually loop through
RNAProtCor_CorVector <- data.frame(Prot = rownames(mrgeProtCorVect), Cor = NA)

for (i in 1:nrow(RNAProtCor_CorVector)) {
  
  rownames(RNAProtCor_CorVector)[i] <- rownames(mrgeRNACorVect)[i]
  RNAProtCor_CorVector[i,2] <- round(cor(mrgeRNACorVect[i,], mrgeProtCorVect[i,], use = "pairwise.complete.obs"),2)
  
}


RNAProtCor_CorVector_noNA <- RNAProtCor_CorVector[!is.na(RNAProtCor_CorVector$Cor),]
combCorrsAboveThresh <- RNAProtCor_CorVector_noNA[which(RNAProtCor_CorVector_noNA$Cor > quantile(RNAProtCor_CorVector_noNA$Cor, c(0,0.25,0.5, 0.75))[3]),]$Prot

mergeRNACorrs_thresh <- mrgeRNACorVect[rownames(mrgeRNACorVect) %in% combCorrsAboveThresh,colnames(mrgeRNACorVect) %in% combCorrsAboveThresh]
mergeProt_corrs_thresh <- mrgeProtCorVect[rownames(mrgeProtCorVect) %in% combCorrsAboveThresh,colnames(mrgeProtCorVect) %in% combCorrsAboveThresh]


RNAProtCor_CorVectorThresh <- data.frame(Prot = rownames(mergeProt_corrs_thresh), Cor = NA)

for (i in 1:nrow(RNAProtCor_CorVectorThresh)) {
  
  rownames(RNAProtCor_CorVectorThresh)[i] <- rownames(mergeRNACorrs_thresh)[i]
  RNAProtCor_CorVectorThresh[i,2] <- round(cor(mergeRNACorrs_thresh[i,], mergeProt_corrs_thresh[i,], use = "pairwise.complete.obs"),2)
  
}

combCorrsAboveThresh_reThresh <- RNAProtCor_CorVectorThresh[which(RNAProtCor_CorVectorThresh$Cor > 0.29),]$Prot

write.table(combCorrsAboveThresh_reThresh, "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/003-alignmentOutputs/geneProducts_fromCorrelationVectorAnalysis_LigerFeatureSpace.csv", sep = "\t", row.names = FALSE)




#### Evaluating our dataset alignment results
## Read in all required files post alignment
all_RNALabels <- read.table("2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/003-alignmentOutputs/RNA_cellTypeLabels_postAlignment.txt", header = T, stringsAsFactors = F) 
all_ProtLabels <- read.table("2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/003-alignmentOutputs/Protein_cellTypeLabels_postAlignment.txt", header = T, stringsAsFactors = F) 
corrVect_pDIA <- read.table(file = "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/003-alignmentOutputs/geneProducts_fromCorrelationVectorAnalysis_LigerFeatureSpace.csv", sep = "\t", header = TRUE, stringsAsFactors = F)


### Correlate Cluster level means and make heatmap (Ext Fig 2B, Left most panel)
## Pre-process data = join to labels, get cluster level means and Z score
colnames(all_SCTrans_integrate_init) <- gsub("evo_","", colnames(all_SCTrans_integrate_init))
testesTwo_SCTrans_int_melt <-reshape2::melt(all_SCTrans_integrate_init, varnames = c("Protein", "id"), value.name = "Intensity") %>%  
  mutate(id = gsub("-1", "", id)) %>% 
  left_join(all_RNALabels[,c("barcode","OGLabel","TransferLabel")], by = c("id" = "barcode")) %>% 
  filter(TransferLabel != "Macrophage")  %>% 
  rename(cellType = TransferLabel)


RNA_meanLev <- testesTwo_SCTrans_int_melt %>% 
  group_by(cellType, Protein) %>% 
  summarise(meanLev = mean(Intensity, na.rm = T))

mergeProt_batch_melt <- reshape2::melt(mergeProt_batch, varnames = c("Protein", "id"), value.name = "Intensity") %>% 
  left_join(all_ProtLabels[c("id", "cellType")], by = "id")

Prot_meanLev <- mergeProt_batch_melt %>% 
  group_by(cellType, Protein) %>% 
  summarise(meanLev = mean(Intensity, na.rm = T))

## Filter to alignment feature space, make matrices, Z score
RNA_meanLev_filt <- RNA_meanLev %>% filter(Protein %in% corrVect_pDIA$x)
Prot_meanLev_filt <- Prot_meanLev %>% filter(Protein %in% corrVect_pDIA$x)


RNA_meanLev_cast<-reshape2::dcast(RNA_meanLev_filt, Protein ~ cellType, value.var = "meanLev", fill=NA) %>% mutate(Protein = as.character(Protein)) %>% arrange(Protein) %>%  column_to_rownames("Protein") %>% as.matrix 
colnames(RNA_meanLev_cast) <- paste0("RNA_",colnames(RNA_meanLev_cast))
RNA_meanLev_cast <- t(apply(RNA_meanLev_cast, 1, FUN = function(x) {(x - mean(x, na.rm = T))/sd(x, na.rm = T)}))


Prot_meanLev_cast<-reshape2::dcast(Prot_meanLev_filt, Protein ~ cellType, value.var = "meanLev", fill=NA) %>%  mutate(Protein = as.character(Protein)) %>% arrange(Protein) %>% column_to_rownames("Protein") %>% select(-Macrophage) %>%  as.matrix
colnames(Prot_meanLev_cast) <- paste0("Prot_",colnames(Prot_meanLev_cast))
Prot_meanLev_cast <- Prot_meanLev_cast[rownames(Prot_meanLev_cast) %in% rownames(RNA_meanLev_cast),]
Prot_meanLev_cast <- t(apply(Prot_meanLev_cast, 1, FUN = function(x) {(x - mean(x, na.rm = T))/sd(x, na.rm = T)}))

col_order_rna <- c("RNA_EC", "RNA_LC", "RNA_PTM", "RNA_SPG", "RNA_SPC", "RNA_St")
col_order_prot <- c("Prot_EC", "Prot_LC","Prot_PTM", "Prot_SPG", "Prot_SPC","Prot_St")

RNA_meanLev_cast <- RNA_meanLev_cast[, col_order_rna]
Prot_meanLev_cast <- Prot_meanLev_cast[, col_order_prot]


## Compute correlations, reframe and plot
CombMatrices <- cbind(RNA_meanLev_cast,Prot_meanLev_cast)
CombMatrices_cor <- Hmisc::rcorr(CombMatrices)$r 

CombMatrices_cor_melt <- reshape2::melt(CombMatrices_cor, varnames = c("cellRow", "cellCol"), value.name = "Cor") %>% filter(!is.na(Cor)) %>% mutate("prep" = ifelse(grepl("RNA_", cellRow), "RNA","Prot")) 

CombMatrices_cor_reCast <- CombMatrices_cor_melt %>%
  pivot_wider(names_from = cellCol, values_from = Cor)

## Annoying roundaboutway to plot
ProteinOnly <- CombMatrices_cor_reCast %>%
  filter(prep == "Prot") %>%
  select(-prep)

# Rename the rows
ProteinOnly_clean <- ProteinOnly %>% column_to_rownames("cellRow") %>% select(-c(7:12)) %>% as.matrix
CombMatrices_cor_Plot <- reshape2::melt(ProteinOnly_clean, varnames = c("cellRow", "cellCol"), value.name = "Cor") 

# Plot
ggplot(CombMatrices_cor_Plot, aes(cellRow, cellCol, fill = Cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1), name="Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), 
        axis.text.y = element_text(size = 15),
        legend.key.size = unit(2, "lines"),
        legend.title = element_text(size = 15),  
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 20) 
  ) +
  labs(x = NULL, y = NULL, title = "Confusion Matrix of sorts") 



### Compute PCA for Protein single cells, plot and color cells by transferred celltype (Ext Fig 2B, middle panel)
mat.sc.imp <- medInt_cr_norm_log(mergeProt_batch)
X.m <- mat.sc.imp[rownames(mat.sc.imp) %in% corrVect_pDIA$x,]

pca.imp.cor<-Hmisc::rcorr(X.m)$r
pca.imp.cor[is.na(pca.imp.cor)] <- 0

# PCA
sc.pca<-eigen(pca.imp.cor)
scx<-as.data.frame(sc.pca$vectors)
colnames(scx)<-paste0("PC",1:ncol(scx))
scx$id<-colnames(pca.imp.cor)

# Percent of variance explained by each principle component
pca_var <- sc.pca$values
percent_var<- pca_var/sum(pca_var)*100
plot(1:length(percent_var), percent_var, xlab="PC", ylab="% of variance explained")

scx$prep <- mergeProt_combat$prep[match(scx$id, mergeProt_combat$id)]
PCA_final <- scx %>% left_join(all_ProtLabels, by = "id")

ggscatter(PCA_final, x = 'PC1', y = 'PC2', color = 'cellType', size = 1, alpha=0.5) + 
  xlab(paste0("PC1 (",round(percent_var[1],2), "%)")) +
  ylab(paste0("PC2 (",round(percent_var[2],2), "%)")) + 
  font("ylab",size=30) +
  font("xlab",size=30) +
  font("xy.text", size=20) +
  annotate("text", x=-0.025, y=-0.045, label=paste0(dim(X.m)[1], " proteins"), size=8) +
  annotate("text", x=-0.025, y=-0.035, label=paste0(dim(X.m)[2], " cells"), size=8)



#### Compute Distance Ratios and plot old vs new alignments (Ext Fig 2B, Right most panel))
# Decide on threshold and compute Mahalonobis of sorts
perVarthresh <- length(percent_var[cumsum(percent_var) < 25])
PCA_ForDists <- scx[, c(1:perVarthresh, which(colnames(scx) == "id"))]
rownames(PCA_ForDists) <- PCA_ForDists$id
#PCA_ForDists <- PCA_ForDists[,c(1:perVarthresh)]
PCA_ForDists <- PCA_ForDists[,c(1:2)]

dist_matrix <- as.matrix(dist(PCA_ForDists)) 

# Some sanity chex and things
nPoP_cellIds_new <- all_ProtLabels[match(colnames(dist_matrix), all_ProtLabels$id),]
nPoP_cellIds_old <- nPoP_cellIds2[match(colnames(dist_matrix), nPoP_cellIds2$id),]

## Define functions to compute distance ratios
# Cleanly collect median distance within each cluster
median_within_cluster_distance <- function(dist_matrix, clusters) {
  median_distances <- sapply(unique(clusters), function(cluster) {
    cluster_indices <- which(clusters == cluster)
    if (length(cluster_indices) > 1) {
      cluster_distances <- dist_matrix[cluster_indices, cluster_indices]
      return(median(as.vector(cluster_distances)))
    } else {
      return(0)  
    }
  })
  names(median_distances) <- unique(clusters)
  return(median_distances)
}

# Cleanly compute median distances between clusters
median_between_cluster_distance <- function(dist_matrix, clusters) {
  cluster_pairs <- combn(unique(clusters), 2, simplify = FALSE)
  median_distances <- sapply(cluster_pairs, function(pair) {
    cluster1_indices <- which(clusters == pair[1])
    cluster2_indices <- which(clusters == pair[2])
    between_distances <- dist_matrix[cluster1_indices, cluster2_indices]
    return(median(as.vector(between_distances)))
  })
  result <- data.frame(cluster1 = sapply(cluster_pairs, function(x) x[1]),
                       cluster2 = sapply(cluster_pairs, function(x) x[2]),
                       median_distance = median_distances)
  return(result)
}

# Function to combine the two and compute actual ratios for one dataset
compute_distance_ratios <- function(dist_matrix, clusters) {
  
  median_within <- median_within_cluster_distance(dist_matrix, clusters)
  median_between <- median_between_cluster_distance(dist_matrix, clusters)
  overall_median_within <- median(median_within)
  overall_median_between <- median(median_between$median_distance)
  
  ratio <- overall_median_within / overall_median_between
  return(list(ratio = ratio,
              median_within = overall_median_within,
              median_between = overall_median_between,
              detailed_within = median_within,
              detailed_between = median_between))

  }

# sort of annoying, but, final function to process through distance ratio results into plot format
compute_ratio_result <- function(dist_matrix, cluster_ids, clustering_label) {
  
  ratio_result <- compute_distance_ratios(dist_matrix, cluster_ids)
  
  ratio_result_within <- ratio_result$detailed_within
  ratio_result_between <- ratio_result$detailed_between
  ratio_result_between$distRatio <- NA
  
  for (i in seq_along(names(ratio_result_within))) {
    cluster_name <- names(ratio_result_within)[i]
    within_ratio <- ratio_result_within[[i]]
    
    ratio_result_between[ratio_result_between$cluster1 == cluster_name, ]$distRatio <- 
      within_ratio / ratio_result_between[ratio_result_between$cluster1 == cluster_name, ]$median_distance
  }
  
  ratio_result_between %>%
    dplyr::mutate(Clustering = clustering_label)
}

ratio_result_between_LIGER <- compute_ratio_result(dist_matrix, nPoP_cellIds_new$Cluster, "LIGER_New")
ratio_result_between_LIGEROld <- compute_ratio_result(dist_matrix, nPoP_cellIds_old$Cluster, "LIGER_Old")

ratioResults_both <- dplyr::bind_rows(ratio_result_between_LIGER, ratio_result_between_LIGEROld)

## So this is the final plot
# Need to revisit this...think something was differetn here
ggplot(ratioResults_both,aes(x=distRatio, fill=Clustering)) + geom_density(alpha=0.25) + geom_vline(xintercept = 1, color = "red", linetype = "solid")





##### ANOVA across mRNA datasets with dataset, human and celltypes as covariates 
## Read in cell type labels = here, I've just loaded the combined table of labels after all the integration as it is simpler to work with...there is enough tedium in here as is
# We used the original/OG labels for this analysis
all_RNALabels[all_RNALabels$Dataset == "neoAdult",]$barcode <- paste0(all_RNALabels[all_RNALabels$Dataset == "neoAdult",]$barcode,"-1")

## Some preprocessing and labelling etc
testicles.integrated_intProtRNA_matrix <- as.matrix(GetAssayData(testicles.integrated_intProtRNA[["integrated"]], slot = "data"))
testicles.integrated_intProtRNA_matrix_startProcess<- testicles.integrated_intProtRNA_matrix
testicles.integrated_intProtRNA_df <- as.data.frame(t(testicles.integrated_intProtRNA_matrix_startProcess)) %>% 
  rownames_to_column("id") %>% 
  mutate("prep" = ifelse(grepl("evo_", id), "new", "old"), "human" = ifelse(grepl("_Human1.", id), "HumanOne", ifelse(grepl("_Human2.", id), "HumanTwo", ifelse(grepl("_Human3.", id), "HumanThree", ifelse(grepl("_Human4.", id), "HumanFour", ifelse(grepl("dOne", id), "H1", ifelse(grepl("dTwo", id), "H2", "")))))))

# we removed human 4 from the second dataset (2020) from our analysis as it's mt reads, nFeatures etc looked significantly worse than other humans and replicates
testicles.integrated_intProtRNA_df_across <- testicles.integrated_intProtRNA_df %>% 
  filter(human != "HumanFour") %>% 
  reshape2::melt()

testicles.integrated_intProtRNA_df_across[testicles.integrated_intProtRNA_df_across$prep == "new",]$id <- substr(testicles.integrated_intProtRNA_df_across[testicles.integrated_intProtRNA_df_across$prep == "new",]$id, 5, nchar(testicles.integrated_intProtRNA_df_across[testicles.integrated_intProtRNA_df_across$prep == "new",]$id))

## join data matrix to cell type labels
testicles.integrated_intProtRNA_df_across <- testicles.integrated_intProtRNA_df_across %>% 
  left_join(all_RNALabels, by = c("id" = "barcode")) 

testicles.integrated_intProtRNA_df_across <- testicles.integrated_intProtRNA_df_across %>% 
  filter(OGLabel %in% c("EC", "LC", "PTM","SPG","SPC","St")) 

### Carry out ANOVA
hereWeJO_integrated_again <- aov(value ~ prep + human + OGLabel, data = testicles.integrated_intProtRNA_df_across)
hereWeJO_interated_again_summary <- summary(hereWeJO_integrated_again)

## Plot results
sum_sq <- hereWeJO_interated_again_summary[[1]]$`Sum Sq`

sum_of_squares <- c(
  prep = sum_sq[1],
  human = sum_sq[2],
  cellType = sum_sq[3]
)

barplot(
  sum_of_squares,
  main = "Evaluating Impact of Sources of Variance",
  ylab = "Sum of Squares Between Groups",
  names.arg = c("Prep", "Human", "Cell Type"),
  col = c("skyblue", "salmon", "lightgreen"), 
  cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.names = 1.5
)




