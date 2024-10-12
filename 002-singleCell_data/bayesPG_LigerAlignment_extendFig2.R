##### Script for LIGER alignment
### Includes: 
## 1) Basic dataset preprocessing and eventual LIGER alignment
## 2) Collecting and transferring celltype labels post alignments
## 3) Figures used in Extended Data Figure 2

### load in all required packages
library(tidyverse)  
library(rliger)
library(seqinr)
library(reshape2)
library(Seurat)
library(ggpubr)
library(data.table)

#### Read in data and preprocess etc
### Read in mRNA data
# the data fron Sohni et al contained the 10x directory and so could be directly read in
dOne <- read10X(sample.dirs =list("2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/001-mRNA/datasetOne_sohni.et.al_2019/adultOne"), sample.names = c("dOne"))
dTwo <-  read10X(sample.dirs =list("2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/001-mRNA/datasetOne_sohni.et.al_2019/adultTwo"), sample.names = c("dTwo"))

# the data from Shami et al consisted of a processed matrix and so required a little processing
evo <- as.matrix(fread("2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/001-mRNA/datasetTwo_shami.et.al_2020/GSE142585_MergedHumanTestis4_DGE.txt"), rownames = 1)
evo_ready <- CreateSeuratObject(counts = evo)
evo_ready[["percent.mt"]] <- PercentageFeatureSet(evo_ready, pattern = "^MT-")
evo_ready[["percent.rb"]] <- PercentageFeatureSet(evo_ready, pattern = "^RP[SL]")
evo_ready <- RenameCells(evo_ready[["RNA"]],new.names = paste0("evo_", colnames(x = evo_ready[["RNA"]])))

## Read in Protein data = Single batch corrected matrix across the five preps
testesAllFiveComb <- read.table(file = "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/002-Protein/Protein_AllFivePreps_batchCorrected.txt", sep = "\t", header = TRUE, stringsAsFactors = F)

# Our protein matrix is log transformed and so, we exponentiate to satisfy non negativity for LIGER
testesAllFiveComb <- 2^testesAllFiveComb

### Read in variable features selected via correlation vector analysis
corrVect_pDIA <- read.table(file = "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/003-alignmentOutputs/geneProducts_fromCorrelationVectorAnalysis_LigerFeatureSpace.csv", sep = "\t", header = TRUE, stringsAsFactors = F)


### Teensy bit of preprocess and just making sure to have failsafes in case of random exports
# important to note here, that we have converted the batch corrected Protein IDs to Gene Names before export
testesOne_Gene <- data.matrix(testesAllFiveComb)
testesOne_Gene <- testesOne_Gene[!rowSums(is.na(testesOne_Gene)),]

featureSpace_inProteindata <- rownames(testesOne_Gene[which(rownames(testesOne_Gene) %in% corrVect_pDIA$x),])
features_rna_prot = intersect(rownames(dTwo), featureSpace_inProteindata) 


#### Dataset alignment via LIGER
# In order to push our data in to LIGER we had to hack it in (so to speak, big thanks to Chao Gao from Welch lab for help with this!)
testes_all <- createLiger(list(dTwo = dTwo[features_rna_prot,], dOne = dOne[features_rna_prot,],evo = evo_ready[features_rna_prot,], testes1 = testesOne_Gene))

testes_all <- normalize(testes_all) 
testes_all@var.genes <- rownames(testes_all@norm.data[[1]])
testes_all <- scaleNotCenter(testes_all) 
testes_all@scale.data$testes1 <- t(testesOne_Gene[testes_all@var.genes,]) 


### now we can carry out our factorization etc
# here, we did play around a lot with the settings not only for number of factors and regularization parameter but also for the louvain clustering, but, eventially settled on what we have here, hehe.
testes_all <- optimizeALS(testes_all, k = 10) 
testes_all <- quantile_norm(testes_all)
testes_all <- louvainCluster(testes_all)

# relatively useful when initially iterating through alignment parameters
calcAlignment(testes_all)
calcAlignment(testes_all, by.dataset = T)
calcAgreement(testes_all)
calcAgreement(testes_all, by.dataset = T)


## we just get this ready for now as we use the output to make our figures (amongst other things) later
testes_all <- runUMAP(testes_all)
plots <- plotByDatasetAndCluster(testes_all, return.plots = T, axis.labels = c("UMAP1", "UMAP2"), legend.fonts.size = 16, pt.size = 0.5)




#### A chunky messy bit to map/munge cell annotations
rnaOne_cellAnnot <- read.csv("2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/001-mRNA/GSE124263_Barcode_Cell_Identity.csv", header=TRUE, stringsAsFactors = F)
rnaTwo_cellAnnot <- fread("2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/001-mRNA/datasetTwo_shami.et.al_2020/GSE142585_MergedHumanTestis4_PerCellAttributes.txt")

initialAnnotations <- c("(SC)","(LC)","(PTM)","(GC)","EC","(PGCL)","(SPG)","(St)","(SPC)","(SSC)","(TC)","Diff SPG","Early SPG","Macrophage")
df <- as.data.frame(rnaOne_cellAnnot)
cleanerAnnotations <- str_extract_all(df$AdultCellAnnotation, str_c(initialAnnotations, collapse="|"), simplify = TRUE)[,1]
df$type <- cleanerAnnotations

adultOne <- df[df$orig.ident=="A1_T",]
adultTwo <- df[df$orig.ident=="A2_T",]

## Subset UMAP data
umap_ByCluster <- plots[[1]]$data 
umap_ByCluster$id <- rownames(umap_ByCluster)
rownames(umap_ByCluster) <- NULL
umap_ByCluster_prot <- umap_ByCluster[which(umap_ByCluster$Dataset == 'testes1'),]
umap_ByCluster_rna <- umap_ByCluster[-which(umap_ByCluster$Dataset == 'testes1'),]

## more munge
rnaAnot_work <- umap_ByCluster_rna
rnaAnot_work[rnaAnot_work$Dataset == "evo",]$id <- substr(rnaAnot_work[rnaAnot_work$Dataset == "evo",]$id, 5, nchar(rnaAnot_work[rnaAnot_work$Dataset == "evo",]$id))
rnaAnot_work[!rnaAnot_work$Dataset == "evo",]$id <- substr(rnaAnot_work[!rnaAnot_work$Dataset == "evo",]$id, 6, nchar(rnaAnot_work[!rnaAnot_work$Dataset == "evo",]$id))
rnaAnot_work <- rename(rnaAnot_work, "barcode" = "id")
rnaAnot_work$barcode<-gsub(str_c(c("^dOne_","^dTwo_"),collapse="|"),"",rnaAnot_work$barcode)
adultOne$Cell.Barcode<-gsub("A1T_","",adultOne$Cell.Barcode)
adultTwo$Cell.Barcode<-gsub("A2T_","",adultTwo$Cell.Barcode)
rnaAnot_work$type <- NA
forType <- data.frame(c(adultOne$type,adultTwo$type),c(adultOne$Cell.Barcode,adultTwo$Cell.Barcode))
colnames(forType) <- c("type","barcode")
rnaAnot_work$type <- forType$type[match(rnaAnot_work$barcode,forType$barcode,incomparables=TRUE)]
rnaAnot_work[rnaAnot_work$Dataset == "evo",]$type <- rnaTwo_cellAnnot$CellType[match(rnaAnot_work[rnaAnot_work$Dataset == "evo",]$barcode,rnaTwo_cellAnnot$V1,incomparables=TRUE)]

## Merge cell type annotations: currently at different levels of differentiation, were collapsed to major celltypes
rnaAnot_work$type <- gsub("Spermatogonia", "SPG", rnaAnot_work$type)
rnaAnot_work$type <- gsub("RoundSpermatid", "St", rnaAnot_work$type)
rnaAnot_work$type <- gsub("Elongating", "St", rnaAnot_work$type)
rnaAnot_work$type <- gsub("Spermatocyte", "SPC", rnaAnot_work$type)
rnaAnot_work$type <- gsub("Myoid", "PTM", rnaAnot_work$type)
rnaAnot_work$type <- gsub("Endothelial", "EC", rnaAnot_work$type)
rnaAnot_work$type <- gsub("ImmLeydig", "LC", rnaAnot_work$type)
rnaAnot_work$type <- gsub("Tcell", NA, rnaAnot_work$type)
rnaAnot_work$type <- gsub("f-Pericyte", NA, rnaAnot_work$type)
rnaAnot_work$type <- gsub("m-Pericyte", NA, rnaAnot_work$type)
rnaAnot_work$type <- gsub("Macrophage", NA, rnaAnot_work$type)


## For each celltype, get fractions per louvain cluster etc
finalCellTypes <- unique(rnaAnot_work$type)
finalCellTypes <- finalCellTypes[!is.na(finalCellTypes)]
cellTypeFrac_DF <- data.frame(matrix(ncol=length(finalCellTypes)+1,nrow=length(unique(rnaAnot_work$Cluster))))
colnames(cellTypeFrac_DF)[1] <- "Cluster"
cellTypeFrac_DF$Cluster <- unique(rnaAnot_work$Cluster)
cellTypeFrac_DF$Cluster <- as.numeric(as.character(cellTypeFrac_DF$Cluster))
cellTypeFrac_DF <- arrange(cellTypeFrac_DF,Cluster)
for(X in 1:length(finalCellTypes)){
  colnames(cellTypeFrac_DF)[X+1] <- as.character(finalCellTypes[X])
}

for(X in cellTypeFrac_DF$Cluster){
  print(X)
  total <- length(which(rnaAnot_work$Cluster==X))
  for(type in 1:length(finalCellTypes)){
    cellTypeFrac_DF[cellTypeFrac_DF$Cluster==X,type+1] <- round((length(which(rnaAnot_work$type[rnaAnot_work$Cluster==X]==finalCellTypes[type]))/total)*100,3)
  }
}

## Get primary/majority cell type per cluster and sort final DF
cellTypeFrac_DF_melt <- melt(cellTypeFrac_DF,id='Cluster', variable.name = "cellType", value.name = "percOfClust")
cellTypeFrac_DF_melt$Cluster <- as.character(cellTypeFrac_DF_melt$Cluster)
colnames(cellTypeFrac_DF_melt)[3] <- "Percent"

clusterToMainCell <- cellTypeFrac_DF_melt %>% group_by(Cluster) %>% slice(which.max(Percent))
clusterToMainCell <- umap_ByCluster_prot %>% left_join(clusterToMainCell, by = 'Cluster') 
clusterToMainCell_forTransferChex <-clusterToMainCell # for use later

### Final bits of processing to get formatted outputs
## mRNA
df_werx <- df
df_werx <- df %>% rename("cellType" = "type")
df_werx$Cell.Barcode<-gsub("A1T_","dOne_",df_werx$Cell.Barcode)
df_werx$Cell.Barcode<-gsub("A2T_","dTwo_",df_werx$Cell.Barcode)

rnaAnot_transfer <- rnaAnot_work
rnaAnot_transfer[rnaAnot_transfer$Dataset == "dOne",]$barcode <- paste0("dOne_",rnaAnot_transfer[rnaAnot_transfer$Dataset == "dOne",]$barcode)
rnaAnot_transfer[rnaAnot_transfer$Dataset == "dTwo",]$barcode <- paste0("dTwo_",rnaAnot_transfer[rnaAnot_transfer$Dataset == "dTwo",]$barcode)

transferred_RNA <- rnaAnot_transfer %>%  
  left_join(rnaTwo_cellAnnot[,c("V1","CellType")], by = c("barcode" = "V1")) %>% 
  mutate(barcode = gsub("-1", "", barcode)) %>% 
  left_join(df_werx[c("Cell.Barcode", "cellType")], by = c("barcode" = "Cell.Barcode")) %>% 
  unite(col = cells, c("cellType","CellType"), na.rm = T) %>% 
  filter(!cells %in% c("Tcell","f-Pericyte","m-Pericyte")) 

transferred_RNA$cells <- gsub("Spermatogonia", "SPG", transferred_RNA$cells)
transferred_RNA$cells <- gsub("RoundSpermatid", "St", transferred_RNA$cells)
transferred_RNA$cells <- gsub("Elongating", "St", transferred_RNA$cells)
transferred_RNA$cells <- gsub("Spermatocyte", "SPC", transferred_RNA$cells)
transferred_RNA$cells <- gsub("Myoid", "PTM", transferred_RNA$cells)
transferred_RNA$cells <- gsub("Endothelial", "EC", transferred_RNA$cells)
transferred_RNA$cells <- gsub("ImmLeydig", "LC", transferred_RNA$cells)

rnaAnot_work[rnaAnot_work$Dataset == "dOne",]$barcode <- paste0("dOne_",rnaAnot_work[rnaAnot_work$Dataset == "dOne",]$barcode)
rnaAnot_work[rnaAnot_work$Dataset == "dTwo",]$barcode <- paste0("dTwo_",rnaAnot_work[rnaAnot_work$Dataset == "dTwo",]$barcode)

rnaAnot_work[rnaAnot_work$Dataset %in% c("dOne","dTwo"),]$Dataset <- "neoAdult"

rnaAnot_work <- rnaAnot_work %>% 
  select(Dataset, Cluster, barcode, type) %>% 
  mutate(Cluster = as.character(Cluster))

RNALabels_clusterLabels <- rnaAnot_work %>% left_join(clusterToMainCell_forTransferChex[,c("Cluster","cellType")], by = "Cluster") %>% rename(OGLabel = type, TransferLabel = cellType)

write.table(rnaAnot_work, "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/003-alignmentOutputs/RNA_cellTypeLabels_postAlignment.txt", sep = "\t")

## Protein

## select and cell up cell id per cluster
clusterToMainCell <- clusterToMainCell %>% 
  dplyr::mutate("Dataset" = ifelse(grepl("dOne_", id), "dOne",ifelse(grepl("dTwo_", id), "dTwo",ifelse(grepl("One_", id), "testes1",ifelse(grepl("Two_", id), "testes2",ifelse(grepl("Three_", id), "testes3", ifelse(grepl("Four_", id), "testes4", "testes5")))))))

# So, there is a difference here, I have not included the workflow that we used for the NN celltypes, as we haven't fully benchmarked and worked that through and so best to leave it for now, plz.ignore
clusterToMainCell_protOut <- clusterToMainCell %>% dplyr::select(Dataset, Cluster,id,cellType)

write.table(clusterToMainCell_protOut, "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/003-alignmentOutputs/Protein_cellTypeLabels_postAlignment.txt", sep = "\t")




#### Make plots for extended data figure 2
rnaAnot_forPlotting <- rnaAnot_work

rnaAnot_forPlotting$type <- gsub("SPG", "S'Gonia", rnaAnot_forPlotting$type)
rnaAnot_forPlotting$type <- gsub("SPC", "S'Cyte", rnaAnot_forPlotting$type)
rnaAnot_forPlotting$type <- gsub("St", "S'tid", rnaAnot_forPlotting$type)
rnaAnot_forPlotting$type <- gsub("PTM", "Myoid", rnaAnot_forPlotting$type)
rnaAnot_forPlotting$type <- gsub("LC", "Endo", rnaAnot_forPlotting$type)
rnaAnot_forPlotting$type <- gsub("EC", "Leydig", rnaAnot_forPlotting$type)


## Figure with cells by data type (mRNA/Prot)
byPrep <- umap_ByCluster %>% mutate(plotByModality = NA)
byPrep[which(byPrep$Dataset == "testes1"),]$plotByModality <- "scProt"
byPrep[-which(byPrep$Dataset == "testes1"),]$plotByModality <- "scRNA"

byPrep$plotByPrep <- NA
byPrep[which(byPrep$Dataset == "testes1"),]$plotByPrep <- "scProt"
byPrep[which(byPrep$Dataset == "dOne"),]$plotByPrep <- "scRNA-1"
byPrep[which(byPrep$Dataset == "dTwo"),]$plotByPrep <- "scRNA-1"
byPrep[which(byPrep$Dataset == "evo"),]$plotByPrep <- "scRNA-2"

## lol, paper colormap
ggplot(byPrep, aes(x=tsne1, y=tsne2, color=plotByModality)) +
  geom_point(size=0.5) +  # Add points
  labs(y = "UMAP2", x = "UMAP1", color = "Dataset") +  # Set axis labels
  theme(axis.title = element_text(size = 18),  # Adjust theme elements
        axis.text = element_text(size = 16),
        legend.key.size = unit(1.2, 'cm'),
        legend.title = element_text(size=17),
        legend.text = element_text(size=16)) +
  # Define color scheme
  scale_color_manual(values = c("#03b1fc", "#edade8")) +  # Pastel blue, purple
  guides(color = guide_legend(override.aes = list(size=2)))  # Adjust legend size


## Figure by transferred celltypes
ggplot(rnaAnot_forPlotting, aes(x=tsne1, y=tsne2, color=type)) +
  geom_point(size=1) + 
  labs(y = "UMAP2", x = "UMAP1", color = "CellTypes") + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.key.size = unit(1.2, 'cm'),
        legend.title = element_text(size=17),
        legend.text = element_text(size=16)) +
  guides(color = guide_legend(override.aes = list(size=2)))