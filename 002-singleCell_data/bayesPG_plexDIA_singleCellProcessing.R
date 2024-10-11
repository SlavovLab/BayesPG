##### Script to process through plexDIA data/npop5 acquired on the TIMS TOF SCP
### Mapping samples prepared via nPoP and processing is done through a hybrid of the QuantQC package and some specific code


## Load libraries
library(QuantQC)
source("functions_parameters_speedUpdate.R")

## locations for raw files
all_cells <- "2024_Khan.Elcheikhali_testes_rPTR/001-MSData/002-auxiliaryFiles/001-datasetMapping/cellenOneAndMSRunMappingMetadata_for.npop5Only/npop5_cellenOneIsolationFile_cellSizesEtc.xls"
data_path <- "2024_Khan.Elcheikhali_testes_rPTR/001-MSData/001-searchedFiles/npop5_DIANN_report.tsv"
link_path <- "2024_Khan.Elcheikhali_testes_rPTR/001-MSData/002-auxiliaryFiles/001-datasetMapping/cellenOneAndMSRunMappingMetadata_for.npop5Only/npop5_annotation.csv"

## assign some variables in advance
all_cells <- list(all = all_cells)
plex_used <- c(0, 4, 8)


### part of Diann to QQC
# essentially, read in tables and set up/munge mostly 
linker <- read.csv(link_path)
linker$Order <- 1:nrow(linker)
columns_to_read <- c("Genes", "Run", "Lib.PG.Q.Value", "RT", 
                     "Precursor.Id", "Stripped.Sequence", "Precursor.Charge", 
                     "Precursor.Quantity", "Ms1.Area", "Protein.Group", "Translated.Q.Value", 
                     "Channel.Q.Value", "Proteotypic")

Raw_data <- data.table::fread(data_path, select = columns_to_read)

Raw_data <- as.data.frame(Raw_data)
Raw_data <- Raw_data %>% filter(Run %in% linker$Run)
Raw_data <- Raw_data %>% left_join(linker, by = c("Run"))
Raw_data$seqcharge <- paste0(Raw_data$Stripped.Sequence, 
                             Raw_data$Precursor.Charge)
Raw_data <- Raw_data %>% filter(Protein.Group != "")
Raw_data$plex <- substr(Raw_data$Precursor.Id[1:nrow(Raw_data)], 
                        10, 10)
Raw_data$ID <- paste0(Raw_data$Well, Raw_data$plate, ".", 
                      Raw_data$plex)
Raw_data$File.Name <- Raw_data$ID
Raw_data$uq <- paste0(Raw_data$File.Name, Raw_data$Protein.Group, 
                      Raw_data$seqcharge)
Raw_data <- Raw_data %>% distinct(uq, .keep_all = T)
Raw_data$uq <- NULL


### Part of cell x peptide
# Again, move processing forward, filter if you'd like, decide which column you want to use
Raw_data_filt <- Raw_data %>% filter(plex %in% plex_used, Channel.Q.Value < 1, Translated.Q.Value < 1)

Raw_data_lim_filt <- Raw_data_filt %>% dplyr::select(Protein.Group, 
                                                     seqcharge, Ms1.Area, File.Name)

Raw_data_lim.d_filt <- reshape2::dcast(Raw_data_lim_filt, Protein.Group + seqcharge ~ File.Name, value.var = "Ms1.Area")
Raw_data_lim.d_filt[Raw_data_lim.d_filt == 0] <- NA

### Do CellenOne to QuantQC object mapping
cellenOne_data <- analyzeCellenONE_mTRAQ(all_cells, plex = 3)

peptide_data <- Raw_data_lim.d_filt[,-c(1:2)]
cellID <- colnames(peptide_data)
cellID <- as.data.frame(cellID)
colnames(cellID) <- "ID"
cellenOne_data_small <- cellenOne_data %>% dplyr::select(any_of(c("ID", 
                                                                  "diameter", "sample", "label", "injectWell", "plate")))
cellenOne_data_small <- as.data.frame(cellenOne_data_small)
cellID <- cellID %>% left_join(cellenOne_data_small, by = c("ID"))
cellID$sample[is.na(cellID$sample) == T] <- "neg"
cellID$prot_total <- log2(colSums(peptide_data[, 1:ncol(peptide_data)], 
                                  na.rm = T))

cellID$WP <- paste0(cellID$plate, cellID$injectWell)
cellID$plate <- NULL
linker$WP <- paste0(linker$plate, linker$Well)
cellID <- cellID %>% left_join(linker, by = c("WP"))
cellID$WP <- NULL


##  Evaluate negative controls
negCells <- cellID[which(cellID$sample == "neg"),]$ID

reportFiltJoined_casted <- Raw_data_lim.d_filt 
rownames(reportFiltJoined_casted) <- Raw_data_lim.d_filt$seqcharge 
oneRing_matrix <- as.matrix(reportFiltJoined_casted %>% select(-c(seqcharge, Protein.Group)))

## Idiotic munging just to use matrix normalisation function
oneRing_matrix_normd <- medInt_cr_norm(oneRing_matrix)

oneRing_remelt <- reshape2::melt(oneRing_matrix_normd, varnames = c("seqcharge", "id"))
oneRing_CVthings <- oneRing_remelt %>% left_join(Raw_data_lim_filt[,c("Protein.Group", "seqcharge","File.Name")], by = c("seqcharge", c("id" = "File.Name")))
oneRing_CVthings <- oneRing_CVthings[!is.na(oneRing_CVthings$value),]

xd4<- oneRing_CVthings %>%
  group_by(id, Protein.Group,) %>%
  mutate(cvq = cv(value)) %>% ungroup()

xd5<- xd4 %>%
  group_by(Protein.Group, id) %>%
  mutate(cvn = cvna(value))  %>% ungroup()

xd6<-xd5 %>% group_by(id) %>% mutate(cvm=median(cvq, na.rm=T)) %>% ungroup()

xd7<-xd6 %>% group_by(id,Protein.Group) %>% mutate(nPep= length(unique(seqcharge))) %>% ungroup()

xd8<-xd7 %>% group_by(id) %>% summarise(nProts = n_distinct(Protein.Group[nPep > 3 & !is.na(value)]), cvm = cvm,signal = sum(value, na.rm = T), totProt = n_distinct(Protein.Group)) %>% distinct()

xd9 <- xd8 %>% left_join(cellID[,c("ID","diameter","sample")], by = c("id" = "ID")) %>% distinct

## This is the final visualization that we use to decide on thresholds
ggscatter(data = xd9, x = 'cvm', y = 'nProts', color = 'sample') +
  geom_vline(xintercept = 0.65) +
  geom_hline(yintercept = 60) + 
  ylab("nProts [>3 peps]") + 
  xlab ("Median CV/Cell") + theme(text = element_text(size = 18)) 

cleanCells <- xd9[which(xd9$nProts > 58 & xd9$cvm < 0.63),]$id
#cleanestCells <- xd9[which(xd9$nProts > 65 & xd9$cvm < 0.57),]$id


### filter for cells per CV and number of proteins with > 3 peptides
# We then normalize the peptide level matrix and also have to batch correct for mTRAQ label
Raw_data_lim_filt_noNeg <- Raw_data_lim.d_filt[,c(1:2,which(colnames(Raw_data_lim.d_filt) %in% cleanCells))]
Raw_data_lim_filt_noNeg <- Raw_data_lim_filt_noNeg[!grepl(";", Raw_data_lim.d_filt$Protein.Group),]

prot_pep <- Raw_data_lim_filt_noNeg[,c(1:2)]
Raw_data_lim_filt_noNeg_t1 <- medInt_cr_norm_log(log2(as.matrix(Raw_data_lim_filt_noNeg[,-c(1:2)])))

trem_premelt_equivalent <- cbind(prot_pep,Raw_data_lim_filt_noNeg_t1)
trem_postMelt_equivalent <- setDT(reshape2::melt(trem_premelt_equivalent))
t4 <- trem_postMelt_equivalent[,list(qp=median(value, na.rm=T)),by = c("Protein.Group","variable")]
t4m<-reshape2::dcast(t4, Protein.Group ~ variable, value.var = "qp", fill=NA)
t4m <- column_to_rownames(t4m, "Protein.Group")

sc.imp_peptides <- hknn(as.matrix(Raw_data_lim_filt_noNeg_t1), k.t)

sc.imp_peptides[is.nan(sc.imp_peptides) == T] <- 0
sc.imp_peptides[is.infinite(sc.imp_peptides) == T] <- 0

batch.covs_peptides <- substr(colnames(sc.imp_peptides), nchar(colnames(sc.imp_peptides)), nchar(colnames(sc.imp_peptides)))

matrix.sc.batch_peptides <- ComBat(sc.imp_peptides, batch=batch.covs_peptides)
matrix.sc.batch_imp_peptides <- matrix.sc.batch_peptides
matrix.sc.batch_peptides[is.na(Raw_data_lim_filt_noNeg_t1)==T] <- NA

peptideLevel_BatchCorrected <- cbind(prot_pep, matrix.sc.batch_peptides)

write.table(peptideLevel_BatchCorrected, "2024_Khan.Elcheikhali_testes_rPTR/002-SingleCellMatrices/002-Protein/npop5/npop5_peptideMatrix_NoImp.BCmTRAQ.txt", sep = "\t", row.names = TRUE)


## Impute protein level data, again for batch correction for mTRAQ labels
imp.input<-as.matrix(t4m)
sc.imp <- hknn(imp.input, k.t)

t5<-sc.imp
sum(is.na(sc.imp))
dim(sc.imp)

sc.imp[(is.na(sc.imp))]<-0

batch.covs <- substr(colnames(sc.imp), nchar(colnames(sc.imp)), nchar(colnames(sc.imp)))

matrix.sc.batch <- ComBat(sc.imp, batch=batch.covs)
matrix.sc.batch_imp <- matrix.sc.batch
matrix.sc.batch[is.na(imp.input)==T] <- NA

## Write out tables
write.table(matrix.sc.batch_imp, "2024_Khan.Elcheikhali_testes_rPTR/002-SingleCellMatrices/002-Protein/npop5/npop5_proteinMatrix_Imputed.BCmTRAQ.txt", sep = "\t", row.names = TRUE)

write.table(matrix.sc.batch, "2024_Khan.Elcheikhali_testes_rPTR/002-SingleCellMatrices/002-Protein/npop5/npop5_proteinMatrix_NoImp.BCmTRAQ.txt", sep = "\t", row.names = TRUE)





