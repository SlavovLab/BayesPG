##### Script to carry out PTM/Phospho analysis from Figure 6 in the Manuscript
### Broadly, each section of the script has to do with making panels b,c,d
library(tidyverse)
library(reshape2)
library(seqinr)
library(ggtext)
library(Hmisc)
library(ggdendro)
library(ggpubr)
library(patchwork)
library(sva)

## Defining our still sterrible hknn function for imputation - it's just very slow; one day, one day I'll come back for you. # K - nearest neighbors imputation

hknn<-function(dat, k){
  
  # Create a copy of the data, NA values to be filled in later
  dat.imp<-dat
  
  # Calculate similarity metrics for all column pairs (default is Euclidean distance)
  dist.mat<-as.matrix( dist(t(dat)) )
  #dist.mat<- 1-as.matrix(cor((dat), use="pairwise.complete.obs"))
  
  #dist.mat<-as.matrix(as.dist( dist.cosine(t(dat)) ))
  
  # Column names of the similarity matrix, same as data matrix
  cnames<-colnames(dist.mat)
  
  # For each column in the data... 
  for(X in cnames){
    
    # Find the distances of all other columns to that column 
    distances<-dist.mat[, X]
    
    # Reorder the distances, smallest to largest (this will reorder the column names as well)
    distances.ordered<-distances[order(distances, decreasing = F)]
    
    # Reorder the data matrix columns, smallest distance to largest from the column of interest
    # Obviously, first column will be the column of interest, column X
    dat.reordered<-dat[ , names(distances.ordered ) ]
    
    # Take the values in the column of interest
    vec<-dat[, X]
    
    # Which entries are missing and need to be imputed...
    na.index<-which( is.na(vec) )
    
    # For each of the missing entries (rows) in column X...
    for(i in na.index){
      
      # Find the most similar columns that have a non-NA value in this row
      closest.columns<-names( which( !is.na(dat.reordered[i, ])  ) )
      
      #print(length(closest.columns))
      
      # If there are more than k such columns, take the first k most similar
      if( length(closest.columns)>k ){
        
        # Replace NA in column X with the mean the same row in k of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns[1:k] ] )
        
      }
      
      
      # If there are less that or equal to k columns, take all the columns
      if( length(closest.columns)<=k ){
        
        # Replace NA in column X with the mean the same row in all of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns ] )
        
      }
      
      
    }
    
    # Populate a the matrix with the new, imputed values
    dat.imp[,X]<-vec
    
  }
  
  return(dat.imp)
  
}


### Read in data
# We read in the peptide level matrix from whichever prep we'd like to do this analysis for, here the example is using the first and second datasets ; npops 1 and 2
npop1_peptideLevel_df <- read.table("2024_Khan.Elcheikhali_testes_rPTR/002-SingleCellMatrices/002-Protein/npop1/npop1_peptideMatrix.PTMs_NoImpNoBCByMSRun.txt", header = TRUE, stringsAsFactors = F, sep = "\t")
npop1_peptideLevelMatrix <- npop1_peptideLevel_df %>% dplyr::select(-c(pep,prot)) %>% as.matrix

npop2_peptideLevel_df <- read.table("2024_Khan.Elcheikhali_testes_rPTR/002-SingleCellMatrices/002-Protein/npop1/npop2_peptideMatrix.PTMs_NoImp.NoBC.txt", header = TRUE, stringsAsFactors = F, sep = "\t")
npop2_peptideLevelMatrix <- npop2_peptideLevel_df %>% dplyr::select(-c(pep,prot)) %>% as.matrix

### Preprocess matrices and batch correct
# impute at peptide level, correct, unimpute
colnames(npop1_peptideLevelMatrix) <- paste("One", colnames(npop1_peptideLevelMatrix), sep = "_")
colnames(npop2_peptideLevelMatrix) <- paste("Two", colnames(npop2_peptideLevelMatrix), sep = "_")

# intersect and subset to have same pepetides
intProt <- Reduce(intersect, list(rownames(npop1_peptideLevelMatrix),rownames(npop2_peptideLevelMatrix)))

testesOne_int_unimp <- npop1_peptideLevelMatrix[which(rownames(npop1_peptideLevelMatrix) %in% intProt),]
testesTwo_int_unimp <- npop2_peptideLevelMatrix[which(rownames(npop2_peptideLevelMatrix) %in% intProt),]
mergePepLevel_unimp <- as.matrix(cbind(testesOne_int_unimp, testesTwo_int_unimp))

# impute
npop1_peptideLevelMatrix_imp <- hknn(npop1_peptideLevelMatrix, 3)
npop2_peptideLevelMatrix_imp <- hknn(npop2_peptideLevelMatrix, 3)

testesOne_int <- npop1_peptideLevelMatrix_imp[which(rownames(npop1_peptideLevelMatrix_imp) %in% intProt),]
testesTwo_int <- npop2_peptideLevelMatrix_imp[which(rownames(npop2_peptideLevelMatrix_imp) %in% intProt),]
mergePepLevel <- as.matrix(cbind(testesOne_int, testesTwo_int))

# Batch correct using ComBat
mergeProt_ref <- mergePepLevel
mergeProt <- mergeProt_ref
mergeProt[is.na(mergeProt)] <- 0
mergeProt[mergeProt == -Inf] <- 0

mergeProt_combat <- reshape2::melt(mergeProt, varnames = c("Protein", "id"))
mergeProt_combat <- mergeProt_combat %>% dplyr::mutate("prep" = ifelse(grepl("One", id), "Npop1", ifelse(grepl("Two", id), "Npop2", ifelse(grepl("Three", id), "Npop3", ifelse(grepl("Four", id), "Npop4", "Npop5")))))

batch.covs_mergeProt <- mergeProt_combat$prep[match(colnames(mergeProt), mergeProt_combat$id)]
mergeProt_batch <- ComBat(mergeProt, batch=batch.covs_mergeProt)
mergeProt_batch_ref <- mergeProt_batch

mergeProt_batch[is.na(mergePepLevel_unimp)] <- NA

#write.table(mergeProt_batch, "2024_Khan.Elcheikhali_testes_rPTR/002-SingleCellMatrices/002-Protein/Protein_npop1And2_peptideMatrix.PTMs_NoImpBCByPrep.txt", sep = "\t", row.names = TRUE)


# Read in LIGER assigned celltypes and assign/process
nPop_cellIds_all <- read.table(file = "2024_Khan.Elcheikhali_testes_rPTR/002-SingleCellMatrices/003-alignmentOutputs/Protein_cellTypeLabels_postAlignment.txt", sep = "\t", header = TRUE, stringsAsFactors = F)
allSSCCells_oneCT <- nPop_cellIds_all %>% filter(cellType %in% c("SPG","SPC","St"), Dataset == "testes1") %>% pull(id) %>% gsub("One_", "", .)

# uniprot DB for conversion
uniProtDB <-  read.fasta(file = "2024_Khan.Elcheikhali_testes_rPTR/001-MSData/002-auxiliaryFiles/004-varModSearches/uniprot-proteome_UP000005640.fasta")

### Preparing convenience maps (pep to prot, prot to gene etc)
## Prot To Gene Map
grabHeads <- getAnnot(uniProtDB)
unlistedHeads <- unlist(lapply(grabHeads, `[[`, 1))
splitHeads <- strsplit(unlistedHeads, " +" )
sp_noGN <- splitHeads[-which(grepl(pattern="GN=.*", splitHeads))]
sp_noGN_upIds <- unlist(lapply(strsplit(unlist(lapply(sp_noGN, `[[`, 1)), "\\|" ), `[[`, 2))
splitHeads <- splitHeads[which(grepl(pattern="GN=.*", splitHeads))]

# Assign to variable
Prot_to_Gene_map <- data.frame(Protein = unlist(lapply(strsplit(unlist(lapply(splitHeads, `[[`, 1)), "\\|" ), `[[`, 2)), Gene = unlist(lapply(splitHeads, grep, pattern="GN=.*", value=TRUE)))
Prot_to_Gene_map$Gene <-gsub("GN=","", Prot_to_Gene_map$Gene)

# Esepcially for shotgun DDA, this tends to be dataset specific
pep_ProtMap1 <- npop1_peptideLevel_df %>% rownames_to_column() %>% dplyr::select(c(pep,prot)) %>% distinct
pep_ProtMap2 <- npop2_peptideLevel_df %>%  rownames_to_column() %>% dplyr::select(c(pep,prot)) %>% distinct
pep_ProtMap <- rbind(pep_ProtMap1, pep_ProtMap2) %>% distinct
# List of kinases...from papers, interwebs etc 
phosphorylation_genes <- c(
  "MAP3K10", "MAP3K12", "PRKAR2A", "PPP1CC", "MAPK3", "PRKACB", "PPP2R1A",
  "CSNK1G3", "STK31", "MAP3K11", "PRKACG", "MAPK1", "PRKDC", "MTOR", "CDK13",
  "CAMK2D", "MAP4K1", "MAPK4", "PRKAR1A", "PRDX2", "PRDX4", "PKN3", "PRKAG3",
  "STK4", "DAPK3", "MAP4K4", "PIK3R1", "PIM1", "RPS6KA1", "WNK2", "GRK2", "STK16",
  "FGFR2", "FLT1", "MAP1B", "TTK", "TRIM62", "CAMK2", "MST1", "TP53BP1", "ATM",
  "MAP3K9", "CDK12", "CAMK4", "TRPM3", "PPP1R16B", "PPP2R3A", "PPP2R5A", "MAPKAPK2",
  "TSSK3", "BUB1", "CDK11B", "SMG1", "GAPDHS", "RPS6KB1"
)

# random convenience Z score function
z_score <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

## First, some preprocessing and populating of variables that are useful for thresholding and getting our main data of phosphopeptides and kinases together

# Get in gene names
t3m_GN <- peptideLevel_df[,colnames(peptideLevel_df) %in% c(allSSCCells_oneCT,"pep","prot")] %>% left_join(Prot_to_Gene_map, by = c("prot" = "Protein"))

# to ensure continuity, *sigh* ; this is an unfortunate bit - make pep column unique
t3m_GN_filt <- t3m_GN %>%
  dplyr::mutate(pep = make.unique(as.character(pep))) %>% 
  column_to_rownames(var = "pep") %>%
  dplyr::select(-c(prot, Gene)) %>%
  as.matrix()

# corelations and number of observations
t3m_GN_filt_corrs <- Hmisc::rcorr(t(t3m_GN_filt))$r %>% 
  {.[upper.tri(., diag = TRUE)] <- NA; .} %>% 
  reshape2::melt() %>% 
  mutate(pepPair = paste0(Var1,Var2)) %>% 
  left_join(pep_ProtMap, by = c("Var1" = "sequence")) %>% 
  left_join(pep_ProtMap, by = c("Var2" = "sequence")) %>% 
  left_join(Prot_to_Gene_map, by = c("protein.x" = "Protein")) %>% 
  left_join(Prot_to_Gene_map, by = c("protein.y" = "Protein")) %>% 
  distinct

t3m_GN_filt_obs <- Hmisc::rcorr(t(t3m_GN_filt))$n %>% 
  reshape2::melt() %>% 
  mutate(pepPair = paste0(Var1,Var2)) %>% 
  rename(obs = value)

t3m_GN_filt_corrs_obs <- left_join(t3m_GN_filt_corrs, t3m_GN_filt_obs[,c("pepPair", "obs")], by = 'pepPair')%>% 
  filter(obs > 30)

# specifically looking at correlations of phosphorylated peptides and thesholding
phosphoStuffs <- t3m_GN_filt_corrs_obs %>%
  filter(grepl("\\(Phospho", Var2), 
         Gene.x %in% phosphorylation_genes, 
         value > 0.35 | value < -0.35) %>%
  pull(Var2) %>%
  unique()

## So setting up our combined peptide and protein table
# so terrible, I'm sorry, idiotrepetitionville.
# First, for kinases
t3m_GN_ourPlayers_kinaseProts <- t3m_GN %>% filter(Gene %in% phosphorylation_genes) %>% pull(pep) %>% unique
t3m_GN_ourPlayers_kinases <-mergeProt_batch[rownames(mergeProt_batch) %in% t3m_GN_ourPlayers_kinaseProts,]
t3m_GN_ourPlayers_kinases_melt <- reshape2::melt(t3m_GN_ourPlayers_kinases) %>% 
  #mutate(Var2 = paste0("One_",Var2)) %>% 
  left_join(nPop_cellIds_all, by = c("Var2" = "id")) %>% 
  na.omit %>% 
  mutate(Var1 = as.character(Var1))

kinases_obsFilt <- t3m_GN_ourPlayers_kinases_melt %>% reframe(count = n(), .by = Var1) %>% filter(count > 80) %>% pull(Var1)
t3m_GN_ourPlayers_kinases_melt_filt <- t3m_GN_ourPlayers_kinases_melt %>% filter(Var1 %in% kinases_obsFilt)

# Work up to Protein level for the kinases
t3m_GN_ourPlayers_kinases_melt_filt_protLevel <- t3m_GN_ourPlayers_kinases_melt_filt %>% 
  left_join(t3m_GN[,c("pep","Gene")], by = c("Var1" = "pep")) %>% 
  reframe(Var2 = Var2, value = median(value,na.rm = T), .by = c(Gene,Var2)) %>% distinct

# Now for PhosphoPeps
t3m_GN_ourPlayers_phosphorylatedPeps <- t3m_GN %>% filter(pep %in% phosphoStuffs) %>% pull(pep) %>% unique 
t3m_GN_ourPlayers_phosPeps <-mergeProt_batch[rownames(mergeProt_batch) %in% t3m_GN_ourPlayers_phosphorylatedPeps,]
t3m_GN_ourPlayers_phosPeps_melt <- reshape2::melt(t3m_GN_ourPlayers_phosPeps) %>% 
  #mutate(Var2 = paste0("One_",Var2)) %>% 
  left_join(nPop_cellIds_all, by = c("Var2" = "id")) %>% 
  na.omit %>% 
  mutate(Var1 = as.character(Var1))

phosPeps_obsFilt <- t3m_GN_ourPlayers_phosPeps_melt %>% reframe(count = n(), .by = Var1) %>% filter(count > 37) %>% pull(Var1)

t3m_GN_ourPlayers_phosPeps_melt_filt <- t3m_GN_ourPlayers_phosPeps_melt %>% filter(Var1 %in% phosPeps_obsFilt)

# some back and forth, but, combining prot Level Kinases and phosphorylated peptides
t3m_GN_ourPlayers_melt_obsFilt_forMerge <- t3m_GN_ourPlayers_phosPeps_melt_filt %>% 
  rename(Gene = Var1) %>% 
  dplyr::select(Gene, Var2, value) %>% 
  rbind(t3m_GN_ourPlayers_kinases_melt_filt_protLevel)






### Panel b, The cluster level abundance dot plots of phosphorylated peptides and kinases
## First, we do some roundabouterrible to count and visualize numbers of pairwise observations to use for filtration
pepKinase_matrix <-reshape2::dcast(t3m_GN_ourPlayers_melt_obsFilt_forMerge, Gene ~ Var2, value.var = "value", fill=NA)

combMatrix <- pepKinase_matrix %>% 
  column_to_rownames("Gene") %>% 
  as.matrix

corrStats_combMatrix   <- rcorr(t(combMatrix))

corrStats_combMatrix_filt <- corrStats_combMatrix$n[grepl("^_", rownames(corrStats_combMatrix$n)),
                                                    !grepl("^_", colnames(corrStats_combMatrix$n))]

corrStats_combMatrix_numObs_peps <- rowSums(corrStats_combMatrix_filt)
corrStats_combMatrix_numObs_kinase <- rowSums(t(corrStats_combMatrix_filt))

par(mfrow = c(1,2))
hist(corrStats_combMatrix_numObs_peps, breaks = 20)
hist(corrStats_combMatrix_numObs_kinase, breaks = 20)

highObs_kinase <- names(corrStats_combMatrix_numObs_kinase[corrStats_combMatrix_numObs_kinase > 999])
highObs_peps <- names(corrStats_combMatrix_numObs_peps[corrStats_combMatrix_numObs_peps > 999])

## Now we can proceed with getting celltype means
abundanceDotPlot_forPVals <- t3m_GN_ourPlayers_melt_obsFilt_forMerge %>%
  left_join(nPop_cellIds_all[,c("id","cellType")], by = c("Var2" = "id")) %>% 
  filter(!cellType == "Macrophage")

abundanceDotPlot_forPVals_ref <- abundanceDotPlot_forPVals
#abundanceDotPlot_forPVals <- abundanceDotPlot_forPVals_ref

abundanceDotPlot_forPVals <- abundanceDotPlot_forPVals %>% filter(Gene %in% c(highObs_kinase, highObs_peps))

average_abundance <- abundanceDotPlot_forPVals %>%
  group_by(Gene, cellType) %>%
  summarise(avg_value = mean(value, na.rm = TRUE)) %>%
  ungroup()

## With celltype means, we head towards testing abundances
# Of course, loop life
# Create an empty list to store p-values and get vector for celltypes as they exist in DF
p_values <- list()
cell_types <- unique(abundanceDotPlot_forPVals$cellType)

# outer loop is at gene level, we set up the data structure etc
for (gene in unique(abundanceDotPlot_forPVals$Gene)) {
  
  gene_data <- abundanceDotPlot_forPVals %>% filter(Gene == gene)
  
  gene_p_values <- data.frame(
    cellType = character(),
    ttest_p = numeric(),
    wilcox_p = numeric(),
    stringsAsFactors = FALSE
  )
  # now for each cell type we test against all others
  for (cellTypes in cell_types) {
    
    data1 <- gene_data %>% filter(cellTypes == cellType) %>% pull(value)
    data2 <- gene_data %>% filter(cellTypes != cellType) %>% pull(value)
    
    # did both, parametric and non parametric
    ttest_p <- t.test(data1, data2)$p.value
    wilcox_p <- wilcox.test(data1, data2)$p.value
    
    gene_p_values <- gene_p_values %>% 
      add_row(
        cellType = cellTypes,
        ttest_p = ttest_p,
        wilcox_p = wilcox_p
      )
  }
  
  gene_p_values <- gene_p_values %>%
    mutate(
      ttest_p_adj = p.adjust(ttest_p, method = "BH"),
      wilcox_p_adj = p.adjust(wilcox_p, method = "BH")
    )
  
  p_values[[gene]] <- gene_p_values
}

# the whole list to df to list to df life
p_values_df <- bind_rows(
  lapply(names(p_values), function(gene) {
    gene_p_values <- p_values[[gene]]
    gene_p_values <- gene_p_values %>% mutate(Gene = gene)
    return(gene_p_values)
  })
)

# Merge with the average abundance data, z score and filter mistaken entries (one is not a kinase, the other is a rare  charge state peptide)
average_abundance_wPVals <- average_abundance %>%
  left_join(p_values_df, by = c("Gene", "cellType"))

average_abundance_wPVals <- average_abundance_wPVals %>%
  group_by(Gene) %>%
  mutate(avg_value_z = z_score(avg_value)) %>%
  ungroup() %>% 
  filter(!Gene %in% c("_EGILKT(Phospho (STY))AK_2", "PRDX4", "PRKACB"))

# back to wide
average_abundance_wPVals_cast <- average_abundance_wPVals %>%
  dplyr::select(Gene, cellType, avg_value_z) %>%
  spread(cellType, avg_value_z)

# Convert to matrix and set row names
average_abundance_wPVals_mat <- as.matrix(average_abundance_wPVals_cast[,-1])
rownames(average_abundance_wPVals_mat) <- average_abundance_wPVals_cast$Gene

## Work towards ordering
# hcluTs
gene_dist <- dist(average_abundance_wPVals_mat, method = "euclidean")
gene_clust <- hclust(gene_dist, method = "complete")

# extract dendogram
gene_dendro <- as.dendrogram(gene_clust)

# Extract for plotting
dendro_data <- ggdendro::dendro_data(gene_dendro)

# Order 
gene_order <- gene_clust$order
ordered_genes <- rownames(average_abundance_wPVals_mat)[gene_order]

# Reorder the original results data frame based on clustering and factor things for plotting
average_abundance_wPVals <- average_abundance_wPVals %>%
  mutate(Gene = factor(Gene, levels = ordered_genes)) %>%
  arrange(Gene)

average_abundance_wPVals <- average_abundance_wPVals %>% 
  mutate(cellType = factor(cellType, levels = c("EC","PTM","LC","SPG","SPC","St")))


## Exceptionally tedious process to programmatically label phosphosites
# Pull phosphopeptide, format and do preliminary processing
phosphoPep_genes <- t3m_GN %>% filter(pep %in% phosphoStuffs) %>% dplyr::select(pep, Gene) %>% distinct
phosphoPep_genes$baseSeq <- gsub("\\(Phospho \\(STY\\)\\)", "", gsub("[0-9_]", "", phosphoPep_genes$pep))

phosphoPep_genes <- phosphoPep_genes %>%
  mutate(phospho_aa = str_extract(pep, "\\w(?=\\(Phospho \\(STY\\)\\))")) %>%
  mutate(phospho_aa_pos = str_locate(baseSeq, phospho_aa)[,1])
phosphoPep_genes[which(phosphoPep_genes$pep == "_TMT(Phospho (STY))ISSK_2"),]$phospho_aa_pos <- 3 # Ugh, manually handled the case of tow phosphosites in 1 peptide. sopoorsopoor. 

# function to clean up and extract gene name
extract_gene_name <- function(x) {
  annot <- attr(x, "Annot")
  gene_name <- sub(".*GN=([A-Z0-9]+) .*", "\\1", annot)
  return(gene_name)
}

# First filter to swissprot only and then generate dataframe with gene name and sequence
filtered_uniProtDB <- uniProtDB[sapply(uniProtDB, function(x) grepl("^sp\\|", attr(x, "name")))]

uniProtDF <- data.frame(
  rowname = names(filtered_uniProtDB),
  gene_name = sapply(filtered_uniProtDB, extract_gene_name),
  protein_seq = sapply(filtered_uniProtDB, function(x) paste(toupper(x), collapse = ""))
)


## Ugh, we had to re-do and format all of this, couple of isoforms were missing from our fasta
# SO, first join and then make manual edits
# After that we can proceed with processing
phosphoPep_genes_interim <- phosphoPep_genes %>%
  left_join(uniProtDF, by = c("Gene" = "gene_name")) 

phosphoPep_genes_interim[phosphoPep_genes_interim$Gene == "METTL6",]$protein_seq <- "MASLQRKGLQARILTSEEEEKLKRDQTLVSDFKQQKLEQEAQKNWDLFYKRNSTNFFKDRHWTTREFEELRSCREFEDQKLTMLEAGCGVGNCLFPLLEEDPNIFAYACDFSPRAIEYVKQNPLYDTERCKVFQCDLTKDDLLDHVPPESVDVVMLIFVLSAVHPDKMHLVLQNIYKCHGCSSELRQPWDKDDFAVTWDPWSPAIRLLGEGLRHVHETLKQALYCTIFTQHLEGTDLAPALEELTSLLWCQ"
phosphoPep_genes_interim[phosphoPep_genes_interim$Gene == "METTL6",]$Gene <- "METTL6-3"

## A nice, simple join and then, get positions and format etc
phosphoPep_genes_final <- phosphoPep_genes_interim %>% 
  rowwise() %>%
  mutate(
    baseSeq_start_pos = str_locate(protein_seq, baseSeq)[1],
    PhosAminoAcidPos = baseSeq_start_pos + phospho_aa_pos - 1, 
    FinalOutput = paste(Gene, paste0(phospho_aa, PhosAminoAcidPos, "p"), sep = " ")
  ) %>%
  ungroup()


## Wonderfully roundabouterribleapproach
average_abundance_wPVals$GeneRename <- as.character(average_abundance_wPVals$Gene) 
average_abundance_wPVals$GeneRename <- ifelse(average_abundance_wPVals$GeneRename %in% phosphoPep_genes_final$pep, phosphoPep_genes_final$FinalOutput, average_abundance_wPVals$GeneRename)

ordered_genes2 <- ifelse(ordered_genes%in% phosphoPep_genes_final$pep,  phosphoPep_genes_final$FinalOutput, ordered_genes)


# Create a named vector for replacement using the FinalOutput from phosphoPep_genes_final
replacement_vector <- setNames(phosphoPep_genes_final$FinalOutput, phosphoPep_genes_final$pep)

# Replace the values using the replacement vector
average_abundance_wPVals$GeneRename <- as.character(average_abundance_wPVals$Gene)
average_abundance_wPVals$GeneRename <- replacement_vector[average_abundance_wPVals$GeneRename]

# Handle NA values by keeping original values where no match was found
average_abundance_wPVals$GeneRename[is.na(average_abundance_wPVals$GeneRename)] <- as.character(average_abundance_wPVals$Gene[is.na(average_abundance_wPVals$GeneRename)])

# Replace values in ordered_genes using the same replacement vector
ordered_genes2 <- replacement_vector[ordered_genes]
ordered_genes2[is.na(ordered_genes2)] <- ordered_genes[is.na(ordered_genes2)]

ordered_genes3 <-str_replace(ordered_genes2, "(\\S+)\\s+(.*)", "\\1 **\\2**")

# Ensure all values in GeneRename are included in ordered_genes2
all_levels <- unique(c(ordered_genes2, average_abundance_wPVals$GeneRename))

average_abundance_wPVals <- average_abundance_wPVals %>%
  mutate(GeneRename = factor(GeneRename, levels = ordered_genes2), isReplaced = GeneRename %in% Gene) %>%
  arrange(GeneRename)


average_abundance_wPVals <- average_abundance_wPVals %>% 
  mutate(GeneRename = str_replace(GeneRename, "(\\S+)\\s+(.*)", "\\1 **\\2**")) %>% 
  mutate(GeneRename = factor(GeneRename, levels = ordered_genes3), isReplaced = GeneRename %in% Gene) %>%
  arrange(GeneRename)


## Finallyfinallyfinally hallelujah, plot
# dotplot
ggplot(average_abundance_wPVals, aes(x = cellType, y = GeneRename, size = pmin(-log10(wilcox_p_adj), 5), color = avg_value_z)) +
  geom_point() +
  scale_size_continuous(name = "-log10(Qvalue)", range = c(1, 10)) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_x_discrete(expand = c(0.2, 0)) +  
  scale_y_discrete(expand = c(0.1, 0)) +  
  coord_fixed(ratio = 1) + 
  labs(
    x = "",
    y = "",
    title = "",
    color = "Z-Score of CT Mean"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 22),   
    axis.title = element_text(size = 24),  
    axis.text = element_text(size = 22),  
    legend.key.size = unit(1, 'cm'),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, margin = margin(t = -40), size = 22),
    axis.text.y = element_markdown(size = 22)
  )

# dendogram, *sigh* yes...I plotted it separately
ggplot(segment(dendro_data)) +
  geom_segment(aes(x = x, y = -y, xend = xend, yend = -yend)) +  
  coord_flip() +  
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),  
    axis.text.x = element_text(),  
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )





### Panel c, Correlation scatters for select phosphopeptide and kinase pair
# some respite/fairly straightforward barring final plot on this one
# selected pair of interest and reformat
kinaseToChex <- "TSSK3"
phosPepProt <- "METTL6-3 T218p"
phosPepProt_forPlot <- str_replace(phosPepProt, "(\\S+)\\s+(.*)", "\\1 **\\2**")
phosPepToChex <- names(replacement_vector[replacement_vector == phosPepProt])

# now, data objects. To make ggplot scatters where going split combine etc
pairOfInterest <- t3m_GN_ourPlayers_melt_obsFilt_forMerge %>% filter(Gene %in% c(kinaseToChex, phosPepToChex))

pairOfInterest_melt <- reshape2::melt(pairOfInterest) %>% 
  left_join(nPop_cellIds_all, by = c("Var2" = "id")) 

x_values <- pairOfInterest_melt %>%
  dplyr::filter(Gene != kinaseToChex) %>%
  dplyr::select(Var2, Gene, x_value = value)

y_values <- pairOfInterest_melt %>%
  dplyr::filter(Gene == kinaseToChex) %>%
  dplyr::select(Var2, y_value = value)

pairOfInterest_backTogefa <- x_values %>%
  inner_join(y_values, by = "Var2") %>%
  inner_join(pairOfInterest_melt %>% dplyr::select(Var2, cellType), by = "Var2") %>%
  distinct() %>%
  mutate(cellType = factor(cellType, levels = c("SPG", "SPC", "St")))

## Actual plots
# mek functionlol
corrPlots <- function(cell_type, data) {
  ggplot(data %>% filter(cellType == cell_type), aes(x = x_value, y = y_value)) +
    geom_point(size = 3) + 
    geom_line(stat = "smooth", method = "lm", color = "blue", se = FALSE, alpha = 0.8, linetype = "dashed", size = 2) + 
    labs(
      x = phosPepProt_forPlot,
      y = kinaseToChex,
      title = cell_type
    ) +
    theme_minimal() +
    theme(
      axis.title.y = element_text(size = 26),
      axis.text = element_text(size = 23),
      axis.title.x = element_markdown(size = 26),
      plot.title = element_text(size = 29, hjust = 0.5, face = "bold")  # Center the title
    ) + 
    ylim(-2,2) + 
    stat_cor(size = 12, label.sep = "\n", label.y = 1.5)
}

cell_types <- c("SPG", "SPC", "St")

# Generate and combine the plots using mapply
combined_plot <- mapply(function(ct) {
  corrPlots(ct, pairOfInterest_backTogefa)
}, cell_types, SIMPLIFY = FALSE)

# Add an empty plot to make the first cell blank
empty_plot <- ggplot() + theme_void()

# Combine all plots into a 2x2 grid
final_plot <- wrap_plots(empty_plot, combined_plot[[1]], combined_plot[[2]], combined_plot[[3]], ncol = 2, nrow = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(plot.margin = margin(10, 10, 10, 10)))

# Display the combined plot with square aspect ratio
final_plot & theme(aspect.ratio = 1)





### Panel d, dotplotlike with correltion differences between celltypes and their significances 
## Our first step is to compute phosPep<>kinase correlation matrices within each cell type
# redundantly assign cells to a vector for each cell type
SPG_only <- nPop_cellIds_all %>% filter(cellType %in% "SPG") %>% pull(id)
SPC_only <- nPop_cellIds_all %>% filter(cellType %in% "SPC") %>% pull(id)
St_only <- nPop_cellIds_all %>% filter(cellType %in% "St") %>% pull(id)

# listem
cell_types <- list(SPG = SPG_only, SPC = SPC_only, St = St_only)

# ready for le loups
corMat_melt <- list()
nPair_melt <- list()  

# loups to compute matrices = for each cell type, make matrix -> corr + pairwise obs -> melt both for later joins
for (cell_type in names(cell_types)) {
  
  cellType_mat <- t3m_GN_ourPlayers_melt_obsFilt_forMerge %>%
    filter(Var2 %in% cell_types[[cell_type]]) %>% 
    dplyr::select(Gene, Var2, value) %>%
    spread(key = Var2, value = value) %>%
    column_to_rownames(var = "Gene")
  
  rcorr_result <- rcorr(t(cellType_mat))
  corMat <- rcorr_result$r
  nPair <- rcorr_result$n
  
  corMat_melt[[cell_type]] <- reshape2::melt(corMat[grepl("^_", rownames(corMat)), !grepl("^_", colnames(corMat))]) %>%
    mutate(cellType = cell_type)
  
  nPair_melt[[cell_type]] <- reshape2::melt(nPair[grepl("^_", rownames(nPair)), !grepl("^_", colnames(nPair))]) %>%
    mutate(cellType = cell_type)
}

# Joins for all 3 celltypes and also correlatons + observations
allCorrs_bound <- bind_rows(corMat_melt) %>%
  mutate(Var1_Var2 = paste(Var2, Var1, sep = "_"), 
         cellType = factor(cellType, levels = c("SPG", "SPC", "St")))

allCounts_bound <- bind_rows(nPair_melt) %>%
  mutate(Var1_Var2 = paste(Var2, Var1, sep = "_"))

allCorrs_with_counts <- left_join(allCorrs_bound, allCounts_bound, by = c("Var1_Var2", "cellType"), suffix = c("_correlation", "_n")) %>% 
  filter(!value_correlation %in% c(-1,1), !is.na(value_correlation), value_n > 15, Var2_correlation != "PRKACB")


## Now, we can go ahead and carry out the Fishers Analytical test for differences in correlations
# very manual, small data so we're going to slip on by, loop is to essentially compute these correlations differences pairwise (for cell types) and get corr + nobs for each, then carry out test. 
analyticalCorrDiff_list <- list()

cellType_pairs <- list(c("SPG", "SPC"), c("SPG", "St"), c("SPC", "St"))

for (var1_var2 in unique(allCorrs_with_counts$Var1_Var2)) {
  
  kinPepPair <- dplyr::filter(allCorrs_with_counts, Var1_Var2 == var1_var2)
  
  for (pair in cellType_pairs) {
    
    cor1 <- kinPepPair$value_correlation[kinPepPair$cellType == pair[1]]
    n1 <- kinPepPair$value_n[kinPepPair$cellType == pair[1]]
    cor2 <- kinPepPair$value_correlation[kinPepPair$cellType == pair[2]]
    n2 <- kinPepPair$value_n[kinPepPair$cellType == pair[2]]
    
    if (length(cor1) > 0 && length(cor2) > 0) {
      
      r_test_result <- r.test(n = as.numeric(n1), r12 = as.numeric(cor1), 
                              n2 = as.numeric(n2), r34 = as.numeric(cor2))
      
      cor_diff <- as.numeric(cor1) - as.numeric(cor2)
      
      analyticalCorrDiff_list[[paste(var1_var2, pair[1], pair[2], sep = "_")]] <- list(
        Var1_Var2 = var1_var2,
        cellType1 = pair[1],
        cellType2 = pair[2],
        cor1 = cor1,
        cor2 = cor2,
        cor_diff = cor_diff,  
        n1 = n1,
        n2 = n2,
        z_value = r_test_result$z,
        p_value = r_test_result$p
      )
    }
  }
}

# pull together loop outputs
analyticalCorrDiff_df <- do.call(rbind, lapply(analyticalCorrDiff_list, as.data.frame))

analyticalCorrDiff_df <- analyticalCorrDiff_df %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH"))

## Plotting toimezs. Preprocess a bit to get labels etc as desired
analyticalCorrDiff_df_ref <- analyticalCorrDiff_df
#analyticalCorrDiff_df <- analyticalCorrDiff_df_ref
analyticalCorrDiff_df <- analyticalCorrDiff_df %>% mutate(Var2 = paste0(cellType1," - ",cellType2)) %>% rename(QVal = p_adjusted, Var1 = Var1_Var2)

# Use our previously worked out phosphopeptide positions, split up pair labels etc
replacement_vector2 <- replacement_vector
names(replacement_vector2) <- gsub("_[0-9]+", "", names(replacement_vector2))
names(replacement_vector2) <- gsub("_", "", names(replacement_vector2))

analyticalCorrDiff_df$pepOnly <- unlist(lapply(str_split(analyticalCorrDiff_df$Var1, "_"), '[[',3))
analyticalCorrDiff_df$protOnly <- unlist(lapply(str_split(analyticalCorrDiff_df$Var1, "_"), '[[',1))

analyticalCorrDiff_df$pepFormatted <- replacement_vector2[analyticalCorrDiff_df$pepOnly]
analyticalCorrDiff_df$pepProtComboFormat <- paste(analyticalCorrDiff_df$protOnly, analyticalCorrDiff_df$pepFormatted, sep = "_")

# Small stuff again, aesthetics and things
analyticalCorrDiff_df <- analyticalCorrDiff_df %>%
  mutate(
    pepProtComboFormat = str_replace(pepProtComboFormat, "_", "<>"),  
    Var2 = str_replace_all(Var2, "_", " - "),  
    pepProtComboFormat = str_replace(pepProtComboFormat, "(\\S+)\\s+(.*)", "\\1 **\\2**"),
    Var2_order = factor(Var2, levels = c("SPG - SPC","SPC - St","SPG - St"))
  ) %>% arrange(Var2_order)

analyticalCorrDiff_df <- analyticalCorrDiff_df %>% mutate(combObs = n1+n2)
result_df_sigPepProts <- analyticalCorrDiff_df %>% filter(QVal <= 0.15) %>% pull(Var1)
analyticalCorrDiff_df <- analyticalCorrDiff_df %>% filter(Var1 %in% result_df_sigPepProts)

analyticalCorrDiff_df <- analyticalCorrDiff_df %>%
  group_by(pepProtComboFormat) %>%
  filter(all(c("SPG - SPC", "SPC - St", "SPG - St") %in% Var2)) %>%
  ungroup()

# *phew* plot
ggplot(analyticalCorrDiff_df, aes(x = Var2_order, y = pepProtComboFormat, size = -log10(QVal), fill = cor_diff)) +
  geom_point(shape = 21) +
  scale_size_continuous(range = c(4, 20), breaks = c(0.3, 0.6,1.2)) +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = -0.25, space = "Lab", breaks = c(-0.5,0,0.5)) +
  theme_minimal() +
  scale_x_discrete(expand = c(0.4, 0)) +  
  scale_y_discrete(expand = c(0.05, 0)) +  
  coord_fixed(ratio = 1) + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 24),
    axis.text.y = element_markdown(size = 24),  
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.key.size = unit(1.25, 'cm'),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 22),
    strip.text = element_text(size = 22, face = "bold"),
    legend.justification = c(0.75, 0.5),  
    legend.margin = margin(t = 10),  
    legend.box.margin = margin(t = 10) 
  ) +
  labs(size = "-Log10(QValue)", fill = "CorrDiff") +
  guides(size = guide_legend(order = 1), fill = guide_colorbar(order = 2)) 
