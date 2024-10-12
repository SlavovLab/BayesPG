##### Script to pre-process dart updated variable mod searches for PTM analysis in single cells
### The script will be comprised of 3 sections:

## Part 1: Running dart updated variable mod evidence files through single cell pipeline to output tables that were used to make figures script: bayesPG_PTM_Analysis
# This script was also used on each of the 3 datasets (npop1-3) to generate the tables used for the third part (filtering on basis of observations)

## Part 2: Calculating and filtering phospho peptides using local/PTM specific FDR

## Part 3: Further filtering PTM peptides based on number of observations across single cells (this is part of fig6 in the manuscript)

############ Part 1: Run variable search results through SCoPE2 pipeline to generate tables used for PTM analysis
source("functions_parameters_speedUpdate.R")
# Reference channel number
ref_channel<-2

# Add your cell type labels, must match those used in experimental design
your_labels<-c("sc", "neg")
your_control_label<-"neg"

# Load raw data 
ev <- as.data.frame(read_delim("2024_Khan.Elcheikhali_testes_rPTR/001-MSData/001-searchedFiles/npop1_varMod_ev_updated.txt"))
colnames(ev) <- gsub(" ", ".", colnames(ev))

# Parse sp|Q00000|HUMAN_XXX into just Uniprot accession: Q00000  

parse_row<-grep("|",ev$Leading.razor.protein, fixed=T)

if(length(parse_row)>0){
  split_prot<-str_split(ev$Leading.razor.protein[parse_row], pattern = fixed("|"))
  split_prot2<-unlist(split_prot)[seq(2,3*length(split_prot),3)]
  ev$Leading.razor.protein[parse_row]<-split_prot2
}
# Load experimental design and batches
design<-as.data.frame(read_csv("2024_Khan.Elcheikhali_testes_rPTR/001-MSData/002-auxiliaryFiles/npop1_annotation.csv"))

# Create unique peptide+charge column:
ev$modseq<-paste0(ev$Modified.sequence,ev$Charge)

# Add X in front of experiment names because R doesn't like column names starting with numbers
ev$Raw.file<-paste0("X",ev$Raw.file)
design$Set<-paste0("X",design$Set)


# Which columns hold the TMT Reporter ion (RI) data
ri.index<-which(colnames(ev)%in%paste0("Reporter.intensity.",1:18))

# Make sure all runs are described in design, if not, print and remove them:
not.described<-unique(ev$Raw.file)[ !unique(ev$Raw.file) %in% paste0(design$Set) ]
ev<-ev[!ev$Raw.file%in%not.described,]

# Filter out reverse hits, contaminants, and contaminated spectra...
ev<-ev[-which(ev$Reverse=="+"),]
if(length(grep("REV", ev$Leading.razor.protein))>0){ ev<-ev[-grep("REV", ev$Leading.razor.protein),] }
if(length(grep("CON", ev$Leading.razor.protein))>0){ ev<-ev[-grep("CON", ev$Leading.razor.protein),] }
if(length(which(ev$Potential.contaminant=="+"))>0){ ev<-ev[-which(ev$Potential.contaminant=="+"),] }
ev<-ev[!is.na(ev$PIF),]
ev<-ev[ev$PIF>0.8,]

# Remove peptides that are more the 10% the intensity of the carrier in the single cell runs (only)
ev<-as.data.frame(ev)
ev$mrri<-0
ev$mrri <- rowMeans(ev[, ri.index[4:length(ri.index)]] / ev[, ri.index[1]], na.rm = T)
ev<-ev[ev$mrri < 0.1, ]

ev<-ev[calc_fdr(ev$dart_PEP)<0.01, ]


# Normalize single cell runs to normalization channel
ev<-as.data.frame(ev)
ev[, ri.index] <- ev[, ri.index] / ev[, ri.index[ref_channel]]

# Organize data into a more convenient data structure:
ev.melt_init <-reshape2::melt(ev[, c("Raw.file","modseq","Leading.razor.protein", colnames(ev)[ri.index]) ],
                              id.vars = c("Raw.file","modseq","Leading.razor.protein"));

ev.melt_init$variable <- gsub("Reporter.intensity.","RI",ev.melt_init$variable)


design_melt <- reshape2::melt(design, id.vars = c("Set"))
design_melt <- design_melt[order(design_melt$Set),]
design_melt$id <- paste0("i", 1:nrow(design_melt))

ev.melt <- left_join(ev.melt_init, design_melt, by = c("Raw.file" = "Set", "variable")) %>% select(-variable)
colnames(ev.melt)<-c("Raw.file","sequence","protein","quantitation","celltype","id")


# Grab the unique number associate to each and every cell, carrier channel, and empty channel
ev.melt$id<-as.factor(ev.melt$id)

# Remove duplicate observations of peptides from a single experiment
ev.melt <- ev.melt %>% distinct(sequence, id, .keep_all = TRUE)
ev.melt<-ev.melt[!is.na(ev.melt$protein), ]

# Create additional meta data matrices
ev.melt.uniqueID<-remove.duplicates(ev.melt,"id")
ev.melt.pep<-remove.duplicates(ev.melt, c("sequence","protein") )

# Not really necessary, but, for convenience output a dataset specific peptide to protein map. 
write.table(ev.melt.pep, "2024_Khan.Elcheikhali_testes_rPTR/001-MSData/002-auxiliaryFiles\004-varModSearches/npop1_pepToProtMap.txt", sep = "\t", row.names = TRUE)

# Create data frame of peptides x cells, populated by quantitation
ev.unmelt<-reshape2::dcast(ev.melt, sequence ~ id, value.var = "quantitation", fill=NA)

# Also create matrix of same shape
ev.matrix<-as.matrix(ev.unmelt[,-1]); row.names(ev.matrix)<-ev.unmelt$sequence

# Replace all 0s with NA
ev.matrix[ev.matrix==0]<-NA
ev.matrix[ev.matrix==Inf]<-NA
ev.matrix[ev.matrix==-Inf]<-NA

# Divide matrix into single cells (including intentional blanks) and carriers
sc_cols<-unique(ev.melt$id[(ev.melt$celltype%in%c(your_labels))])
ev.matrix.sc<-ev.matrix[, colnames(ev.matrix)%in%sc_cols]


### Output tables to get PTM observation thresholds
# load in LIGER assigned celltypes/overview table of cell IDs and process
nPop_cellIds_all <- read.table(file = "2024_Khan.Elcheikhali_testes_rPTR/002-SingleCellMatrices/003-alignmentOutputs/Protein_cellTypeLabels_postAlignment.txt", sep = "\t", header = TRUE, stringsAsFactors = F) %>% filter(Dataset == "testes1") %>%
  mutate(id = sub("^.*_", "", id)) %>%
  pull(id)

# subset matrix to previously filtered cells via npop1 and also to be for PTM containing peptides only
ev.matrix.sc_filt <- ev.matrix.sc[grepl("\\(Acetyl|\\(Methyl|\\(Phospho", rownames(ev.matrix.sc)), 
                              colnames(ev.matrix.sc) %in% nPop_cellIds_all]

numObsPerModPep <- rowSums(!is.na(ev.matrix.sc_filt))
# Convert named vector to data frame
numObsPerModPep_df <- data.frame(modSeq = names(numObsPerModPep), count = as.numeric(numObsPerModPep), row.names = NULL)

# Write data frame to a text file
write.table(numObsPerModPep_df, file = "2024_Khan.Elcheikhali_testes_rPTR/001-MSData/002-auxiliaryFiles\004-varModSearches/npop1_varModPeps_numCellObs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Filter single cells ----------------------------------------------------------------------

sc.melt<-ev.melt

xd<-as_tibble( sc.melt )
xd$protein<- gsub("-[0-9]*","",xd$protein)

xd <- xd %>% group_by(id) %>% mutate(med_per_c = median(quantitation, na.rm=T))

xd$quantitation[(xd$quantitation)==Inf]<-NA
xd$quantitation[(xd$quantitation)==0]<-NA

xd <- xd %>% mutate(across(where(is.factor), as.character))

xd1 <- xd %>%
  group_by(id) %>%
  mutate(norm_q1 = quantitation / median(quantitation, na.rm=T))


xd2 <- xd1 %>%
  group_by(sequence, Raw.file) %>%
  mutate(norm_q = quantitation / mean(norm_q1, na.rm=T))


xd3<- xd2 %>%
  filter(celltype%in%c(your_labels))


xd3 <- setDT(xd3)
xd4 <- xd3[, cvq := cv(norm_q), by = c("protein","id")]


xd5<- xd4 %>%
  group_by(protein, id) %>%
  mutate(cvn = cvna(norm_q))


xd6<- xd5 %>%
  filter(cvn > 5)

xd7<-xd6 %>% group_by(id) %>% mutate(cvm=median(cvq, na.rm=T))

xdf<-as_tibble(xd7)

# Visualize distributions and assign CV threshold

hist(unique(xdf$cvm[xdf$celltype!=your_control_label]), col=rgb(0,1,0,1/4), prob=T, breaks=50, main = "X single cells ", xlab="CV")
hist(unique(xdf$cvm[xdf$celltype==your_control_label]), col=rgb(1,0,0,1/4), prob=T, add=T, breaks=40)
abline(v=0.52)

# Filter out variable wells and controls
cvPar <- 0.52
sc_kept<-unique( xdf$id[xdf$celltype!=your_control_label & xdf$cvm < cvPar])
sc0_kept<-unique( xdf$id[xdf$celltype==your_control_label & xdf$cvm > cvPar])

# Which wells to keep
keep_these<-unique( xdf$id)

sc_total<-unique( xdf$id[xdf$celltype!=your_control_label])
sc0_total<-unique( xdf$id[xdf$celltype==your_control_label])
scrate<-round(length(sc_kept) / length(sc_total),2)*100

ev.matrix.sc.f<-ev.matrix.sc[,colnames(ev.matrix.sc)%in%sc_kept]; dim(ev.matrix.sc.f)
ev.matrix.sc.f[ev.matrix.sc.f==Inf]<-NA
ev.matrix.sc.f[ev.matrix.sc.f==-Inf]<-NA
ev.matrix.sc.f[ev.matrix.sc.f==0]<-NA

xdf$control<-"sc"
xdf$control[xdf$celltype==your_control_label]<-"ctl"


# Data transformations ----------------------------------------------------

# Perform normalizations / transformations in multiple steps with visual sanity checks:
b.t<-"FD"
xlim.t<-c(-2,2)
par(mfrow=c(3,3))

# Original data, normalized to reference channel, filtered for failed wells:
t0<-ev.matrix.sc.f

hist(c(t0), breaks=b.t, xlim=xlim.t)

# Column then row normalize by intersected median (see source functions):
t1<-medInt_cr_norm(t0)
hist(c(t1), breaks=b.t, xlim=xlim.t)

# Filter for missing data:
t2<-filt.mat.rc(t1, na.row, na.col)
hist(c(t2), breaks=b.t, xlim=xlim.t)

# Log2 transform:
t3<-log2(t2)
t3[t3==Inf]<-NA
t3[t3==-Inf]<-NA
t3[t3==0]<-NA
hist(c(t3), breaks=b.t, xlim=xlim.t)


#### Batch correct
## Details outlined in paper, but, observed a sequence of MS runs that had significantly higher intensity, batch corrected for this, cells were previously selected and we just read them in here
secondClust_PTMDart_mclust <- read.csv("2024_Khan.Elcheikhali_testes_rPTR/001-MSData/002-auxiliaryFiles\004-varModSearches/npop1_secondClusterCells_DART.PTM.csv")
secondClusterCells <- secondClust_PTMDart_mclust$x
secondClusterCells <- secondClust_PTMDart_mclust

## impute peptide level missing values
t3_imp <- hknn(t3, k.t)
sum(is.na(t3_imp))
dim(t3_imp)
t3_imp[(is.na(t3_imp))]<-0

## Batch correct using ComBat
batch.covs_secondClust_t3 <- ev.melt.uniqueID$Raw.file[match(colnames(t3_imp), secondClusterCells)]
batch.covs_secondClust_t3[!is.na(batch.covs_secondClust_t3)] <- "secondClust"
batch.covs_secondClust_t3[is.na(batch.covs_secondClust_t3)] <- "mainPop"
t3_bc <- ComBat(t3_imp, batch=batch.covs_secondClust_t3)

## Remove imputed values post batch correction
t3_bc[is.na(t2)==T] <- NA
t3 <- t3_bc


## Collapse to protein level by median:
t3m<-data.frame(t3)
t3m$pep<-rownames(t3)
t3m$prot <- ev.melt.pep$protein[match(t3m$pep, ev.melt.pep$sequence)]
t3m$prot <- gsub("-[0-9]*","",t3m$prot)
trem_preMelt <- t3m

## Writing out batch corrected matrix that contains proteins matched to peptides
write.table(trem_preMelt, "2024_Khan.Elcheikhali_testes_rPTR/001-MSData/001-searchedFiles/npop1_peptideMatrix.PTMs_NoImpNoBCByMSRun.txt", sep = "\t", row.names = TRUE)


############ Part 2: Calculating and filtering phospho peptides using local/PTM specific FDR

# Define a quick function for basic MQ output file processing
MQ_base_CalcFDROnly_dart <- function(ev) {
  
  ev <- ev[-which(ev$Reverse == "+"),]
  
  if(length(grep("REV", ev$`Leading razor protein`)) > 0) { 
    ev <- ev[-grep("REV", ev$`Leading razor protein`),] 
  }
  
  if(length(grep("CON", ev$`Leading razor protein`)) > 0) { 
    ev <- ev[-grep("CON", ev$`Leading razor protein`),] 
  }
  
  if(length(which(ev$`Potential contaminant` == "+")) > 0) { 
    ev <- ev[-which(ev$`Potential contaminant` == "+"),] 
  }
  
  ev <- ev[!is.na(ev$PIF),]
  
  ev <- ev[ev$PIF > 0.8,]
  
  ev$FDR <- calc_fdr(ev$dart_PEP)
  
  return(ev)
}


## Read in DART updated variable modification searches
npop1_allMods <- as.data.frame(read_delim("2024_Khan.Elcheikhali_testes_rPTR/001-MSData/001-searchedFiles/npop1_varMod_ev_updated.txt")) # yesyes, I know, redundancies.
npop2_allMods <- as.data.frame(read_delim("2024_Khan.Elcheikhali_testes_rPTR/001-MSData/001-searchedFiles/npop2_varMod_ev_updated.txt"))
npop3_allMods <- as.data.frame(read_delim("2024_Khan.Elcheikhali_testes_rPTR/001-MSData/001-searchedFiles/npop3_varMod_ev_updated.txt"))


### Process so we can process this all cleanly, but, looping
# List with the names of modification names outputted by MQ mapped to broad modification
modifications_list <- list(
  "Methyl" = c("2 Methyl (KR)", "Methyl (KR)", "3 Methyl (KR)"),
  "Acetyl" = c("2 Acetyl (K)", "Acetyl (K)"),
  "Phospho" = c("2 Phospho (STY)", "Phospho (STY)", "3 Phospho (STY)")
)

# simply compile datasets in list for processing ease
npop_list <- list(
  npop1 = npop1_allMods,
  npop2 = npop2_allMods,
  npop3 = npop3_allMods
)

# predefine output list
out_list <- list(Acetyl = character(), Methyl = character(), Phospho = character())

## Good ol Loopsesz ; nice and patent
for (npop in names(npop_list)) {
  npop_data <- npop_list[[npop]]
  
  for (mod in names(modifications_list)) {
    mod_group <- modifications_list[[mod]]
    
    
    unique_proteins <- npop_data %>%
      MQ_base_CalcFDROnly_dart() %>%
      filter(Modifications %in% mod_group) %>%
      arrange(dart_PEP) %>%
      mutate(ModOnlyFDR = cumsum(dart_PEP) / seq_along(dart_PEP)) %>%
      filter(ModOnlyFDR <= 0.05) %>%
      pull(`Leading razor protein`) %>%
      sub("-.*", "", .) %>%
      unique()
    
    
    out_list[[mod]] <- unique(c(out_list[[mod]], unique_proteins))
  }
}

# split data up again
acetylProts_DART <- results[["Acetyl"]]
methylProts_DART <- results[["Methyl"]]
phosphoProts_DART <- results[["Phospho"]]


write.csv(acetylProts_DART, "2024_Khan.Elcheikhali_testes_rPTR/001-MSData/002-auxiliaryFiles\004-varModSearches/acetylPepContainingProts_5PercModFDR_npop1.2.3.csv")
write.csv(methylProts_DART, "2024_Khan.Elcheikhali_testes_rPTR/001-MSData/002-auxiliaryFiles\004-varModSearches/methylPepContainingProts_5PercModFDR_npop1.2.3.csv")
write.csv(phosphoProts_DART, "2024_Khan.Elcheikhali_testes_rPTR/001-MSData/002-auxiliaryFiles\004-varModSearches/phosphoPepContainingProts_5PercModFDR_npop1.2.3.csv")


############ Part 3: Further filtering PTM peptides based on number of observations across single cells (this is part of fig6 in the manuscript)

## Read in the observations per PTM peptide by dataset tables

npop1_allSCP <- as.data.frame(read_delim("2024_Khan.Elcheikhali_testes_rPTR/001-MSData/002-auxiliaryFiles\004-varModSearches/npop1_varModPeps_numCellObs.txt"))
npop2_allSCP <- as.data.frame(read_delim("2024_Khan.Elcheikhali_testes_rPTR/001-MSData/002-auxiliaryFiles\004-varModSearches/npop1_varModPeps_numCellObs.txt"))
npop3_allSCP <- as.data.frame(read_delim("2024_Khan.Elcheikhali_testes_rPTR/001-MSData/002-auxiliaryFiles\004-varModSearches/npop1_varModPeps_numCellObs.txt"))

## Preprocess
allMods_obsCount <- rbind(npop1_allSCP, npop2_allSCP, npop3_allSCP) %>% 
  group_by(modSeq) %>% 
  dplyr::summarise(totalCount = sum(count))

allMods_obsCount <- allMods_obsCount %>%
  mutate(mod = case_when(
    grepl("\\(Acetyl", modSeq) ~ "Acetyl",
    grepl("\\(Methyl", modSeq) ~ "Methyl",
    grepl("\\(Phospho", modSeq) ~ "Phospho",
    TRUE ~ "Other"
  ))

# Create cumulative observation table grouped by modification
observation_thresholds <- sort(unique(allMods_obsCount$totalCount))

cumulative_observation_table <- observation_thresholds %>%
  purrr::map_dfr(function(threshold) {
    allMods_obsCount %>%
      group_by(mod) %>%
      dplyr::summarise(num_peps = sum(totalCount >= threshold)) %>%
      mutate(observations = threshold)
  })

# Ugh, data, bins etc so, getting the value closest to selected threshold (x_intercept_value) and get number of PTMdPeps
x_intercept_value <- 300
closest_threshold <- min(cumulative_observation_table$observations[cumulative_observation_table$observations >= x_intercept_value])
num_peps_at_closest_threshold <- cumulative_observation_table %>%
  filter(observations == closest_threshold) %>%
  summarise(num_peps = sum(num_peps))

# Get pastable and formatted text for annotations
annotation_df <- data.frame(
  x = x_intercept_value,
  y = num_peps_at_closest_threshold$num_peps,
  label = paste("# at Threshold =", num_peps_at_closest_threshold$num_peps)
)

# Plothingsz
ggplot(cumulative_observation_table, aes(x = observations, y = num_peps, color = mod)) +
  geom_line(size = 1.5) +  # Adjust line thickness here
  labs(title = "",
       x = "# Single Cells",
       y = "# Modified Peptides", 
       color = "Modification") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 27), 
        axis.text.y = element_text(size = 27),
        axis.title = element_text(size = 29), 
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=29),
        legend.text = element_text(size=27),
        legend.position = c(0.3, 0.72)) +   
  geom_vline(xintercept = x_intercept_value, linetype = "dotted", color = "black", size = 1.5) +  
  geom_text(data = annotation_df, aes(x = x, y = y, label = label),
            vjust = -1, hjust = -0.1, color = "black", size = 11) + 
  xlim(0, 500) +
  guides(color = guide_legend(override.aes = list(size = 1.5)))

