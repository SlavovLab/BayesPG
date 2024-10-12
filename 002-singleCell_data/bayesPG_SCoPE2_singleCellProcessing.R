##### SCoPE2 pipeline to process MS output for single cell sets acquired in SCoPE2 style (Datasets/npops 1-4)
### This is iteration where we implemented significant speedups to the code (so we could, for some testing, process all four datasets together, the hknn function remains a to do)
### Everything remains the same parring the CV threshold we used, I've added them all in this script



source("functions_parameters_speedUpdate.R") # found in same folder in github
# User specific

# Reference channel number (1-11, or 1-16)
ref_channel<-2

# Add your cell type labels, must match those used in experimental design
your_labels<-c("sc", "neg")
your_control_label<-"neg"

# Import ------------------------------------------------------------------

# Load raw data 
ev <- as.data.frame(read_delim("2024_Khan.Elcheikhali_testes_rPTR/001-MSData/001-searchedFiles/npop1_ev_updated.txt"))

# Parse sp|Q00000|HUMAN_XXX into just Uniprot accession: Q00000  

parse_row<-grep("|",ev$Leading.razor.protein, fixed=T)

if(length(parse_row)>0){
  split_prot<-str_split(ev$Leading.razor.protein[parse_row], pattern = fixed("|"))
  split_prot2<-unlist(split_prot)[seq(2,3*length(split_prot),3)]
  ev$Leading.razor.protein[parse_row]<-split_prot2
}
# Load experimental design 
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


# Filter by PEP or FDR: CHOOSE ONE

# ev<-ev[ev$dart_PEP<0.02, ]
# ev<-ev[ev$PEP<0.02, ]
ev<-ev[calc_fdr(ev$dart_PEP)<0.01, ]
# ev<-ev[calc_fdr(ev$PEP)<0.01, ]

# Normalize single cell runs to normalization channel
ev<-as.data.frame(ev)
ev[, ri.index] <- ev[, ri.index] / ev[, ri.index[ref_channel]]

# Organize data into a more convenient data structure:
#ev <- ev %>% filter(Modifications != "Unmodified")
ev.melt_init <-reshape2::melt(ev[, c("Raw.file","modseq","Leading.razor.protein", colnames(ev)[ri.index]) ],
                              id.vars = c("Raw.file","modseq","Leading.razor.protein",));

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

# Filter single cells ----------------------------------------------------------------------

sc.melt<-ev.melt

xd<-as_tibble( sc.melt )
xd$protein<- gsub("-[0-9]*","",xd$protein)


xd <- xd %>% group_by(id) %>% mutate(med_per_c = median(quantitation, na.rm=T))

length(unique(xd$id))

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


## NOTE. We're using mean here as oppsoed to median in OG pipeline
xd7<-xd6 %>% group_by(id) %>% mutate(cvm=median(cvq, na.rm=T))

xdf<-as_tibble(xd7)


print("Number of unique proteins used in calculation:", length(unique(xdf$protein)))

hist(unique(xdf$cvm[xdf$celltype!=your_control_label]), col=rgb(0,1,0,1/4), prob=T, breaks=50, main = "X single cells ", xlab="CV")
hist(unique(xdf$cvm[xdf$celltype==your_control_label]), col=rgb(1,0,0,1/4), prob=T, add=T, breaks=40)


# USER TUNED

hist(unique(xdf$cvm[xdf$celltype!=your_control_label]), col=rgb(0,1,0,1/4), prob=T, breaks=50, main = "X single cells ", xlab="CV")
hist(unique(xdf$cvm[xdf$celltype==your_control_label]), col=rgb(1,0,0,1/4), prob=T, add=T, breaks=40)
abline(v=0.52)
# Filter out variable wells and controls
cvPar <- 0.52 # npop1
#cvPar <- 48 # npop2
#cvPar <- 0.47 # npop3
#cvPar <- 0.48 # npop4

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

my_col3<-c( "black", "purple2")

# Plot!
ggplot(data=xdf, aes(x=cvm)) + geom_density(aes(fill=control, alpha=0.5), adjust=4) + theme_pubr() +
  scale_fill_manual(values=my_col3[c(1,2)]) +
  xlab("Quantification variability") + ylab("Density") + rremove("y.ticks") + rremove("y.text") +
  font("xylab", size=35) +
  font("x.text", size=30) +
  coord_cartesian(xlim=c(0,1))+
  annotate("text", x=0.172, y= 14, label=paste0(length(sc_kept)," single cells"), size=10, color=my_col3[c(2)])+
  annotate("text", x=0.6, y= 12, label=paste0(length(sc0_kept)," control wells"), size=10, color=my_col3[c(1)])+
  annotate("text", x=0.6, y= 14, label=paste0(length(sc_total) -length(sc_kept)," single cells"), size=10, color=my_col3[c(2)])+
  annotate("text", x=0.165, y= 12, label=paste0(length(sc0_total) - length(sc0_kept)," control wells"), size=10, color=my_col3[c(1)])+
  rremove("legend") + geom_vline(xintercept=0.47, lty=2, size=2, color="gray50")





# Data transformations ----------------------------------------------------

# Perform normalizations / transformations in multiple steps with visual sanity checks:
b.t<-"FD"
xlim.t<-c(-2,2)
par(mfrow=c(3,3))

# Original data, normalized to reference channel, filtered for failed wells:
t0<-ev.matrix.sc.f

hist(c(t0), breaks=b.t, xlim=xlim.t)

# Column then row normalize by intersected median or mean (see source functions):
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



## Collapse to protein level by median:
# Important little interlude here, so, while processing this dataset we found a sequence of mass spec runs that had a pretty large batch effect (much higher intensities), we batch corrected for it and have read the file in below, so ordinarily you wouldn't read the file in and proceed. We've outlined the MS runs (more specifically cells) in another file for the variable modification PTM evidence file: 2024_Khan.Elcheikhali_testes_rPTR/001-MSData/002-auxiliaryFiles/004-varModSearches/npop1_secondClusterCells_DART.PTM.csv

t3m<-data.frame(t3)
t3m$pep<-rownames(t3)
t3m$prot <- ev.melt.pep$protein[match(t3m$pep, ev.melt.pep$sequence)]
#t3m <- read.table("2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/002-Protein/npop1/npop1_peptideMatrix_NoImp.NoBCByMSRun.txt", header = TRUE, stringsAsFactors = F, sep = "\t")
t3m$prot <- gsub("-[0-9]*","",t3m$prot)
trem_preMelt <- t3m
t3m<-reshape2::melt(t3m, variable.names = c("pep", "prot"))
colnames(t3m) <-c("pep","prot","id","quantitation")
trem_postMelt <- t3m
t3m <- setDT(t3m)
t3m2 <- t3m[,list(qp=median(quantitation, na.rm=T)),by = c("prot","id")]

t4m<-reshape2::dcast(t3m2, prot ~ id, value.var = "qp", fill=NA)
t4<-as.matrix(t4m[,-1]); row.names(t4)<-t4m$prot
hist(c(t4), breaks=b.t, xlim=xlim.t)


# Re-column and row normalize:
t4b<-medInt_cr_norm_log(t4)
hist(c(t4b), breaks=b.t, xlim=xlim.t)

# Assign to a final variable name:
ev.matrix.sc.f.n<-t4b

#write.table(t4b, "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/002-Protein/npop1/npop1_proteinMatrix_NoImp.NoBCByMsRun.txt", sep = "\t", row.names = TRUE)


## Impute single celldata
imp.input<-ev.matrix.sc.f.n
sc.imp <- hknn(imp.input, k.t)

t5<-sc.imp
sum(is.na(sc.imp))
dim(sc.imp)

sc.imp[(is.na(sc.imp))]<-0

# Batch correction with ComBat
batch.N<-table(ev.melt.uniqueID$Raw.file[ev.melt.uniqueID$id%in%colnames(sc.imp)])
sc.imp<-sc.imp[,!colnames(sc.imp)%in%ev.melt.uniqueID$id[ev.melt.uniqueID$Raw.file%in%names(batch.N)[batch.N==1]]]

# Define the batches and model:
batch.covs <- ev.melt.uniqueID$Raw.file[match(colnames(sc.imp), ev.melt.uniqueID$id)]
mod<-data.frame(ev.melt.uniqueID$celltype[match(colnames(sc.imp), ev.melt.uniqueID$id)]); colnames(mod)<-"celltype"
mod<-model.matrix(~as.factor(celltype), data=mod)

matrix.sc.batch <- ComBat(sc.imp, batch=batch.covs)
t6<-matrix.sc.batch

# visual sanity checks post-imputation:
hist(c(t5), breaks=b.t, xlim=xlim.t)
hist(c(t6), breaks=b.t, xlim=xlim.t)

#write.table(matrix.sc.batch, "2024_Khan.Elcheikhali_testes_rPTR/002-singleCellMatrices/002-Protein/npop1/npop1_proteinMatrix_Imputed.NoBC.txt", sep = "\t", row.names = TRUE)

### We generally use the unimputed matrix for most of everything we do, so you can reassign it down here
#matrix.sc.batch <- t4b




### Getting quick dataset summary
prepSummarySetup <- t3m2[-which(t3m2$qp=="" | is.na(t3m2$qp)),]
proteinsPerCell <- as_tibble(prepSummarySetup) %>% dplyr::group_by(id) %>% dplyr::summarise(protCount = length(unique(prot)))
obsPerProtein <- as_tibble(prepSummarySetup) %>% dplyr::group_by(prot) %>% dplyr::summarise(cellCount = length(unique(id)))

# PCA ------------------------------------------------------------------------

# PC's to display:
PCx<-"PC1"
PCy<-"PC2"

# Data to use:
mat.sc.imp<-medInt_cr_norm_log(matrix.sc.batch)

# Dot product of each protein correlation vector with itself
r1<-fastCor(t(matrix.sc.batch), optBLAS = T)
rsum<-rowSums(r1^2)

# Calculate the weighted data matrix:

X.m <- mat.sc.imp
#X.m <- diag(rsum) %*%  X.m # in case you'd like to use imputed data and calculated weighted PCA
pca.imp.cor <- Hmisc::rcorr(X.m)$r
#pca.imp.cor <- fastCor(X.m, optBLAS = T) # when in rome, exceptionally fast correlation with complete matrix

# PCA
sc.pca<-eigen(pca.imp.cor)
scx<-as.data.frame(sc.pca$vectors)
colnames(scx)<-paste0("PC",1:ncol(scx))
scx$id<-colnames(pca.imp.cor)

# Percent of variance explained by each principle component
pca_var <- sc.pca$values
percent_var<- pca_var/sum(pca_var)*100
plot(1:length(percent_var), percent_var, xlab="PC", ylab="% of variance explained")

add.cols<-colnames(ev.melt)[4:7]
pca.display <- left_join(scx, ev.melt.uniqueID[,c('id',add.cols)], by = 'id')

clustToCell <- read.table(file = "2024_Khan.Elcheikhali_testes_rPTR/002-SingleCellMatrices/003-alignmentOutputs/Protein_cellTypeLabels_postAlignment.txt", sep = "\t", header = TRUE, stringsAsFactors = F) %>% 
  filter(Dataset == "testes1") %>% 
  mutate(id = gsub("One_","",id)) 


refPca.disp <- pca.display
#pca.display <- refPca.disp
pca.display <- pca.display %>% left_join(clustToCell, by = 'id')

# Display 
ggscatter(pca.display,color = 'cellType' ,x =PCx, y = PCy, size = 1, alpha=0.5) +
  xlab(paste0(PCx,"  (", round(percent_var[1],0),"%)")) +
  ylab(paste0(PCy,"  (", round(percent_var[2],0),"%)")) +
  font("ylab",size=30) +
  font("xlab",size=30) +
  font("xy.text", size=20)+
  annotate("text", x=0.05-0.02, y=-0.05, label=paste0(dim(mat.sc.imp)[1], " proteins"), size=8) +
  annotate("text", x=0.062-0.03, y=-0.04, label=paste0(dim(mat.sc.imp)[2], " cells"), size=8)
