##### Generation of inclusion lists for pSCoPE
### Script to generate inclusion list that will implemented via MaxQuant.Live for our Prioritized SCoPE acquisiton. 
### Broadly, consists of three sections
## 1) Processing DIANN search outputs of DIA runs to record accurate retention times
## 2) Intersecting with and processing all DDA data from previous datasets to expand inclusion list, this just shows an example of iterative process employed
## 3) Generating final inclusion list use results from pSCoPE implemented, for final inclusion list


## Load required libraries
library(stringr)
library(dplyr)
library(reshape2)
library(ggplot2)


### Define functions
## Handling DIANN files/formatting for inclusion list
# Function to munge and preprocess DIANN output to be more amenable for final inclusion list
preprocess_bulkDIANN <- function(df) {
  df$SeqCharge <- paste0(df$Stripped.Sequence, df$Precursor.Charge)
  df$modSeqProp <- gsub("\\(UniMod:730\\)", "(TMT)", df$Modified.Sequence)
  df$modSeqProp <- gsub("\\(UniMod:35\\)", "(Oxidation (M))", df$modSeqProp)
  df$modSeqPropMatch <- gsub("\\(TMT\\)", "", df$modSeqProp)
  df$modseqcharge <- paste0(df$modSeqPropMatch, df$Precursor.Charge)
  
  df$Charge_double <- as.double(df$Precursor.Charge)
  df$Mass <- ((df$Precursor.Mz) * (df$Charge_double)) - (1.0072826748 * (df$Charge_double))
  df$MassRound <- round(df$Mass, 7)
  
  return(df)
}


# Function to adjust for when our actual MS run starts/when voltage is turned on (around 25 mins, this becomes time 0 for inclusion list), convert time to seconds  
adjust_rt <- function(df) {
  df$RTAdjusted <- df$RT - 25
  df$RT_sec <- round(df$RTAdjusted * 60)
  df$ind <- rownames(df)
  df <- df %>% ungroup()
  return(df)
}


## Processing/Munging MqxQuant output
# Define a quick function for basic MQ output file processing
calc_fdr <- function(pep) {
  
  return( (cumsum(pep[order(pep)]) / seq(1, length(pep)))[order(order(pep))] )
  
}

MQ_base_FDRFilter <- function(ev, FDR) {
  
  parse_row <- grep("|", ev$Leading.razor.protein, fixed = TRUE)
  
  if (length(parse_row) > 0) {
    split_prot <- str_split(ev$Leading.razor.protein[parse_row], pattern = fixed("|"))
    split_prot2 <- unlist(split_prot)[seq(2, 3 * length(split_prot), 3)]
    ev$Leading.razor.protein[parse_row] <- split_prot2
  }
  
  ev <- ev[-which(ev$Reverse == "+"),]
  
  rev_index <- grep("REV", ev$Leading.razor.protein)
  if (length(rev_index) > 0) {
    ev <- ev[-rev_index, ]
  }
  
  con_index <- grep("CON", ev$Leading.razor.protein)
  if (length(con_index) > 0) {
    ev <- ev[-con_index, ]
  }
  
  contaminant_index <- which(ev$Potential.contaminant == "+")
  if (length(contaminant_index) > 0) {
    ev <- ev[-contaminant_index, ]
  }
  
  ev <- ev[!is.na(ev$PIF), ]
  ev <- ev[ev$PIF > 0.8, ]
  
  ev <- ev[calc_fdr(ev$PEP) < FDR, ]
  
  return(ev)
}



#### Processing DIA search results (in this case DIANN reports)
# We use the results output at 1% global FDR for prioritsation and less confident identifications for retention time alignment only (report_5perc)
# As a reminder, here, we're using our combined library made in FragPipe that combines bulk DIA data with single cell data across our first 3 datasets

DIA_report_5perc <- read.delim("2024_Khan.Elcheikhali_testes_rPTR/001-MSData/002-auxiliaryFiles/003-pSCoPE/RTGenDIARun_unfilteredDIANNReport_RTAlignOnly.tsv")
DIA_report <- read.delim("2024_Khan.Elcheikhali_testes_rPTR/001-MSData/002-auxiliaryFiles/003-SCoPE/RTGenDIA_filteredDIANNReport.tsv")


### Preprocess
## Have the right column names etc in place
DIA <- preprocess_bulkDIANN(DIA)
DIA_5perc <- preprocess_bulkDIANN(DIA_5perc)

# Adjust RT
DIA_filt_single <- adjust_rt(DIA)
DIA_unfilt_single <- DIA_5perc %>% 
  filter(!modseqcharge %in% DIA$modseqcharge) %>%   # we ensure that only less confident are used only for RT Align
  adjust_rt()

# As mentioned, unfiltered is for RT alignment only and not targetting
DIA_filt_single$TargBool <- TRUE
DIA_unfilt_single$TargBool <- FALSE

DIA_input <- rbind(DIA_filt_single, DIA_unfilt_single)

DIA_input$RTBool <- TRUE # use everything for RT alignment
DIA_input$Masses <- "376.27" # cannot be empty for MQL but isn't actually utilized...the things, you geddit
DIA_input$Leading.razor.protein <- DIA_input$Protein.Group
DIA_input$newMod <- paste0(DIA_input$Stripped.Sequence,"0",rownames(DIA_input)) # a unique identifier for each row

FinalDF <- DIA_input %>% ungroup() %>% dplyr::select(c("newMod","SeqCharge","Precursor.Charge","MassRound","Ms1.Area","Precursor.Mz","RT","RTAdjusted","RT_sec","Masses","RTBool", "TargBool"))
colnames(FinalDF) <-c("uniSeq","Seqcharge","Charge","MassRound","PeakInt","PrecMz","RT","RTAdjust","RTAdjustSec","FakeMass","RTBool","TargBool")

FinalDF$Priority <- NA # to be filled later
FinalDF$fillTime <- 512 # specifically for this work, small cells etc

# ppms etc
FinalDF$RTAdjust <- round(FinalDF$RTAdjust,3)
FinalDF$PeakInt <- round(FinalDF$PeakInt,6)
FinalDF <- as.data.frame(FinalDF)

#### Processing previous shotgun data datasets to expand inclusion list

# Read in evidence files
evOne<-as.data.frame(read_delim("2024_Khan.Elcheikhali_testes_rPTR/001-MSData/001-searchedFiles/npop1_ev_updated.txt", header = TRUE, stringsAsFactors = F)
evTwo<-as.data.frame(read_delim("2024_Khan.Elcheikhali_testes_rPTR/001-MSData/001-searchedFiles/npop2_ev_updated.txt", header = TRUE, stringsAsFactors = F))
evThree<-as.data.frame(read_delim("2024_Khan.Elcheikhali_testes_rPTR/001-MSData/001-searchedFiles/npop3_ev_updated.txt", header = TRUE, stringsAsFactors = F))
ev <- rbind(evOne, evTwo, evThree)

# Assign a variable that contains a list of finalized markers for different cell type markers we'd separately assigned using our proteomic data
Finalised_markers <- c("P09874","P31947-2","P49321","Q09028-3","Q12905","Q9UKA9","P05204","O43830","P05141","P09622","P09622-2","P18124","P41567","P42765","P50213","P50213-2","P99999","Q8TEX9","Q9UFV3","Q9UQE7","P11177-3","P14406","P20674","P21912","P22061","P84074","Q16543","Q7Z2X7","Q9UN72-2","P22314-2","P35637-2","P62081","P63313","Q15366-7","P00387-2","P35580","P62913-2","P52788","P52788-2")
Finalised_markers_DIANN <- gsub("-.*","", Finalised_markers)

# Do basic SCoPE2 preprocessing and filter
# Filtering at PIF > 0.8, feel free to change in function above etc.
ev <- MQ_base_FDRFilter(ev, FDR = 0.01)

#testiclesFromEv <- ev %>% filter(Leading.razor.protein %in% Finalised_markers) %>% distinct(modseq, .keep_all = TRUE)

### For our initial scout run, let's do the in DDA vs not in DDA separation first
# given that we iterate through, we first set confident identifications to the lowest tier so we could explore more of the space/identifiability of peptides
inDDA <- FinalDF %>% filter(Seqcharge %in% ev$seqcharge)
inDDA$Priority <- 1
notInDDA <- FinalDF %>% filter(!Seqcharge %in% ev$seqcharge & TargBool == TRUE)

### Pushing peptides identified in DIA only to see how well we can identify them in our MQL run
# Increasing our odds by pushing more intense into a higher tier

DDAIntQuant <-quantile(FinalDF$PeakInt, c(.26,.49))

notInDDA[which(notInDDA$PeakInt < DDAIntQuant[1]),]$Priority <- 2
notInDDA[which(notInDDA$PeakInt > DDAIntQuant[1] & notInDDA$PeakInt < DDAIntQuant[2]),]$Priority <- 3
notInDDA[which(notInDDA$PeakInt > DDAIntQuant[2]),]$Priority <- 4

FinalDF_act <- rbind(notInDDA, inDDA)

## So this was our first inclusion list
#write.table(FinalDF_act, file = "/Users/Saad/Desktop/testis_scoutRunDIACheck_inclusionList.tsv", sep = "\t",
#            row.names = FALSE, col.names = TRUE, quote = FALSE)

## After this, we iterated a few times over move peptides between tiers and finally settled on a final inclusion list in which we included our marker proteins in the top tier i.e tier 4 and increased the fill times for the ones that had lower data completeness. 





