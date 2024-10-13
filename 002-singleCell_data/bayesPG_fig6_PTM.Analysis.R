##### Script to carry out PTM/Phospho analysis from Figure 6 in the Manuscript
### Broadly, each section of the script has to do with making panels b,c,d
library(tidyverse)
library(seqinr)
library(ggtext)


### Read in data
# We read in the peptide level matrix from whichever prep we'd like to do this analysis for, here the example is using the first data/npop1
peptideLevel_df <- read.table("2024_Khan.Elcheikhali_testes_rPTR/001-MSData/001-searchedFiles/npop1_peptideMatrix.PTMs_NoImpNoBCByMSRun.txt", header = TRUE, stringsAsFactors = F, sep = "\t")
peptideLevelMatrix <- trem_preMelt %>% dplyr::select(-c(pep,prot)) %>% as.matrix

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
pep_ProtMap <- trem_preMelt %>% dplyr::select(c(pep,prot)) %>% unique

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

# to ensure continuity, *sigh*
t3m_GN_filt <- t3m_GN %>% 
  column_to_rownames(var = "pep") %>% 
  dplyr::select(-c(prot,Gene)) %>% 
  as.matrix

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
t3m_GN_ourPlayers_kinases <-t3[rownames(t3) %in% t3m_GN_ourPlayers_getSeqs,]
t3m_GN_ourPlayers_kinases_melt <- reshape2::melt(t3m_GN_ourPlayers_kinases) %>% 
  mutate(Var2 = paste0("One_",Var2)) %>% 
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
t3m_GN_ourPlayers_phosPeps <-t3[rownames(t3) %in% t3m_GN_ourPlayers_phosphorylatedPeps,]
t3m_GN_ourPlayers_phosPeps_melt <- reshape2::melt(t3m_GN_ourPlayers_phosPeps) %>% 
  mutate(Var2 = paste0("One_",Var2)) %>% 
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
  filter(! Gene %in% c("_EGILKT(Phospho (STY))AK_2", "PRDX4"))

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
ordered_genes <- rownames(results_matrix)[gene_order]

# Reorder the original results data frame based on clustering and factor things for plotting
average_abundance_wPVals <- average_abundance_wPVals %>%
  mutate(Gene = factor(Gene, levels = ordered_genes)) %>%
  arrange(Gene)

average_abundance_wPVals <- average_abundance_wPVals %>% 
  mutate(cellType = factor(cellType, levels = c("EC","PTM","LC","SPG","SPC","St")))


## Exceptionally tedious process to programmatically label phosphosites
# Pull phosphopeptide, format and do preliminary processing
phosphoPep_genes <- t3m_GN %>% filter(pep %in% phosphoStuffs) %>% dplyr::select(pep, Gene) %>% distinct
phosphoPep_genes[is.na(phosphoPep_genes$Gene),]$Gene <- "ELAC1" # had to manually add this as was missing
phosphoPep_genes$baseSeq <- gsub("\\(Phospho \\(STY\\)\\)", "", gsub("[0-9_]", "", phosphoPep_genes$pep))

phosphoPep_genes <- phosphoPep_genes %>%
  mutate(phospho_aa = str_extract(pep, "\\w(?=\\(Phospho \\(STY\\)\\))")) %>%
  mutate(phospho_aa_pos = str_locate(baseSeq, phospho_aa)[,1])

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


phosphoPep_genes_interim[phosphoPep_genes_interim$Gene == "SMARCD2",]$protein_seq <-  "RPGMSPGNRMPMAGLQVGPPAGSPFGAAAPLRPGMPPTMMDPFRKRLLVPQAQPPMPAQRRGLKRRKMADKVLPQRVSVKRTPGRPAQWPPWDPGACSRVSGVHGSLGF"
phosphoPep_genes_interim[phosphoPep_genes_interim$Gene == "SMARCD2",]$Gene <- "SMARCD2_part"


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
kinaseToChex <- "STK31"
phosPepProt <- "SEMA3F T522p"
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
cor_matrices_melted <- list()

# loups to compute matrices
for (cell_type in names(cell_types)) {
  
  data_matrix <- t3m_GN_ourPlayers_melt_obsFilt_forMerge %>%
    filter(Var2 %in% cell_types[[cell_type]]) %>% 
    dplyr::select(Gene, Var2, value) %>%
    spread(key = Var2, value = value) %>%
    column_to_rownames(var = "Gene")
  
  cor_matrix <- rcorr(t(data_matrix))$r
  
  cor_matrices_melted[[cell_type]] <- reshape2::melt(cor_matrix[grepl("^_", rownames(cor_matrix)), !grepl("^_", colnames(cor_matrix))]) %>%
    mutate(cellType = cell_type) 
}

allCorrs_bound <- bind_rows(cor_matrices_melted) %>%
  mutate(Var1_Var2 = paste(Var2, Var1, sep = "_"), 
         cellType = factor(cellType, levels = c("SPG", "SPC", "St")))

# Reshape and convert to matrix
cellTypeCorMat <- reshape2::dcast(allCorrs_bound, Var1_Var2 ~ cellType, value.var = "value", fill = NA) %>%
  column_to_rownames("Var1_Var2") %>%
  as.matrix()



### OK, time for significants tests
# Function to calculate Fisher's Z-transformation
fisher_z_transform <- function(cor) {
  return(0.5 * log((1 + cor) / (1 - cor)))
}

# Function to test the significance of the difference in correlations
test_correlation_difference_analytical <- function(cor1, cor2, n1, n2) {
  z1 <- fisher_z_transform(cor1)
  z2 <- fisher_z_transform(cor2)
  
  diff_z <- z1 - z2 
  empirical_diff <- cor1 - cor2 
  
  se_diff <- sqrt((1 / (n1 - 3)) + (1 / (n2 - 3)))
  
  z_score <- abs(diff_z) / se_diff
  
  p_value <- 2 * (1 - pnorm(z_score))
  
  return(list(p_value = p_value, diff_z = diff_z, empirical_diff = empirical_diff))
}

# Function to apply the significance test for each row
apply_significance_test_analytical <- function(mat, n1, n2, n3) {
  results <- apply(mat, 1, function(row) {
    cor1 <- row[1]
    cor2 <- row[2]
    cor3 <- row[3]
    
    test_12 <- test_correlation_difference_analytical(cor1, cor2, n1, n2)
    test_23 <- test_correlation_difference_analytical(cor2, cor3, n2, n3)
    test_13 <- test_correlation_difference_analytical(cor1, cor3, n1, n3)
    
    return(c(test_12$p_value, test_12$diff_z, test_12$empirical_diff,
             test_23$p_value, test_23$diff_z, test_23$empirical_diff,
             test_13$p_value, test_13$diff_z, test_13$empirical_diff))
  })
  
  results_df <- as.data.frame(t(results))
  colnames(results_df) <- c("p_value_12", "diff_z_12", "empirical_diff_12",
                            "p_value_23", "diff_z_23", "empirical_diff_23",
                            "p_value_13", "diff_z_13", "empirical_diff_13")
  
  return(results_df)
}

# For our analytical test, we set sample size for each condition to 100, there are 61 phoshosite, kinase pairs
n1 <- n2 <- n3 <- 100

# Do eeet.
significance_results_analytical <- apply_significance_test_analytical(cellTypeCorMat, n1, n2, n3)

significance_results_analytical$Q_value_12 <- p.adjust(significance_results_analytical$p_value_12, method = "BH")
significance_results_analytical$Q_value_23 <- p.adjust(significance_results_analytical$p_value_23, method = "BH")
significance_results_analytical$Q_value_13 <- p.adjust(significance_results_analytical$p_value_13, method = "BH")

# Post testing, manually assign significant out, etc
phosPepKinase_corrDiffSeqs <- c(
  "STK31__LLGEGLRHVHET(Phospho (STY))LK_3",
  "STK31__TMT(Phospho (STY))ISSK_2",
  "PRKDC__AEEDEILNRS(Phospho (STY))PR_3",
  "PRDX4__QPT(Phospho (STY))MPILK_2",
  "STK16__TMT(Phospho (STY))ISSK_2",
  "TSSK3__LLGEGLRHVHET(Phospho (STY))LK_3",
  "STK16__EGILKT(Phospho (STY))AK_3",
  "STK31__SSS(Phospho (STY))PAPADIAQTVQEDLR_3",
  "STK31__LVLRIAT(Phospho (STY))DDSK_2"
)

# our favorite munging
significance_results_analytical_subset_QVal <- significance_results_analytical[rownames(significance_results_analytical) %in% phosPepKinase_corrDiffSeqs,] %>% dplyr::select(Q_value_12, Q_value_23, Q_value_13) %>% dplyr::rename(SPG_SPC = Q_value_12, SPC_St = Q_value_23, SPG_St = Q_value_13) %>%  as.matrix() %>% reshape2::melt() %>% dplyr::rename(QVal = value)

significance_results_analytical_subset_Diffs <- significance_results_analytical[rownames(significance_results_analytical) %in% phosPepKinase_corrDiffSeqs,] %>% dplyr::select(empirical_diff_12,empirical_diff_23,empirical_diff_13) %>% dplyr::rename(SPG_SPC = empirical_diff_12, SPC_St = empirical_diff_23, SPG_St = empirical_diff_13) %>% as.matrix() %>% reshape2::melt() %>% dplyr::rename(corrDiff = value)

significance_results_analytical_subset <- left_join(significance_results_analytical_subset_QVal, significance_results_analytical_subset_Diffs, by = c("Var1" = "Var1", "Var2" = "Var2"))


## bring in proper labels
replacement_vector2 <- replacement_vector
names(replacement_vector2) <- gsub("_[0-9]+", "", names(replacement_vector2))
names(replacement_vector2) <- gsub("_", "", names(replacement_vector2))

significance_results_analytical_subset$pepOnly <- unlist(lapply(str_split(significance_results_analytical_subset$Var1, "_"), '[[',3))
significance_results_analytical_subset$protOnly <- unlist(lapply(str_split(significance_results_analytical_subset$Var1, "_"), '[[',1))

significance_results_analytical_subset$pepFormatted <- replacement_vector2[significance_results_analytical_subset$pepOnly]

significance_results_analytical_subset$pepProtComboFormat <- paste(significance_results_analytical_subset$protOnly, significance_results_analytical_subset$pepFormatted, sep = "_")


# Further plotting formatting
significance_results_analytical_subset <- significance_results_analytical_subset %>%
  mutate(
    pepProtComboFormat = str_replace(pepProtComboFormat, "_", "<>"),  
    Var2 = str_replace_all(Var2, "_", " - "),  
    pepProtComboFormat = str_replace(pepProtComboFormat, "(\\S+)\\s+(.*)", "\\1 **\\2**"),
    Var2_order = factor(Var2, levels = c("SPG - SPC","SPC - St","SPG - St"))
  ) %>% arrange(Var2_order)


# *phew* plot
ggplot(significance_results_analytical_subset, aes(x = Var2_order, y = pepProtComboFormat, size = -log10(QVal), fill = corrDiff)) +
  geom_point(shape = 21) +
  scale_size_continuous(range = c(4, 20)) +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, space = "Lab", breaks = c(-0.5,0,0.5)) +
  theme_minimal() +
  scale_x_discrete(expand = c(0.2, 0)) +  
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
    legend.justification = c(0.85, 0.5),  
    legend.margin = margin(t = 10),  
    legend.box.margin = margin(t = 10) 
  ) +
  labs(size = "-Log10(QValue)", fill = "CorrDiff") +
  guides(size = guide_legend(order = 1), fill = guide_colorbar(order = 2))


