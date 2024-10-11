#### Little mini script to make plot for Figure 4
## Honestly, this isn't the most helpul mini script as figures like the one used are not super generalizable...but, anywho. Here it be. 

### load libraries
library(tidyverse)
library(gplots)
library(reshape2)
library(Hmisc)

### Load GO testing results object, available in google drive, here: 
# the object name is: test_res
load(".../filtered_go_rptr_test.RData")

### Filter/preprocess 
# for the figure it was: respiratory chain complex containing GO terms across the Spermatogonial cell types only
GOTermsOfInterest <- test_res %>%
  filter(grepl("mitochondrial respiratory chain complex", TERM, ignore.case = TRUE)) %>%
  filter(ct %in% c("SPG", "SPC", "St")) %>%
  filter(TERM != "mitochondrial respiratory chain complex assembly") %>% 
  ungroup

### Assign letters for plotting to each of the GO terms
terms <- unique(GOTermsOfInterest$TERM)
letters <- LETTERS[1:length(terms)]
term_labels <- setNames(letters, terms)
GOTermsOfInterest$term_label <- term_labels[GOTermsOfInterest$TERM]

### Make the lovely convoluted plot
ggplot(GOTermsOfInterest, aes(x = prot_av, y = mu_av, color = r_av)) +
  geom_point(size = 10) +
  geom_text(aes(label = term_label), vjust = -1, size = 8) +
  labs(x = "Consensus Protein",
       y = "Consensus mRNA",
       color = "rPTR",
       shape = "Cell Type") +
  scale_color_gradient(low = "blue", high = "red",
                       breaks = c(-0.35, 0.25, 0.9),  
                       labels = c("-0.3", "0.25", "0.9")) +  
  theme_minimal() +
  theme(axis.title = element_text(size = 26),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22),  
        axis.text = element_text(size = 24)) +
  guides(color = guide_colorbar(title = "rPTR")) +
  xlim(-0.6, 0.6) + ylim(-0.6, 0.6) +
  ggforce::geom_mark_ellipse(aes(group = ct, label = ct), 
                             size = 0,  
                             fill = "grey70",  
                             alpha = 0.2,  
                             label.fontsize = 20) +
  annotate("text", x = -0.1, y = 0.45, label = paste(paste(letters, terms, sep = ": ", collapse = ";\n")),
           size = 8, fontface = "bold", hjust = 0.5)   +
  geom_curve(aes(x = 0.3, y = -0.6, xend = -0.45, yend = -0.15),
             curvature = 0.3,   
             color = "grey70",  
             size = 2,         
             arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  geom_text(aes(x = -0.050, y = -0.375, label = "Spermatogenesis"), 
            size = 10,  angle = 320, color = "black")

