# Description -------------------------------------------------------------

# Script for running the PrimerTree package
# Audrey Bourret
# 2020-06-02


# Packages ----------------------------------------------------------------

#library("ape")
#library("dplyr")
library("ggplot2")

library("readxl")
#library(rentrez)
#library(devtools)
#install_github("mvesuviusc/primertree", ref = "fixNcbiRateLimitIssue")

library("primerTree")
#library("PrimerMiner")
library("tidyverse")

#install.packages("spider")
#library(spider)


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("msa")
#library(msa)


# Modification of a PrimerTree function to include a path to clustalo
source("01_Code/00_Functions/fct_calc_rank_dist_ave.R")
source("01_Code/00_Functions/fct_filterResults.R")

# Path --------------------------------------------------------------------

# Path to Clustal-omega (for primerTree)

clustalo.path <- 'C:\\Users\\bourreta\\Documents\\Programs\\clustal-omega-1.2.2-win64\\clustalo.exe'
clustalo.path

# Data --------------------------------------------------------------------

# Amorces cibles

Primers <- read_excel(file.path("00_Data/Primers.xlsx"))
Primers


# Step 1: PrimerTree - Initial primer screening ------------------------------------

# This script section is a modification of Bylemans et al. (2018) Ecology and Evolution 
# https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.4387

# You should get clustal-omega : http://www.clustal.org/omega/#Download

# Also with an NCBI API key, it will run faster
# Info on NCBI API key : https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/


# 1. - Query each primer pair against the NCBI database and construct a primertree object.

# Can be very long if the primer is degenerated (will query every combination)

# Loop for each primer 
for(i in 1:nrow(Primers)) {
  
  print(Primers$Primer[i])
  
  # Create a "PT ..," objects
  assign(paste("PT", Primers$Primer[i], sep = "."),
         # PrimerTree core function
         search_primer_pair(name = Primers$Primer[i], 
                            Primers$Sequence.F[i], 
                            Primers$Sequence.R[i], 
                            num_aligns = 500, 
                            #api_key = "123",
                            clustal_options=c(exec=clustalo.path)
         )
  )
} # End of the loop

# Check the file that were created - one by primer pair
ls(pattern = "PT.")


# 2. - Inspect the sequence length distribution for each primer pair and remove any sequence records with a length deviating from the
#        majority of the sequences.

seq_lengths(PT.COIe)          
PT.COIe     <- filter_seqs(PT.COIe, min_length = 650, max_length = 800, exec = clustalo.path)

seq_lengths(PT.Echino)      
PT.Echino     <- filter_seqs(PT.Echino, min_length = 650, max_length = 800, exec = clustalo.path)


# 3. - Evaluate the taxonomic coverage and the specificity of the primers within a specific taxonomic group.
#        on the NCBI nomenclature). Also evaluate the taxonomic resolution of the primers at the genus level and correct for length of the
#        barcode to allow for comparisons between primers. 

# This is was has been created for each primer

PT.COIe$taxonomy
PT.Echino$taxonomy

# Check which level of taxo and which group is pertinent

taxo.level <- "phylum"
taxo.name  <- "Echinodermata"

Primers.PT <- Primers %>% add_column(N.TAXID = NA,
                                     N.TAXID.OK = NA,
                                     Specificity = NA, 
                                     PWDistance = NA,
                                     SEQ.Length = NA,
                                     Resolution = NA,
                                     N.Order = NA,
                                     N.Family = NA,
                                     N.Genus = NA,
                                     .name_repair = "check_unique") 
Primers.PT

# Compute specificity

for(i in 1:nrow(Primers.PT)) {
  # Select a "PT ..." object
  tmp1 <- paste("PT", Primers$Primer[i], sep = ".")
  # Compute the number of unique taxID (kind of species)
  Primers.PT$N.TAXID[i] <- get(tmp1)$taxonomy %>% select(taxId) %>% 
    unique() %>% 
    nrow()
  # Compute de number of right group
  Primers.PT$N.TAXID.OK[i] <- get(tmp1)$taxonomy %>% select(Group = all_of(taxo.level), taxId) %>% 
    filter(Group == taxo.name) %>% 
    select(taxId) %>% 
    unique() %>% 
    nrow()
}

Primers.PT <- Primers.PT %>% mutate(Specificity = round((N.TAXID.OK/N.TAXID)*100, digits = 2))

# Compute rank distance, mean sequence length and resolution

for(i in 1:nrow(Primers.PT)) {
  # Select a "PT ..." object
  tmp1 <- paste("PT", Primers$Primer[i], sep = ".")
  # Rank distance
  Primers.PT$PWDistance[i]  <- as.numeric(calc_rank_dist_ave(get(tmp1), ranks = c("genus"), exec = clustalo.path) )
  # Mean sequence length
  Primers.PT$SEQ.Length[i]  <- as.integer(mean(get(tmp1)$BLAST_result$product_length))
}

Primers.PT <- Primers.PT %>% mutate(Resolution = as.numeric(PWDistance)/as.numeric(SEQ.Length))

# Compute the number of unique order, family and genus within the taxonomic group of interest

for(i in 1:nrow(Primers.PT)) {
  # Select a "PT ..." object
  tmp1 <- paste("PT", Primers$Primer[i], sep = ".")
  # Order
  Primers.PT$N.Order[i] <- get(tmp1)$taxonomy %>% select(Group = all_of(taxo.level), order) %>% 
    filter(Group == taxo.name) %>% 
    select(order) %>% 
    unique() %>% 
    nrow()
  # Family
  Primers.PT$N.Family[i] <- get(tmp1)$taxonomy %>% select(Group = all_of(taxo.level), family) %>% 
    filter(Group == taxo.name) %>% 
    select(family) %>% 
    unique() %>% 
    nrow()
  # Genus
  Primers.PT$N.Genus[i] <- get(tmp1)$taxonomy %>% select(Group = all_of(taxo.level), genus) %>% 
    filter(Group == taxo.name) %>% 
    select(genus) %>% 
    unique() %>% 
    nrow()
  
}

Primers.PT
View(Primers.PT)


# Save the results

write.csv(Primers.PT, file = file.path("02_Results/01_PrimerTree", "PrimerTree.COI.Echino.csv"), row.names = F)

#Primers.PT <- read.csv(file.path("02_Results/01_PrimerTree", "PrimerTree.COI.Echino.csv"))

# 4. - Construct bar plots for: A. The taxonomic resolution of the primers expressed as the average number of nucleotide difference for a 
#        100bp barcode fragment; B. The specificity of the primers and C. The taxonomic coverage of the primers within the Actinopterygii.


#   A. Taxonomic resolution (i.e. the average number of nucleotide differences between species for a 100bp barcode fragment)

gg.resolution <- ggplot(Primers.PT, aes(x = factor(Primer, levels = Primer), y = as.numeric(Resolution)*100)) +
  geom_bar(stat="identity", fill="gray50") +
  #geom_hline(aes(yintercept = 5), linetype = "dashed", colour = "red3") +
  #coord_cartesian(ylim = c(0,15)) +
  ggtitle("A") +
  ylab("Average no. of bp difference per 100 bp") +
  theme_classic() +
  theme(plot.title =  element_text(size = 8, colour="black", face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7, colour="black"),
        axis.text.x = element_text(size = 7, colour="black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 7, colour="black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size=0.25, colour="black"),
        axis.ticks.length = unit(0.1, "cm"))
gg.resolution


#   B. Primer specificity (i.e. the percentage of unique species of interest out of the total number of unique species recovered)

gg.specificity <- ggplot(Primers.PT, aes(x = factor(Primer, levels = Primer), y = as.numeric(Specificity))) +
  geom_bar(stat="identity", fill="gray50") +
  #geom_hline(aes(yintercept = 90), linetype = "dashed", colour = "red3") +
  #coord_cartesian(ylim = c(30,100)) +
  ggtitle("B") +
  ylab("% of unique species of interest") +
  theme_classic() +
  theme(plot.title =  element_text(size = 8, colour="black", face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7, colour="black"),
        axis.text.x = element_text(size = 7, colour="black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 7, colour="black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size=0.25, colour="black"),
        axis.ticks.length = unit(0.1, "cm"))
gg.specificity

#   C. Taxonomic coverage (i.e. no. of orders of interest for which sequences were obtained)   

gg.coverage <- ggplot(Primers.PT, aes(x = factor(Primer, levels = Primer), y = as.numeric(N.Order))) +
  geom_bar(stat="identity", fill="gray50") + 
  #geom_hline(aes(yintercept = 30), linetype = "dashed", colour = "red3") +
  #coord_cartesian(ylim = c(0, 40)) +
  ggtitle("C") +
  ylab("No. of orders") +
  theme_classic() +
  theme(plot.title =  element_text(size = 8, colour="black", face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7, colour="black"),
        axis.text.x = element_text(size = 7, colour="black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 7, colour="black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size=0.25, colour="black"),
        axis.ticks.length = unit(0.1, "cm"))
gg.coverage

ggsave("02_Results/01_PrimerTree/Figure1_PrimerTree.COI.Echino.png", grid.arrange(gg.resolution, gg.specificity, gg.coverage),
       device = "png", width = 6.68, height = 16.5, units = "cm", dpi = 500)

# Tree

for(i in 1:nrow(Primers.PT)) {
  tmp1 <- paste("PT", Primers$Primer[i], sep = ".")
  tmp2 <- paste("gg.tree", Primers.PT$Primer[i], sep = ".")
  assign(tmp2, plot(get(tmp1), taxo.level, size =0.5) +
           theme_classic() +
           labs(title = LETTERS[i], x = paste("Focal taxonomic group: ", Primers.PT$N.Order[i], " orders, ",
                                              Primers.PT$N.Family[i], " families, ", Primers.PT$N.Genus[i], " genera", sep = "")) +
           theme (plot.title = element_text(size = 8, colour = "black", face = "bold"),
                  axis.title.x = element_text(size = 7, colour="gray50"),
                  axis.title.y = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.line = element_blank(),
                  legend.position = "none"))
  }

ls(pattern = "gg.tree")

gg.class.Echino
gg.class.COIe

ggsave("02_Results/01_PrimerTree/Figure2_PrimerTree.COI.Echino.png",
       grid.arrange(gg.tree.Echino, gg.tree.COIe, 
                    nrow = 1, ncol = 2),
       device = "png", width = 19.8, height = 9.9, units = "cm", dpi = 500)



