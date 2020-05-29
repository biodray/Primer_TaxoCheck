# Description -------------------------------------------------------------

# Code pour créer des jeux de données pour tester les amorces
# Audrey Bourret
# 2020-05-27

# Library -----------------------------------------------------------------

library(readxl)
library(tidyverse)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("msa")

library(Biostrings)
library(msa)

# Data --------------------------------------------------------------------

# Especes cibles

SP.target <- read_excel(file.path("00_Data/TargetSpecies.xlsx"))
SP.target

# Amorces cibles

Primers <- read_excel(file.path("00_Data/Primers.xlsx"))
Primers

# Autres infos

LOCUS <- "COI" # Locus
NMIS <-  5     # Nombre de mismatches permis

# Sequences disponibles

files.seq      <- c(list.files("00_Data/00_REFSequences/"))
files.seq[1:6]

# La même liste, mais avec le chemin d'accès
files.seq.wPATH <- c(list.files("00_Data/00_REFSequences/", full.names = T))
files.seq.wPATH[1:6]

# Fonction pour merger les sequences
source(file.path("01_Code/00_Functions/fct_merge_sequences.R"))
       


# Merge sequences ---------------------------------------------------------

merge_sequences(LOCUS = c("ADNmt", "mtDNA", LOCUS), 
                GROUP = "Group", 
                LIST = files.seq, 
                LISTwPATH = files.seq.wPATH, 
                PATH = "./00_Data/01_AlignSequences/COI/", 
                FN = paste(LOCUS, "cut", sep ="_"),
                BUFFER = 0, 
                REF = SP.target,
                #REF = REF.cibles %>% filter(Groupe == "Envahissante"), 
                LAB = Primers, 
                KEEP =  FALSE, # Pour garder ceux qui fit pas
                NMIS = NMIS,
                verbose = FALSE)




list.files("./00_Data/01_AlignSequences/COI/")


file <- file.path("./00_Data/01_AlignSequences/COI/",
               "Echinodermata_COI_cut.fasta")


#primer.seq <- DNAstrin

DNA <- readDNAStringSet(file)

DNA


hist(DNA@ranges@width)



new.file <- str_replace(file, ".fasta", "_align.fasta") 

res <- msa(DNA, method = "ClustalW")

str(res)

DNA.alig <- res@unmasked



# Search for F primers

for(x in Primers$Sequence.F){
    cat("\n", x, "\n")
    cat("Sequence length = ", nchar(x), "\n")
  
  
    res <- vmatchPattern(DNAString(x), DNA.alig, max.mismatch = NMIS)
  
    print(res@ends %>% unlist() %>% table())
  
}

# This is the End position
# Start position : End - lenght + 1


# Search for R primers

for(x in Primers$Sequence.R){
  print(x)
  cat("Sequence length = ", nchar(x))
  
  
  res <- vmatchPattern(reverseComplement(DNAString(x)), DNA.alig, max.mismatch = NMIS)
  
  print(res@ends %>% unlist() %>% table())
  
}

# This is the End position
# Start position : End - lenght + 1


# So this information need to be add in your file