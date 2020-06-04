# Description -------------------------------------------------------------

# Code pour créer des alignements de séquence et trouver la position 
# des amorces (préalable à certains tests des amorces)
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

# Fonction pour merger les sequences a loader
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
                LAB = Primers, 
                KEEP =  FALSE, # Pour garder ceux qui fit pas ou non
                NMIS = NMIS,
                verbose = FALSE)




# Do the aligment ---------------------------------------------------------

# Check the name of the aligment

list.files("./00_Data/01_AlignSequences/COI/")


file <- file.path("./00_Data/01_AlignSequences/COI/",
                  "Echinodermata_COI_cut.fasta")
file

# Upload the alignment as a DNAStringSet object

DNA <- readDNAStringSet(file)
DNA

# Check the length distribution
hist(DNA@ranges@width)

# Find a name for the new aligment

new.file <- str_replace(file, ".fasta", "_align.fasta") 

# Do the alignment

alig.res <- msa(DNA, method = "ClustalW")

# To see the entire range
print(alig.res, show="complete")

# Then transform it as a DNAStringSet object and save it

DNA.alig <- alig.res@unmasked
DNA.alig



writeXStringSet(DNA.alig, file.path(new.file), 
                append=FALSE, format="fasta")


start(DNA.alig)


# Rename labels in aligment -----------------------------------------------

DNA.alig <- readDNAStringSet(new.file)
new.name <- paste(sapply(str_split(names(DNA.alig), pattern = " "), "[",2),
                  sapply(str_split(names(DNA.alig), pattern = " "), "[",3),
                  sep = "_")

new.name
names(DNA.alig) <- new.name

DNA.alig

writeXStringSet(DNA.alig, file.path(new.file) %>% str_replace("cut_align", "cut_align_rename"), 
                append=FALSE, format="fasta")

# Find primer position ----------------------------------------------------

# Search for F primers

for(x in Primers$Sequence.F){
    cat("\n", x, "\n")
    cat("Sequence length = ", nchar(x), "\n")
  
  
    res <- vmatchPattern(DNAString(x), DNA.alig, max.mismatch = NMIS)
  
    print(res@ends %>% unlist() %>% table())
  
}

# This is the End position
# Start position : End - lenght + 1
# Write down this information in the Primer spreadsheet

# Search for R primers

for(x in Primers$Sequence.R){
  print(x)
  cat("Sequence length = ", nchar(x))
  
  
  res <- vmatchPattern(reverseComplement(DNAString(x)), DNA.alig, max.mismatch = NMIS)
  
  print(res@ends %>% unlist() %>% table())
  
}

# This is the End position
# Start position : End - lenght + 1
# Write down this information in the Primer spreadsheet

# And its the end!


# Graphical representation ------------------------------------------------


Primers

# Try to do something to reproduce what genious is doing ...

test <- Primers %>% select(-c(F.Stop, R.Stop)) %>% 
            pivot_longer(c(Sequence.F, Sequence.R), names_to = "Direction", values_to = "Sequence") %>% 
            pivot_longer(c(F.Start, R.Start), names_to = "temp", values_to = "Position") %>% 
            mutate(Direction = Direction %>% str_remove("Sequence."),
                   temp = temp %>% str_remove(".Start")) %>% 
            filter(Direction == temp) %>% 
            select(-temp)
  

test

# loop to create a table by nucleotide

res <- data.frame(Primers = character(),
                  Locus = character(),
                  Direction = character(),
                  Group = character(),
                  Position = numeric(),
                  Nuc = character(),
                  stringsAsFactors = F)

res

for(x in 1:nrow(test)){
  
  for(y in 1:str_length(test$Sequence[x])){
    nuc <- str_sub(test$Sequence[x], y, y)
    pos = test$Position[x] + y - 1
    vec = data.frame(Primers = test$Primer[x],
                     Locus = test$Locus[x],
                     Direction = test$Direction[x],
                     Group = test$Group[x],
                     Position = pos,
                     Nuc = nuc,
                     stringsAsFactors = F)
      
      c(test$Primer[x], test$Locus[x], test$Direction[x],test$Group[x])
    res = rbind(res, vec) 
  }

  
}

res

res %>% ggplot(aes(x=Position, y = Primers, fill = Nuc)) +
  geom_tile() +
  theme_bw()
  

