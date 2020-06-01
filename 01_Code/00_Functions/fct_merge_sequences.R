# INFO --------------------------------------------------------------------

# Fonction pour regrouper les sequences d'un dossier a la fois par locus et par
# groupe taxonomique. 

# Audrey Bourret
# novembre 2018  


# merge_sequences ---------------------------------------------------------

library(tidyverse)
library(Biostrings)

merge_sequences <- function(LOCUS = c("ADNmt", "12S"), 
                            GROUP = "All", 
                            LIST = files.seq, 
                            LISTwPATH = files.seq.wPATH, 
                            PATH = path.group, 
                            FN   = LOCUS[2],
                            BUFFER = 50, 
                            REF = REF, 
                            LAB = Labo.amorces,
                            NMIS = 5,
                            KEEP = TRUE,
                            verbose = FALSE){
  
  # @LOCUS     : Doit toujours avoir ADNmt en premier, puis le locus d'intérêt ensuite, jusqu'a ce qu'on change le mode de nom de fichier .fasta  
  # @GROUP     : Doit être le nom d'une des colones de REF
  # @LIST      : Liste des fichiers fasta à traiter
  # @LISTwPATH : Chemin d'accès COMPLET de ces fichiers
  # @PATH      : Chemin d'accès des fichiers de sortie
  # @FN        : Ce qu'on veut comme info dans le fichier de sortie, par default le nom du second locus
  # @BUFFER    : Buffer autour de la zone à couper
  # @REF       : Fichier des references taxonomiques
  # @LAB       : Fichier avec les informations sur les amorces
  # @NMIS      : Nombre max de mistmatches acceptes entre les amorces et les sequences
  # @KEEP      : TRUE garde les sequences complete quand l'amorce ne colle pas, FALSE l'enleve
  # @verbose   : TRUE pour voir le nom complet de tous les fichiers traités
  
  # Traiter les amorces
  
  amorcesF <- LAB %>% filter(Locus %in% LOCUS, !is.na(Sequence.F)) %>% select(Sequence.F)  
  amorcesR <- LAB %>% filter(Locus %in% LOCUS, !is.na(Sequence.R)) %>% select(Sequence.R)      
  
  if(nrow(amorcesF) ==0 | nrow(amorcesR) ==0) {
    stop(call. = FALSE, "Aucune amorce reconnue pour ce locus!!!")
  }
  
  # Creer la liste de fichier reduite au locus
  LIST.red <- vector()
  for(x in 1:length(LOCUS)){
    LIST.red <- c(LIST.red, files.seq[str_detect(files.seq, pattern = LOCUS[x])])
  }
  
  
  LIST.red.wPATH <- vector()
  for(x in 1:length(LOCUS)){
    LIST.red.wPATH  <- c(LIST.red.wPATH , files.seq.wPATH[str_detect(files.seq.wPATH, pattern = LOCUS[x])])
  }  
  
  # Compter les sequences traites
  N.total <- length(LIST.red)
  N.int   <- 0
  
  REF$All <-  "All"
  
  # Faire la liste des groupes
  DF <- data.frame(SP =  unique(sapply(strsplit( LIST.red, "_"), `[`, 1))) %>% 
    left_join(REF %>% select(Species, GROUP),
              by = c("SP" = "Species")) 
  
  SP.reject <- DF[which(is.na(DF[,GROUP])), "SP"]

  DF <- DF[which(!is.na(DF[,GROUP])),]
  
  cat(paste0(N.total," fichiers FASTA ont ete trouves \n",length(unique(DF[[GROUP]]))," groupe(s) d'interets ont ete trouves\n"))
  
  # Faire les fichiers par groupes
  
  for(x in unique(DF[[GROUP]])){
    # Trouver les especes
    SP.GROUP <- DF[which(DF[[GROUP]] == x), "SP"]
    
    # Creer un fichier fasta par groupe
    
    fn <- paste0(x, "_", FN ,".fasta")
    file.create(file.path(PATH,fn))
    
    # Trouver les fichiers
    LIST.group <- vector()
    for(s in 1:length(SP.GROUP)){
      LIST.group <- c(LIST.group, LIST.red.wPATH[str_detect(LIST.red.wPATH, pattern = SP.GROUP[s])])
    }    
    
    # Compter
    N.int <- N.int + length(LIST.group)
    
    cat(paste0("\n",x, ": ", length(SP.GROUP), " espece(s) (",N.int,"/",N.total," sequences)- ", fn,"\n"))
    
    # Prendre les sequences une a la fois, dans chacun des groupes, puis les couper
    
   # Partir un progress bar 
     pb <- txtProgressBar(min = 0, max = length(LIST.group), style = 3)
     
     for(i in 1:length(LIST.group)){
     
      y <-  LIST.group[i]
        
      # Pour imprimer le nom complet du fichier traité s'il y a des problèmes
      if(verbose == TRUE){
        cat("\n",y)
      }

      # update progress bar
      setTxtProgressBar(pb, i)
    

      
      
      
      SEQ <- readDNAStringSet(file.path(y))
      
      min.F <- vector(length = nrow(amorcesF))
      min.R <- vector(length = nrow(amorcesF))
      
      # Boucle pour l'amorce F
      for (l in 1:length(min.F)){
        aF <- amorcesF[l,] %>% pull
        info.F <- vmatchPattern(pattern = aF, subject = SEQ,  max.mismatch = NMIS)
        
        # Boucle pour l'argument KEEP
        if(KEEP == TRUE){
          min.F[l] <- ifelse(is.null(info.F@ends[[1]]), 0, (info.F@ends[[1]]-info.F@width0+1))
        } else if(KEEP ==  FALSE){
          min.F[l] <- ifelse(is.null(info.F@ends[[1]]), NaN, (info.F@ends[[1]]-info.F@width0+1))
        }
        
      } # Fin de la boucle pour l'amorce F
      
      # Boucle pour l'amorce R      
      for (l in 1:length(min.R)){
        aR <- amorcesR[l,] %>% pull
        info.R <- vmatchPattern(pattern = as.character(reverseComplement(DNAStringSet(aR))), subject = SEQ,  max.mismatch = NMIS)
        
        # Boucle pour l'argument KEEP
        if(KEEP == TRUE){
          min.R[l] <- ifelse(is.null(info.R@ends[[1]]), SEQ@ranges@width, info.R@ends[[1]])       
        } else if(KEEP ==  FALSE){
          min.R[l] <- ifelse(is.null(info.R@ends[[1]]), NaN, info.R@ends[[1]])       
        }
        
      } # Fin de la boucle pour l'amorce R
      
      # Trouver le min et max absolu
      min.all <- suppressWarnings(min(min.F, min.R, na.rm = TRUE)) # on peut enlever le warning car donne -Inf 
      max.all <- suppressWarnings(max(min.F, min.R, na.rm = TRUE))
      
      # Si il n'y a pas de min ou max, on n'ecrit pas la sequence
      if(min.all == Inf | max.all == Inf | min.all == max.all) {
        if(verbose == TRUE){
          cat(paste0("Sequence non intetegree :", y,"\n"))
        }
      } else { # Si aucun probleme on ecrit la sequence
        
        # Corriger les min/max selon la sequence
        start <- ifelse((min.all - BUFFER)<=0, 0, (min.all - BUFFER))
        end   <- ifelse((max.all + BUFFER)>= SEQ@ranges@width, SEQ@ranges@width, (max.all + BUFFER))      
        
        # Couper la sequence
        SEQ.FINAL <- DNAStringSet(SEQ[[1]][start:end])
        names(SEQ.FINAL) <- names(SEQ)
        
        # Ecrire le fichier
        writeXStringSet(SEQ.FINAL, file.path(PATH,fn), 
                        append=TRUE, format="fasta")
      } # Fin de l'ecriture ou non d'une sequence
      
    } # Fin du traitement d'une sequence 
    
     close(pb)
  } # Fin d'un groupe
  
} # Fin de la fonction
