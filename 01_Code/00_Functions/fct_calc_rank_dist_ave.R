# calc_rank_dist_ave function from PrimerTree, but including a "clustalo" path

if (!requireNamespace("plyr", quietly = TRUE))
  install.packages("plyr")

if (!requireNamespace("reshape2", quietly = TRUE))
  install.packages("reshape2")

if (!requireNamespace("ape", quietly = TRUE))
  install.packages("ape")


library("plyr")
library("reshape2")
library("ape")

#x <- PT.COIe
#ranks <- "genus"
#exec <- clustalo.path


calc_rank_dist_ave <- function(x, ranks = common_ranks, exec = "clustalo") {
  used_ranks <- grep("species", ranks, invert = T, value = T)
  rank_dist_mean <- data.frame(matrix(nrow = 1, ncol = 0))
  
  # Raw taxonomy data
  taxa <- as.data.frame(x$taxonomy)
  
  # Randomize the order of the taxa data frame
  taxa <- taxa[sample(nrow(taxa)), ]
  rownames(taxa) <- taxa$gi
  
  # Pick random example per species and add back in the taxa info
  unique_factors <- as.data.frame(unique(taxa$species))
  colnames(unique_factors) <- "species"
  unique_factors <- join(unique_factors, taxa, type = "left", match = "first", by = "species")
  
  # Get sequences for randomly selected species
  seqs <- x$sequence
  seqs <- seqs[names(seqs) %in% unique_factors$gi]
  
  # Align and calculate pairwise distances and convert dists to dataframe
  align <- clustalo(seqs, exec = exec)
  dists <- as.data.frame(as.matrix(dist.dna(align, model = "N", pairwise.deletion = T)))
  dists$gi <- row.names(dists)
  
  # Melt the dists dataframe so I can drop any distance that isn't within the (rank)
  melted <- melt(dists, id = "gi", variable.name = "gi2", value.name = "dist")
  
  for(rank in used_ranks) {
    
    # Gather only the needed taxa data
    unique_factors_sub <- unique_factors[ , colnames(unique_factors) %in% c("gi", "species", rank)]
    
    # Drop any row in (rank) where there is only one species represented
    # Any instance of this leads to a distance within that rank of 0, skewing the results downward
    counts <- as.data.frame(table(unique_factors_sub[[rank]]))
    colnames(counts) <- c(rank, "count")
    unique_factors_sub <- join(unique_factors_sub, counts, by = rank)
    unique_factors_sub <- unique_factors_sub[unique_factors_sub$count > 1, ]
    
    # Pull the nucleotide distance data in
    # Replace the rank1 gi with the rank1 taxa
    melted_sub <- join(melted, unique_factors_sub,by = "gi")
    melted_sub$rank1 <- as.factor(melted_sub[[rank]])
    
    # Drop all columns except the three needed so the next join doesn't get messed up
    melted_sub <- melted_sub[, colnames(melted_sub) %in% c("gi2", "dist", "rank1", "species")]
    colnames(melted_sub)[1] <- "gi"
    
    # Replace the rank2 gi with the rank2 taxa
    melted_sub <- join(melted_sub, unique_factors_sub, by = "gi")
    melted_sub$rank2 <- as.factor(melted_sub[[rank]])
    
    # Drop all columns except the three needed
    melted_sub <- melted_sub[ , colnames(melted_sub) %in% c("rank2", "dist", "rank1", "species")]
    
    # Drop all rows with missing information
    melted_sub <- na.omit(melted_sub)
    
    # We only want distances within a taxa, so drop all comparisons between taxa
    # We also want to drop any comparisons of a species to itself, which will have dist == 0
    melted_sub <- melted_sub[melted_sub$rank1 == melted_sub$rank2 & melted_sub$species != melted_sub$species.1, ]
    
    # Calculate the mean distance for each taxa compared
    #   We calculate each separately to avoid any one taxa with lots of hits (like human seqs)
    #   from skewing the mean
    melted_sub$group <- paste(melted_sub$rank1, melted_sub$rank2)
    
    # Plug the means into the storage dataframe
    rank_dist_mean[[rank]] <- mean(melted_sub$dist)
  }
  message("\nAverage number of nucleotide differences between sequences within a given taxonomic group")
  message("See function description for further details")
  rank_dist_mean
}
