# Description -------------------------------------------------------------

# Test taxonomic specificity and compile results
# Audrey Bourret
# 2020-06-05


# Packages ----------------------------------------------------------------

library("readxl")
library("tidyverse")

#install.packages("spider")
library("spider")
library("ape")


`%nin%` = Negate(`%in%`)


# Data --------------------------------------------------------------------


# Especes cibles

SP.target <- read_excel(file.path("00_Data/TargetSpecies.xlsx"))
SP.target

# Amorces cibles

Primers <- read_excel(file.path("00_Data/Primers.xlsx"))
Primers

# Locus

Locus <- "COI"



# Compute Taxonomic specificity ---------------------------------------------------

SP.assign <- data.frame("Primer" = character(), "Group" = character(),
                        "Species" = character(), "assign.sp" = character(), stringsAsFactors=F)


Taxo.assign <- data.frame("Primer" = character(), "Group" = character(),  
                          "Cat" = character(), 
                          "Perc.OK" = integer(), stringsAsFactors=F)

Taxo.gap <- data.frame("Primer" = character(), "Group" = character(),  
                       "Cat" = character(), 
                       "Val" = integer(), stringsAsFactors=F)


for(x in unique(Primers$Group)){
  print(x)
  
  align <- paste0("./00_Data/01_AlignSequences/", Locus, "/", x,"_", Locus, "_cut_align_rename.fasta")
  
  DNA <- read.dna(file = align,
                  format = "fasta", as.matrix = T)
  
  #dimnames(DNA)[[1]] <- data.frame(Espece_initial = dimnames(DNA)[[1]]) %>% left_join(REF %>% select(Espece_initial, Espece)) %>% 
  #                                                                          pull(Espece)
  #Vecteur pour les analyses plus tard
  vec.sp <- dimnames(DNA)[[1]]
  vec.genre <- sapply(str_split(dimnames(DNA)[[1]], pattern = "_"), `[`, 1)
  
  #vec.genre <-  data.frame(Espece_initial = dimnames(DNA)[[1]]) %>% left_join(REF %>% select(Espece_initial, Genus)) %>% 
  #                                                                   pull(Genus)
  
  for(a in unique(Primers$Primer)){
    print(a)
    
    Am.sub <- Primers %>% filter(Group == x, Primer == a) %>% select(F.Stop, R.Start)
    
    DNA.sub <- DNA[,  (Am.sub[[1]] + 1):(Am.sub[[2]] - 1)]
    
    dist.sub <- dist.dna(DNA.sub, model = "N") 
    
    # pê à faire en gardant 1 espece seulement
    bs.tab <- bestCloseMatch(dist.sub, vec.sp, threshold = 2)
    
    bs.res <- data.frame(SP = vec.sp, bs= bs.tab) %>% 
      mutate(assign.sp = ifelse(bs %in% c("correct", "no id"), 1, 0)) %>% 
      group_by(SP, assign.sp) %>% 
      summarise(N = n()) %>% 
      arrange(SP, assign.sp) %>% 
      distinct(SP, .keep_all = TRUE) %>% 
      mutate(Primer = a, Group = x) %>% 
      select(Primer, Group, Species = SP, assign.sp) %>% 
      as.data.frame()    
    
    SP.assign <- rbind(SP.assign, bs.res)
    
    
    
    sp.tab <- table(bs.tab)
    
    res <- (ifelse("correct" %in% dimnames(sp.tab)[[1]], sp.tab[["correct"]],0) + ifelse("no id" %in% dimnames(sp.tab)[[1]], sp.tab[["no id"]],0)) / sum(sp.tab)
    
    # res
    Taxo.assign <- rbind(Taxo.assign ,
                         data.frame("Primer" = a, "Group" = x,  
                                    "Cat" = "Species", 
                                    "Perc.OK" = res, stringsAsFactors=F))
    
    
    
    genre.tab <- table(bestCloseMatch(dist.sub, vec.genre, threshold = 2))
    
    res <- (ifelse("correct" %in% dimnames(genre.tab)[[1]], genre.tab[["correct"]],0) + ifelse("no id" %in% dimnames(genre.tab)[[1]], genre.tab[["no id"]],0)) / sum(genre.tab)
    
    
    
    # res
    Taxo.assign <- rbind(Taxo.assign ,
                         data.frame("Primer" = a, "Group" = x,  
                                    "Cat" = "Genus", 
                                    "Perc.OK" = res, stringsAsFactors=F))
    
    # Variation intra vs inter
    inter <- nonConDist(dist.sub, vec.sp)
    intra <- maxInDist(dist.sub, vec.sp)    
    
    
    Taxo.gap <- rbind(Taxo.gap,
                      data.frame("Primer" = a, "Group" = x,  
                                 "Cat" = c(rep("Inter", length(inter)), rep("Intra", length(intra))), 
                                 "Val" = c(inter, intra), stringsAsFactors=F))
    
  }
}


Taxo.gap %>% group_by(Primer, Group, Cat, Val) %>% 
  summarise (n = n()) %>%
  mutate(freq = n / sum(n) *100) %>% 
  ggplot(aes(x = Val, y = freq,fill = Cat)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_vline(xintercept = 2, lty = "dashed") +
  labs(x= "genetic distance (pb)", y = "% observations") +
  facet_grid(Primer ~ Group) +
  theme_bw()

Taxo.assign %>% ggplot(aes(y = Perc.OK, x = Primer, fill = Primer)) + 
  geom_bar(position = "dodge", stat = "identity") +
  facet_grid(Cat ~ Group, scale = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




# Compile info by species ----------------------------------------


SP.assign

SP.amplif <- data.frame("Primer" = character(), "Group" = character(),
                        "Species" = character(), "seuil.10" = character(), "seuil.75" = character(),
                        stringsAsFactors=F)


for(x in unique(Primers$Group)){
  print(x)
  for(a in unique(Primers$Primer)){
    print(a)
    
    # load PrimerMiner results
    
    temp1 <- read.csv(paste0("./02_Results/02_PrimerMiner/", Locus, "/Evaluation_table_",x,"_", a, "_F.csv"))
    temp2 <- read.csv(paste0("./02_Results/02_PrimerMiner/", Locus, "/Evaluation_table_",x,"_", a, "_R.csv"))
    temp3 <- read.csv(paste0("./02_Results/02_PrimerMiner/", Locus, "/Evaluation_table_",x,"_", a, "_BOTH.csv"))
    
    temp <- cbind(temp1[,c("Template", "sequ")], temp2[,c("sequ")], temp3[,c("fw", "rw", "sum")])
    names(temp) <- c("Species", "F.seq", "R.seq", "F.biais", "R.biais", "Total.bias")
    
    res <- temp %>% mutate(seuil.10 = ifelse(Total.bias <= 10, 1,0),
                           seuil.75 = ifelse(Total.bias <= 75, 1,0),
                           Primer = a,
                           Group = x) %>% 
      select(Primer, Group, Species, seuil.10, seuil.75) %>% 
      as.data.frame()
    
    
    SP.amplif <- rbind(SP.amplif, res)
    
  }
  
}



SP.amplif %>% head()

SP.assign %>% head()

Sp.decision <- SP.amplif %>% full_join(SP.assign, by = c("Primer", "Group", "Species")) %>% 
  mutate(res.10 = ifelse(seuil.10 == 0, "No amplif",
                         ifelse(assign.sp == 0, "Amplif and no resolution", "Amplif and resolution")),
         
         res.75 = ifelse(seuil.75 == 0, "No amplif",
                         ifelse(assign.sp == 0, "Amplif and no resolution", "Amplif and resolution"))
  )


graph.final <- Sp.decision %>% gather(res.10, res.75, key="seuil", value = "res") %>% 
  group_by(Primer, Group, seuil, res) %>% 
  summarise(N = n()) %>% 
  filter(!is.na(res)) %>% 
  mutate(res = factor(res, levels = rev(c("Amplif and resolution", "Amplif and no resolution", "No amplif")))) %>% 
  ggplot(aes(y=N, x=seuil, fill = res)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("darkgray", "blue1", "green3")) +
  facet_grid(Group ~ Primer, scale = "free") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0))
graph.final

ggsave("./02_Results/Figure_PrimerDecision.png", graph.final, 
       units = "in", dpi = 300, width = 7.5, height = 5)   
