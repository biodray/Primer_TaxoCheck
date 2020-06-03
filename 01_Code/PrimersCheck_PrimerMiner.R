# Description -------------------------------------------------------------

# Script for running the PrimerTree package
# Audrey Bourret
# 2020-06-03


# Packages ----------------------------------------------------------------

library("ggplot2")
library("readxl")
library("tidyverse")

library("PrimerMiner")

# Data --------------------------------------------------------------------


# Target primer including group and primer position

Primers <- read_excel(file.path("00_Data/Primers.xlsx"))
Primers

# Locus

Locus <- "COI"


# PrimerMiner: Evalutate each primer by group ------------------------------------------------


for(x in unique(Primers$Group)){
  align <- paste0("./00_Data/01_AlignSequences/", Locus, "/", x,"_", Locus, "_cut_align_rename.fasta")
  
  print(x)
  
  for(a in unique(Primers$Primer)){
    
    print(a)
    
    # Foward evaluation
    evaluate_primer(align, 
                    Primers %>% filter(Primer == a) %>% pull(Sequence.F) %>% unique(), 
                    Primers %>% filter(Group == x, Primer == a) %>% pull(F.Start), 
                    Primers %>% filter(Group == x, Primer == a) %>% pull(F.Stop), 
                    forward = T, gap_NA = T,
                    save=paste0("./02_Results/02_PrimerMiner/", Locus, "/Evaluation_table_",x,"_", a, "_F.csv"), 
                    mm_position ="Position_v1", adjacent=2, mm_type="Type_v1")
    
    # Reverse evaluation
    evaluate_primer(align, 
                    Primers %>% filter(Primer == a) %>% pull(Sequence.R) %>% unique(), 
                    Primers %>% filter(Group == x, Primer == a) %>% pull(R.Start), 
                    Primers %>% filter(Group == x, Primer == a) %>% pull(R.Stop),  
                    forward = F, gap_NA = T,
                    save=paste0("./02_Results/02_PrimerMiner/", Locus, "/Evaluation_table_",x,"_", a, "_R.csv"), 
                    mm_position ="Position_v1", adjacent=2, mm_type="Type_v1")
    
  }
}

# Set thresholds

Thresholds <- seq(10, 300, 10)

# Create a table to extract the output

PM.Output <- data.frame("Primer" = character(), "Group" = character(),  
                        "InputThreshold" = integer(), 
                        "OK" = integer(), "FAIL" = integer(), "MISSING" = integer(), stringsAsFactors=F)


# Loop to extract the output

for(x in unique(Primers$Group)){
  print(x)
  for(a in unique(Primers$Primer)){
    print(a)
    
    # Rewrite without duplicates - F
    temp <- read.csv(paste0("./02_Results/02_PrimerMiner/", Locus, "/Evaluation_table_",x,"_", a, "_F.csv"))
    
    temp <- temp %>% select(-X) %>% arrange(Template, desc(sum)) %>% distinct(Template, .keep_all = T)
    
    write.csv(temp, paste0("./02_Results/02_PrimerMiner/", Locus, "/Evaluation_table_",x,"_", a, "_F.csv"))
    
    # Rewrite without duplicates - R
    temp <- read.csv(paste0("./02_Results/02_PrimerMiner/", Locus, "/Evaluation_table_",x,"_", a, "_R.csv"))
    
    temp <- temp %>% select(-X) %>% arrange(Template, desc(sum)) %>% distinct(Template, .keep_all = T)
    
    write.csv(temp, paste0("./02_Results/02_PrimerMiner/", Locus, "/Evaluation_table_",x,"_", a, "_R.csv"))
    
    
    for (m in Thresholds) {
      
      # Test various threshold
      res <- primer_threshold(fw = paste0("./02_Results/02_PrimerMiner/", Locus, "/Evaluation_table_",x,"_", a, "_F.csv"),
                              rw = paste0("./02_Results/02_PrimerMiner/", Locus, "/Evaluation_table_",x,"_", a, "_R.csv"),
                              threshold = m, 
                              file = paste0("./02_Results/02_PrimerMiner/", Locus, "/Evaluation_table_",x,"_", a, "_BOTH.csv"))
      
      
      PM.Output <- rbind(PM.Output, 
                         data.frame("Primer" = a, "Group" = x,  
                                    "InputThreshold" = m, 
                                    "OK" = res[1], "FAIL" = res[2], "MISSING" = res[3], stringsAsFactors=F)
      )
      
      
      
      
    }
    
  }
}



write.csv(PM.Output, "./02_Results/02_PrimerMiner/PrimerMiner_results_COI.csv",  row.names = FALSE)

PM.Output <- read.csv("./02_Results/02_PrimerMiner/PrimerMiner_results_COI.csv")


# Graphical representation 

gg.PrimerMiner <- ggplot(PM.Output, aes(fill = InputThreshold)) +
  geom_bar(aes(x = factor(InputThreshold), y = ((OK/(OK + FAIL))*100)),
           stat = "identity", width = 1) +
  scale_fill_gradient(low = "grey90", high = "grey0", guide = "colourbar") +
  labs(x = NULL, y = "% species amplified") +
  facet_grid(Primer ~ Group) +
  theme_bw() + 
  theme(axis.title = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"),
        axis.text.x = element_blank(),
        axis.line = element_line(size = 0.25, colour = "black"),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(0, "cm"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, colour = "black", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 7, colour = "black"),
        legend.margin = margin(0, 0.1, 0, 0.1, "cm"))
gg.PrimerMiner


ggsave("./02_Results/02_PrimerMiner/Figure_PrimerMiner.COI.png", gg.PrimerMiner,
       device = "png", width = 16.5, height = 13.2, units = "cm", dpi = 500)
