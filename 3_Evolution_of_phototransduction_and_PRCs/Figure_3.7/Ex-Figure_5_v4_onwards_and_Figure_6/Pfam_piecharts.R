#### Pfam domains of the TFs shared across PRCs of different species ####



setwd("/data/evassvis/aa1176/Scratch_Backup_and_substitute_from_2022-04-25/scRNAseq_analysis_working_dir/VennDiagrams_species-specific_TFs/")


# install.packages("ggplot2")
# install.packages("ggrepel")
# install.packages("tidyverse")
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

# Pie chart of how many OGs have Pfam domain vs those that don't.

# Import data file
Pfam_perc <- read.csv("/data/evassvis/aa1176/Scratch_Backup_and_substitute_from_2022-04-25/scRNAseq_analysis_working_dir/VennDiagrams_species-specific_TFs/Pfam_vs_no-pfam.csv")
# Get the positions
Pfam_perc2 <- Pfam_perc %>% 
  mutate(csum = rev(cumsum(rev(Perc))), 
         pos = Perc/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Perc/2, pos))

# Make the plot
Pfam_perc_pie <- ggplot(Pfam_perc, aes(x = "" , y = Perc, fill = fct_inorder(Category))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Purples") +
  geom_label_repel(data = Pfam_perc2,
                   aes(y = pos, label = paste0(Perc, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "")) +
  theme_void()



## Pie chart for pfam codes in all OGs.

# Import data file
Pfam_all <- read.csv("/data/evassvis/aa1176/Scratch_Backup_and_substitute_from_2022-04-25/scRNAseq_analysis_working_dir/VennDiagrams_species-specific_TFs/Pfam_all.csv")

# Get the positions
Pfam_all2 <- Pfam_all %>% 
  mutate(csum = rev(cumsum(rev(Perc))), 
         pos = Perc/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Perc/2, pos))

# Make the plot
Pfam_all_pie <- ggplot(Pfam_all, aes(x = "" , y = Perc, fill = fct_inorder(Pfam))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set3") +
  geom_label_repel(data = Pfam_all2,
                   aes(y = pos, label = paste0(Perc, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Pfam")) +
  theme_void()



## Pie chart for pfam codes in OGs present in 3 or more phyla.

# Import data file
Pfam_top <- read.csv("/data/evassvis/aa1176/Scratch_Backup_and_substitute_from_2022-04-25/scRNAseq_analysis_working_dir/VennDiagrams_species-specific_TFs/Pfam_top.csv")

# Get the positions
Pfam_top2 <- Pfam_top %>% 
  mutate(csum = rev(cumsum(rev(Perc))), 
         pos = Perc/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Perc/2, pos))

# Make the plot
Pfam_top_pie <- ggplot(Pfam_top, aes(x = "" , y = Perc, fill = fct_inorder(Pfam))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set3") +
  geom_label_repel(data = Pfam_top2,
                   aes(y = pos, label = paste0(Perc, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Pfam")) +
  theme_void()


## Pie chart for PFams in OGs present in 3 or more phyla. EXCLUDING NAs.

# Import data file
Pfam_top_no_NA <- read.csv("/data/evassvis/aa1176/Scratch_Backup_and_substitute_from_2022-04-25/scRNAseq_analysis_working_dir/VennDiagrams_species-specific_TFs/Pfam_top_no_NA.csv")

# Get the positions
Pfam_top_no_NA2 <- Pfam_top_no_NA %>% 
  mutate(csum = rev(cumsum(rev(Perc))), 
         pos = Perc/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Perc/2, pos))

# Make the plot
Pfam_top_no_NA_pie <- ggplot(Pfam_top_no_NA, aes(x = "" , y = Perc, fill = fct_inorder(Pfam))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Reds") +
  geom_label_repel(data = Pfam_top_no_NA2,
                   aes(y = pos, label = paste0(Pfam)),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Pfam")) +
  theme_void()



## Make Figure with all and top plots

Pfam_all_and_top <- ggarrange(Pfam_all_pie, Pfam_top_pie, 
                               labels = c("A", "B"), 
                               ncol = 2, nrow = 1) 

## Save Figure:
# As png:
ggsave2("/data/evassvis/aa1176/Scratch_Backup_and_substitute_from_2022-04-25/scRNAseq_analysis_working_dir/VennDiagrams_species-specific_TFs/Fig_Pfam_pies.png", plot = Pfam_all_and_top, height = 15, width = 45, units = c("cm"), dpi = 600)
# As pdf:
ggsave2("/data/evassvis/aa1176/Scratch_Backup_and_substitute_from_2022-04-25/scRNAseq_analysis_working_dir/VennDiagrams_species-specific_TFs/Fig_Pfam_pies.pdf", plot = Pfam_all_and_top, height = 15, width = 45, units = c("cm"), dpi = 600)






#### ATFDB Families of the TFs shared across PRCs of different species ####


# Pie chart of how many OGs have ATFDB Families vs those that don't.

# Import data file
ATFDB_perc <- read.csv("/data/evassvis/aa1176/Scratch_Backup_and_substitute_from_2022-04-25/scRNAseq_analysis_working_dir/VennDiagrams_species-specific_TFs/ATFDB_vs_no-ATFDB.csv")
# Get the positions
ATFDB_perc2 <- ATFDB_perc %>% 
  mutate(csum = rev(cumsum(rev(Perc))), 
         pos = Perc/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Perc/2, pos))

# Make the plot
ATFDB_perc_pie <- ggplot(ATFDB_perc, aes(x = "" , y = Perc, fill = fct_inorder(Category))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Reds") +
  geom_label_repel(data = ATFDB_perc2,
                   aes(y = pos, label = paste0(Perc, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "")) +
  theme_void()



## Pie chart for ATFDB Families in OGs present in 3 or more phyla.
## Note: Data for ATFDB pies should be updated!!!


# Import data file
ATFDB_top <- read.csv("/data/evassvis/aa1176/Scratch_Backup_and_substitute_from_2022-04-25/scRNAseq_analysis_working_dir/VennDiagrams_species-specific_TFs/ATFDB_top.csv")

# Get the positions
ATFDB_top2 <- ATFDB_top %>% 
  mutate(csum = rev(cumsum(rev(Perc))), 
         pos = Perc/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Perc/2, pos))

# Make the plot
ATFDB_top_pie <- ggplot(ATFDB_top, aes(x = "" , y = Perc, fill = fct_inorder(ATFDB))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set3") +
  geom_label_repel(data = ATFDB_top2,
                   aes(y = pos, label = paste0(Perc, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "ATFDB")) +
  theme_void()




## Pie chart for ATFDB Families in OGs present in 3 or more phyla. EXCLUDING NAs.

# Import data file
ATFDB_top_no_NA <- read.csv("/data/evassvis/aa1176/Scratch_Backup_and_substitute_from_2022-04-25/scRNAseq_analysis_working_dir/VennDiagrams_species-specific_TFs/ATFDB_top_no_NA.csv")

# Get the positions
ATFDB_top_no_NA2 <- ATFDB_top_no_NA %>% 
  mutate(csum = rev(cumsum(rev(Perc))), 
         pos = Perc/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Perc/2, pos))

# Make the plot
ATFDB_top_no_NA_pie <- ggplot(ATFDB_top_no_NA, aes(x = "" , y = Perc, fill = fct_inorder(ATFDB))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Greens") +
  geom_label_repel(data = ATFDB_top_no_NA2,
                   aes(y = pos, label = paste0(ATFDB)),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "ATFDB")) +
  theme_void()



#####################
## Make Figure with Pfam perc pie and ATFDB perc pie

Pfam_and_ATFDB <- ggarrange(Pfam_perc_pie, ATFDB_perc_pie, 
                              labels = c("A", "B"), 
                              ncol = 2, nrow = 1) 
# Save Figure
svg(file="Pfam_and_ATFDB.svg")
Pfam_and_ATFDB
dev.off()



#####################
## Make Figure with Pfam top pie and ATFDB top pie EXCLUDING NA version.

Pfam_and_ATFDB_tops_no_NA <- ggarrange(Pfam_top_no_NA_pie, ATFDB_top_no_NA_pie, 
                            labels = c("A", "B"), 
                            ncol = 2, nrow = 1) 
# Save Figure
svg(file="Pfam_and_ATFDB_top_no_NA.svg")
Pfam_and_ATFDB_tops_no_NA
dev.off()




###################################################################

## Pie to compare how many TFs there are compared to cofactors and other genes.
## Using here data updated up to july 6 2023:

# Doing the pie for the top dataset (OGs in >3 phyla)
# Import data file
Main_data_top <- read.csv("C:/Users/ale_a/OneDrive - University of Leicester/Third_Year_Report_Paper/1_Paper/My_Working_dir/Figure_5_v4/Main_top.csv")

# Get the positions
Main_data_top2 <- Main_data_top %>% 
  mutate(csum = rev(cumsum(rev(Perc))), 
         pos = Perc/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Perc/2, pos))

# Choose colours
my_colors <- brewer.pal(9, "Set1")[c(9, 5, 2)]

# Make the plot
Main_data_top_pie <- ggplot(Main_data_top, aes(x = "" , y = Perc, fill = fct_inorder(Category))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = my_colors) +   # Manually choosing to use the colours I defined.
  geom_label_repel(data = Main_data_top2,
                   aes(y = pos, label = paste0(Category)),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Category")) +
  theme_void()



