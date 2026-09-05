#### Supplemental_Figures_2,3,& 5####
###Created on: 3/10/26   
###Emma Thompson 
# Creating Supplemental Figure 2: Effect of time on composition in SynComs containing Bd (A) or Bsal (B) was determined via PERMANOVA on weighted UniFrac distances. Shape depicts sampling timepoint.
# Figure 3: Bacterial community composition over time subset by temperature. Points represent composition in individual SynComs containing Bd or Bsal (represented by shape) over the 18-day experimental period (Day represented by color). Temperature treatments: (A) 29°C, (B) 25°C, (C) 20°C, (D) 18°C, and (E) 15°C. Data from this figure represent the same data from Figure 2C-D.
# Figure 5: Mean dissimilarity in SynComs at Day 18 compared to Day 0. Dissimilarity was calculated as the mean Bray-Curtis distance between individual Day 18 and Day 0 SynComs within each temperature treatment. SynComs containing Bd (A) or Bsal (B).

#### Packages ####
library(ggplot2)
library(vegan)
library(dplyr)
library(phyloseq)
library(ape)
library(tidyr)
library(microbiome)
library(vegan)
library(patchwork)
library(cowplot)
library(microViz)
library(ggpubr)
library(tibble)

#### Reading in Data ####
sv<- read.csv("featuretable16S 2.csv")
row.names(sv) <- sv$OTU_ID
sv$OTU_ID <- NULL

#Taxonomy#
taxa <- read.csv("taxa_columns.csv")
row.names(taxa) <- taxa$Feature_ID
taxa$Feature_ID <- NULL

#Metadata#
metadata <- read.csv("Metadata16sThesis.csv")

# Convert dataframes into phyloseq format

sv = otu_table(sv, taxa_are_rows = TRUE)

taxa = tax_table(as.matrix(taxa))

metadata <- sample_data(metadata)
rownames(metadata) <- metadata$SampleID
metadata$SampleID <- NULL

metadata.2 <- metadata[!is.na(metadata$Temperature), ]

# merge SV & taxa data into phyloseq object
physeq <- phyloseq(sv, taxa)

# create random phylo tree with ASV names as tips
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))

# create new phyloseq object with all info
physeq16s <- phyloseq(sv, taxa, metadata.2, random_tree)

physeq_filtered <- subset_samples(physeq16s, Contents != "Synth")

physeq_filtered1 <- subset_samples(physeq_filtered, Contents != "Bd")

physeq_filtered2 <- subset_samples(physeq_filtered1, Contents != "Bsal")

physeq_filtered3 <- subset_samples(physeq_filtered2, Contents != "Ext_cntrl")

physeq_filtered4 <- subset_samples(physeq_filtered3, Contents != "PCR_cntrl")

physeq_filtered5 <- subset_samples(physeq_filtered4, Contents != "Media")

####Pruning Low Abundance Taxa ####

taxa_to_keep <- taxa_sums(physeq_filtered5) > 800
ps_pruned <- prune_taxa(taxa_to_keep, physeq_filtered5)

sample_data(physeq_filtered5)
ps_sorted <- ps_arrange(physeq_filtered5, Temperature)
sample_data(ps_sorted)

####Mean Reads Per Sample ####
total_reads_per_sample <- sort(sample_sums(physeq_filtered5))
mean_reads <- mean(total_reads_per_sample)
print(mean_reads)

#### Seperating Bd and Bsal####

physeq_filtered_Bd_Only <- subset_samples(ps_pruned, Contents != "Bsal_Synth")

sample_info <- sample_data(physeq_filtered_Bd_Only)

physeq_filtered_Bsal_Only <- subset_samples(ps_pruned, Contents != "Bd_Synth")

ps_sorted_Bsal <- ps_arrange(physeq_filtered_Bsal_Only, Temperature)
sample_data(ps_sorted_Bsal)

####Rarification####

#Bd#
set.seed(123)
ps.rarefied.Bd=rarefy_even_depth(physeq_filtered_Bd_Only, rngseed=1, sample.size =8000, replace=F)
ps.rarefied.Bd=rarefy_even_depth(physeq_filtered_Bd_Only, rngseed=1, sample.size =5000, replace=F)

sample_info_rare <- sample_data(ps.rarefied.Bd)

physeq16s.log.Bd <- transform_sample_counts(ps.rarefied.Bd, function(x) log(1 + x))
physeq16s.ord.Bd <- ordinate(physeq16s.log.Bd, method = "NMDS", distance = "bray")

sample_data(physeq16s.log.Bd)$Temperature <- as.factor(sample_data(physeq16s.log.Bd)$Temperature)
sample_data(physeq16s.log.Bd)$Day <- as.factor(sample_data(physeq16s.log.Bd)$Day)
#sample_data(physeq16s.log.Bd)

#Bsal#
set.seed(1)
ps.rarefied.Bsal=rarefy_even_depth(physeq_filtered_Bsal_Only, rngseed=1, sample.size =8000, replace=F) #try 5000, or run on un rare data 

physeq16s.log.Bsal <- transform_sample_counts(ps.rarefied.Bsal, function(x) log(1 + x))
physeq16s.ord.Bsal <- ordinate(physeq16s.log.Bsal, method = "NMDS", distance = "bray")

sample_data(physeq16s.log.Bsal)$Temperature <- as.factor(sample_data(physeq16s.log.Bsal)$Temperature)
sample_data(physeq16s.log.Bsal)$Day <- as.factor(sample_data(physeq16s.log.Bsal)$Day)

ps_sorted_Bsal_rare <- ps_arrange(ps.rarefied.Bsal, Temperature)
sample_data(ps_sorted_Bsal_rare)

####Supplemental Figure 2 ####

Bd_Day_Only <- plot_ordination(physeq16s.log.Bd, physeq16s.ord.Bd, shape = "Day") +
  geom_point(size = 2) +
  theme_bw() 
print(Bd_Day_Only)

Bsal_Day_Only <- plot_ordination(physeq16s.log.Bsal, physeq16s.ord.Bsal, shape = "Day") +
  geom_point(size = 2) +
  theme_bw() 
print(Bsal_Day_Only)

Day_Only_Bd_Bsal<- Bd_Day_Only / Bsal_Day_Only
Day_Only_Bd_Bsal_w_Tags<-Day_Only_Bd_Bsal + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = 'bold'))

####Supplemental Figure 3 ####

set.seed(123)
ps.rarefied.both=rarefy_even_depth(physeq_filtered5, rngseed=1, sample.size =5000, replace=F)

physeq16s.log.both <- transform_sample_counts(ps.rarefied.both, function(x) log(1 + x))
physeq16s.ord.both <- ordinate(physeq16s.log.both, method = "NMDS", distance = "bray")

sample_data(physeq16s.log.both)$Temperature <- as.factor(sample_data(physeq16s.log.both)$Temperature)
sample_data(physeq16s.log.both)$Day <- as.factor(sample_data(physeq16s.log.both)$Day)


set.seed(123)
ps_T29_both <- subset_samples(physeq16s.log.both, Temperature == 29)
ord_T29_both <- ordinate(ps_T29_both, method = "NMDS", distance = "bray")


Bsal_T29_both <- plot_ordination(ps_T29_both, ord_T29_both, shape = "Contents", color = "Day") +
  geom_point(size = 5) +
  labs(title = "29") +
  scale_color_manual(values = c("0" = "#Bdd7ee", "18" = "#272838", "12" = "#1f78b4", "6" = "#68AEd6"))+
  theme_classic()  + 
  theme(legend.position = "right")+
  labs(title = "29", shape = "Chytrid")
print(Bsal_T29_both)

set.seed(123)
ps_T25_both <- subset_samples(physeq16s.log.both, Temperature == 25)
ord_T25_both <- ordinate(ps_T25_both, method = "NMDS", distance = "bray")

Bsal_T25_both <- plot_ordination(ps_T25_both, ord_T25_both, shape = "Contents", color = "Day") +
  geom_point(size = 5) +
  labs(title = "25") +
  scale_color_manual(
    values = c("0" = "#Bdd7ee", "6" = "#68AEd6", "12" = "#1f78b4", "18" = "#272838")
  ) +
  theme_classic() + 
  theme(legend.position = "right")+
  labs(title = "25", shape = "Chytrid")

print(Bsal_T25_both)

set.seed(123)
ps_T20_both <- subset_samples(physeq16s.log.both, Temperature == 20)
ord_T20_both <- ordinate(ps_T20_both, method = "NMDS", distance = "bray")


Bsal_T20_both <- plot_ordination(ps_T20_both, ord_T20_both, shape = "Contents", color = "Day") +
  geom_point(size = 5) +
  labs(title = "20") +
  scale_color_manual(values = c("0" = "#Bdd7ee", "18" = "#272838", "12" = "#1f78b4", "6" = "#68AEd6"))+
  theme_classic()  + 
  theme(legend.position = "right")+ 
  labs(title = "20", shape = "Chytrid")
print(Bsal_T20_both)

set.seed(123)
ps_T18_both <- subset_samples(physeq16s.log.both, Temperature == 18)
ord_T18_both <- ordinate(ps_T18_both, method = "NMDS", distance = "bray")


Bsal_T18_both <- plot_ordination(ps_T18_both, ord_T18_both, shape = "Contents", color = "Day") +
  geom_point(size = 5) +
  labs(title = "18") +
  scale_color_manual(values = c("0" = "#Bdd7ee", "18" = "#272838", "12" = "#1f78b4", "6" = "#68AEd6"))+
  theme_classic()  + 
  theme(legend.position = "right")+ 
  labs(title = "18", shape = "Chytrid")
print(Bsal_T18_both)

set.seed(123)
ps_T15_both <- subset_samples(physeq16s.log.both, Temperature == 15)
ord_T15_both <- ordinate(ps_T15_both, method = "NMDS", distance = "bray")


Bsal_T15_both <- plot_ordination(ps_T15_both, ord_T15_both, shape = "Contents", color = "Day") +
  geom_point(size = 5) +
  labs(title = "15") +
  scale_color_manual(values = c("0" = "#Bdd7ee", "18" = "#272838", "12" = "#1f78b4", "6" = "#68AEd6"))+
  theme_classic()  + 
  theme(legend.position = "right")+
  labs(title = "15", shape = "Chytrid")
print(Bsal_T15_both)

Bsal_DayBy1Temp_both <- (Bsal_T29_both | Bsal_T25_both | Bsal_T20_both) / (Bsal_T18_both | Bsal_T15_both) + plot_layout(guides = "collect") &
  theme(legend.position = "right")
Bsal_DayBy1Temp_wTags_both <- Bsal_DayBy1Temp_both + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = 'bold'))
print(Bsal_DayBy1Temp_wTags_both)

####Supplemental Figure 5 ####

##Bd 

bray.dist.bd <- phyloseq::distance(physeq_filtered_Bd_Only, method="bray")

bray_df_Bd <- broom::tidy(bray.dist.bd) %>%
  rename(sample_a = item1, sample_b = item2, bray.dist = distance)

env.df <- as.data.frame(as.matrix(sample_data(ps_pruned)))

bray_combined_bd <- bray_df_Bd %>%
  left_join(env.df %>% rownames_to_column("sample_a"), by = "sample_a") %>%
  left_join(env.df %>% rownames_to_column("sample_b"), by = "sample_b", suffix = c("_a", "_b"))

bray_baseline_bd <- bray_combined_bd %>%
  filter(Day_a == Day_b) %>%
  mutate(Baseline_Is_Day_0 = if_else(Day_a == "0", Day_b, Day_a)) %>% 
  group_by( Temperature = Temperature_b, Day = Baseline_Is_Day_0) %>% # group data by the time x treatment
  summarize(mean_dist_from_start = mean(bray.dist),
            se_dist = sd(bray.dist) / sqrt(n()),
            groups = "drop")

bray_baseline_bd_18<- bray_baseline_bd[bray_baseline_bd$Day == 18, ]

Bd_Baseline_Plot<- ggplot(bray_baseline_bd_18, aes(x = Day, y = mean_dist_from_start, fill = Temperature)) +
  geom_col( stat = "identity", position = position_dodge(width = 0.9), color = "black") +
  geom_errorbar(aes(ymin = mean_dist_from_start - se_dist, 
                    ymax = mean_dist_from_start + se_dist), width = 0.1, position = position_dodge(0.9)) +
  labs( x = "Day 18",
        y = "Mean Bray-Curtis Distance (Relative to T1)") +
  theme_minimal()+
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = PalleteAdivGTDB, name = "Temperature") +
  theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank())
print(Bd_Baseline_Plot)

##Bsal 
bray.dist.bsal <- phyloseq::distance(physeq_filtered_Bsal_Only, method="bray")

bray_df_bsal <- broom::tidy(bray.dist.bsal) %>%
  rename(sample_a = item1, sample_b = item2, bray.dist = distance)

env.df <- as.data.frame(as.matrix(sample_data(ps_pruned)))

bray_combined_bsal <- bray_df_bsal %>%
  left_join(env.df %>% rownames_to_column("sample_a"), by = "sample_a") %>%
  left_join(env.df %>% rownames_to_column("sample_b"), by = "sample_b", suffix = c("_a", "_b"))
bray_baseline_bsal <- bray_combined_bsal%>%
  filter(Day_a == Day_b) %>%
  mutate(Day_0 = if_else(Day_a == "0", Day_b, Day_a)) %>% 
  group_by( Temperature = Temperature_b, Day = Day_0) %>% # group data by the time x treatment
  summarize(mean_dist_from_start = mean(bray.dist, na.rm = TRUE),
            se_dist = sd(bray.dist, na.rm = TRUE) / sqrt(n()),
            groups = "drop")

bray_baseline_bsal_18<- bray_baseline_bsal[bray_baseline_bsal$Day == 18, ]


Bsal_Baseline_Plot<- ggplot(bray_baseline_bsal_18, aes(x = Day, y = mean_dist_from_start, fill = Temperature)) +
  geom_col( stat = "identity", position = position_dodge(width = 0.9), color = "black") +
  geom_errorbar(aes(ymin = mean_dist_from_start - se_dist, 
                    ymax = mean_dist_from_start + se_dist), width = 0.1, position = position_dodge(0.9)) +
  labs(y = "Mean Bray-Curtis Distance (Relative to T1)",
    x = "Day 18") +
  theme_minimal()+
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = PalleteAdivGTDB, name = "Temperature") +
  theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank())
print(Bsal_Baseline_Plot)

## Combine Bd and Bsal 
Baseline_Eco_Graphs<- Bd_Baseline_Plot + Bsal_Baseline_Plot + plot_layout(guides = "collect") &
  theme(legend.position = "right")
Baseline_Eco_Graphs_w_Tags<-Baseline_Eco_Graphs + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = 'bold'))
print(Baseline_Eco_Graphs_w_Tags)

