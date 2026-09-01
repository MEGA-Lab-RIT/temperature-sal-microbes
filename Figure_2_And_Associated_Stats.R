#### Figure_2####
###Created on: 12/14/25 
###Emma Thompson 
# Creating Figure 2: Effect of temperature on pathogen inhibition and community composition. Percent inhibition on Day 18 for Bd (A) and Bsal (B) in each temperature trial as measured by competition assay. Median and interquartile range (IQR) with whiskers extending to 1.5X IQR, with outliers denoted as solid black points and all points as gray points. Inhibition was analyzed via one -way ANOVA Bd: F(4,20) = 2.2, p = 0.11 & Bsal: F(4,20) = 4.5, p = 0.01. Community Composition was examined using a PERMANOVA analysis on weighted UniFrac distances for communities with Bd (C) or Bsal (D). Time is represented by shape and temperature represented by color. Time (Bd: F(3,20) =58.8, p= 0.001 & Bsal: F(3,57)= 79.35, p= 0.001), temperature (Bd: F(4,20) =2.2, p= 0.106 & Bsal: F(4,57)= 4.3, p= 0.02), and their interaction (Bd: F(11,20) =2.35, p= 0.006 & Bsal: F(12,57)= 2.39, p= 0.001) significantly influence community composition. 

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





#### A & B: Pathogen Inhibition ####

## Bd Inhibition Data ##
Bd_data <- read.csv("Thesis_BD_Metabolites_R.csv")

Bd_data1 <- Bd_data %>%
  filter(Day == 18)

Bd_data$Temperature <- factor(Bd_data$Temperature)

PalleteAdivGTDB <- c("#f2aa84","#83e291", "#96dcf8", "#d86ecc", "#ffa1e7")

Bd_Inhibition <-
  ggplot( data = Bd_data, aes(
    x = Temperature,
    y = Mean,
    fill = Temperature
  )) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 1.5, alpha = 0.5, position = position_jitter(0.2)) +
  labs(title= "A", x = "", y = "Pathogen Inhibition (%)") +
  theme_bw() +
  coord_cartesian(ylim = c(0, 105)) +
  theme(
    legend.position = "right",
    legend.key.size = unit(1, 'cm'),
    axis.text.x = element_text()  # or element_blank() if you prefer hidden labels
  ) +
  scale_fill_manual(values = PalleteAdivGTDB, name = "Temperature") 

Bsal_data <- read.csv("Thesis_Bsal_Metabolites_R.csv")

Bsal_data1 <- Bsal_data %>%
  filter(Day == 18)

Bsal_data$Temperature <- factor(Bsal_data$Temperature)

Bsal_Inhibition <-
  ggplot( data = Bsal_data, aes(
    x = Temperature,
    y = Mean,
    fill = Temperature
  )) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 1.5, alpha = 0.5, position = position_jitter(0.2)) +
  labs(title= "", x = "", y = "Pathogen Inhibition (%)") +
  theme_bw() +
  coord_cartesian(ylim = c(0, 105)) +
  theme(
    legend.position = "right",
    legend.key.size = unit(1, 'cm'),
    axis.text.x = element_text()  # or element_blank() if you prefer hidden labels
  ) +
  scale_fill_manual(values = PalleteAdivGTDB, name = "Temperature") +
  scale_y_continuous(limits = c(0, 100))

#### D & C: Comunity Composition ####

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

## Bd plot ##
Day_Temp_Bd_Only <- plot_ordination(physeq16s.log.Bd, physeq16s.ord.Bd, shape = "Day", color = "Temperature") +
  geom_point(size = 2) +
  scale_color_manual(values = c("15" = "#f2aa84", "18" = "#83e291", "20" = "#96dcf8", "25" = "#d86ecc", "29" = "#ffa1e7"))+
  theme_bw() 
print(Day_Temp_Bd_Only)

##Bsal Plot ##
Day_Temp_Bsal_Only <- plot_ordination(physeq16s.log.Bsal, physeq16s.ord.Bsal, shape = "Day", color = "Temperature") +
  geom_point(size = 2) +
  scale_color_manual(values = c("15" = "#f2aa84", "18" = "#83e291", "20" = "#96dcf8", "25" = "#d86ecc", "29" = "#ffa1e7"))+
  theme_bw() 
print(Day_Temp_Bsal_Only)

#### Combining A,B,C, & D ####

Figure_2 <- (Bd_Inhibition + Bsal_Inhibition) / (Day_Temp_Bd_Only + Day_Temp_Bsal_Only + plot_layout(guides = "collect") &
                                                       theme(legend.position = "right"))

ggsave("Figure_2.png", plot= Figure_2, height = 10, width = 10, units = "in")


#### Supplemental Table 4: PERMANOVA of Day and Temperature's Effect on Community Comp####

Distance_Matrix_Bsal <- distance(physeq16s.log.Bsal, method = "bray")
metadataPERM_Bsal <- as(sample_data(physeq16s.log.Bsal), "data.frame")

permanova_Bsal <- adonis2(Distance_Matrix_Bsal ~ Day + Temperature, data = metadataPERM_Bsal, permutations = 999, by = "terms")
print(permanova_Bsal)

Distance_Matrix_Bd <- distance(physeq16s.log.Bd, method = "bray")
metadataPERM_Bd <- as(sample_data(physeq16s.log.Bd), "data.frame")

permanova_Bd <- adonis2(Distance_Matrix_Bd ~ Day + Temperature, data = metadataPERM_Bd, permutations = 999, by = "terms")
print(permanova_Bd)

####Supplemental Table 5 Beta dispersion of Day and Temperature####

Day_Factor_Bsal <- factor(metadataPERM_Bsal$Day)

BetaDisDay_Bsal <- betadisper(Distance_Matrix_Bsal, Day_Factor_Bsal)
anova(BetaDisDay_Bsal)

Temp_Factor_Bsal <- factor(metadataPERM_Bsal$Temperature)

BetaDisTemp_Bsal <- betadisper(Distance_Matrix_Bsal, Temp_Factor_Bsal)
anova(BetaDisTemp_Bsal)

Day_Factor_Bd <- factor(metadataPERM_Bd$Day)

BetaDisDay_Bd <- betadisper(Distance_Matrix_Bd, Day_Factor_Bd)
anova(BetaDisDay_Bd)

Temp_Factor_Bd <- factor(metadataPERM_Bd$Temperature)

BetaDisTemp_Bd <- betadisper(Distance_Matrix_Bd, Temp_Factor_Bd)
anova(BetaDisTemp_Bd)

#### Supplemental Table 6: Pathogen Inhibition ANOVA ####

Bd_18_Only <- Bd_data %>%
  filter(Day == 18)

Bd_18_Only <- Bd_18_Only %>%
  mutate(Temperature = factor(Temperature))

Bd_anova <- aov(Mean ~ Temperature, data = Bd_18_Only)
summary(Bd_anova)

Bsal_18_Only <- Bsal_data %>%
  filter(Day == 18)

Bsal_anova <- aov(Mean ~ Temperature, data = Bsal_18_Only)
summary(Bsal_anova)

#### Supplemental Table 7: PostHoc Test of Temperatures for Bsal ####

Bsal_PostHoc<-TukeyHSD(Bsal_anova, "Temperature")
print(Bsal_PostHoc)







