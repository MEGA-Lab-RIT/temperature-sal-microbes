#### Figure_3####
###Created on: 12/14/25 
###Emma Thompson 
# Creating Figure 2: Relative abundance of bacteria in Bd (A) and Bsal (B) microcosms by both temperature and time. Samples are ordered by time within each temperature subset, with each bar representing individual replicates. The family Mycobacteriaceae disappears over time; the time point at which Mycobacteriaceae abundance reached zero is denoted with an asterisk. Change in community composition between experimental temperatures vs. pathogen thermal optima (as measured by Bray-Curtis dissimilarity) in microcosms containing Bd (C) or Bsal (D). Based on the literature, the thermal optimum for Bd was 20°C and for Bsal was 15°C.

#### Packages ####
library(ggplot2)
library(vegan)
library(dplyr)
library(phyloseq)
library(ape)
library(tidyr)
library(microbiome)
library(stringr)
library(patchwork)
library(cowplot)
library(microViz)
library(lme4)
library(lmerTest)
library(tidyverse)
library(emmeans)


#### Reading In Data####
sv <- read.csv("featuretable16S 2.csv")
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

#Removing controls 
physeq_filtered <- subset_samples(physeq16s, Contents != "Synth")

physeq_filtered1 <- subset_samples(physeq_filtered, Contents != "Bd")

physeq_filtered2 <- subset_samples(physeq_filtered1, Contents != "Bsal")

physeq_filtered3 <- subset_samples(physeq_filtered2, Contents != "Ext_cntrl")

physeq_filtered4 <- subset_samples(physeq_filtered3, Contents != "PCR_cntrl")

physeq_filtered5 <- subset_samples(physeq_filtered4, Contents != "Media")

#### Pruning low abundance reads ####

taxa_to_keep <- taxa_sums(physeq_filtered5) > 800
ps_pruned <- prune_taxa(taxa_to_keep, physeq_filtered5)

#### Seperating Bd and Bsal####

physeq_filtered_Bd_Only <- subset_samples(ps_pruned, Contents != "Bsal_Synth")

physeq_filtered_Bsal_Only <- subset_samples(ps_pruned, Contents != "Bd_Synth")

####Rarification####

#Bd#
set.seed(123)
ps.rarefied.Bd=rarefy_even_depth(physeq_filtered_Bd_Only, rngseed=1, sample.size =8000, replace=F)

physeq16s.log.Bd <- transform_sample_counts(ps.rarefied.Bd, function(x) log(1 + x))
physeq16s.ord.Bd <- ordinate(physeq16s.log.Bd, method = "NMDS", distance = "bray")

sample_data(physeq16s.log.Bd)$Temperature <- as.factor(sample_data(physeq16s.log.Bd)$Temperature)
sample_data(physeq16s.log.Bd)$Day <- as.factor(sample_data(physeq16s.log.Bd)$Day)

physeq16s.rel.GDTB.Bd = transform_sample_counts(physeq16s.log.Bd, function(x) x/sum(x)*100)


#Bsal#
set.seed(1)
ps.rarefied.Bsal=rarefy_even_depth(physeq_filtered_Bsal_Only, rngseed=1, sample.size =8000, replace=F)

physeq16s.log.Bsal <- transform_sample_counts(ps.rarefied.Bsal, function(x) log(1 + x))
physeq16s.ord.Bsal <- ordinate(physeq16s.log.Bsal, method = "NMDS", distance = "bray")

sample_data(physeq16s.log.Bsal)$Temperature <- as.factor(sample_data(physeq16s.log.Bsal)$Temperature)
sample_data(physeq16s.log.Bsal)$Day <- as.factor(sample_data(physeq16s.log.Bsal)$Day)

physeq16s.rel.GDTB.Bsal = transform_sample_counts(physeq16s.log.Bsal, function(x) x/sum(x)*100)

#### Palette 4 ####

palette4 <- c( "#CAB2D6","#6A3D9A", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "grey", "#FFFF99","#FF7F00", "#FDBF6F", "#FB9A99")

####Relative Abundance Graph Bd####
glom <- tax_glom(physeq16s.rel.GDTB.Bd, taxrank = 'Family', NArm = FALSE)
ps.melt10 <- psmelt(glom)

ps.melt10$Family <- as.character(ps.melt10$Family)

ps.melt10 <- ps.melt10 %>%
  group_by(Family) %>%
  mutate(median=median(Abundance))

keep <- unique(ps.melt10$Family[ps.melt10$Abundance > 3.75]) 

ps.melt10$Family[!(ps.melt10$Family %in% keep)] <- "Other"

ps.melt20 <- ps.melt10 %>%
  group_by(Temperature, Day, Sample, Family) %>%
  dplyr::summarise(Abundance=sum(Abundance))

ps.melt20 <- ps.melt20 %>%
  arrange(Temperature, Day, Sample) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

ps.melt20 <- ps.melt20 %>%
  mutate(Sample_ordered = interaction(Temperature, Day, Sample, sep = "_")) %>%
  arrange(Temperature, Day) %>%
  mutate(Sample_ordered = factor(Sample_ordered, levels = unique(Sample_ordered)))


GTDB_Bd <- ggplot(ps.melt20, aes(x = Sample_ordered, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity", aes(fill=factor (Family))) + 
  labs(title= "A", x="", y="Relative Abundance (%)") +
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right", 
        legend.key.size = unit(1,'cm'),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())+
  scale_fill_manual(values= c ("Xanthomonadaceae" = "#FB9A99",
                               "Weeksellaceae" = "#FDBF6F",
                               "Sphingobacteriaceae" = "#FF7F00", 
                               "Pseudomonadaceae"= "#FFFF99",
                               "Other"= "gray",
                               "Mycobacteriaceae" = "#2E6F40", 
                               "Moraxellaceae" = "#B2DF8A",
                               "Microbacteriaceae" = "#1F78B4",
                               "Enterobacteriaceae" = "#6A3D9A", 
                               "Burkholderiaceae" = "#CAB2D6"), 
                    name="Family")+
  facet_grid(.~factor(Temperature, levels = c('15', '18',
                                              '20',
                                              '25', '29')), 
             scales = "free", switch = "x", space = "free_x")
print(GTDB_Bd)

relabund_sorted <- ps.melt20 %>% arrange(desc(Abundance))

relabund_sorted <- relabund_sorted %>% filter(Family != "Microbacteriaceae")

mycobacteriaecea_relab <- relabund_sorted %>% 
  filter(Family == "Mycobacteriaceae")

Moraxellaceae_relab <- relabund_sorted %>% 
  filter(Family == "Moraxellaceae")

####Relative Abundance Graph Bsal ####

glom_Bsal <- tax_glom(physeq16s.rel.GDTB.Bsal, taxrank = 'Family', NArm = FALSE)
ps.melt10_Bsal <- psmelt(glom_Bsal)

ps.melt10_Bsal$Family <- as.character(ps.melt10_Bsal$Family)

ps.melt10_Bsal <- ps.melt10_Bsal %>%
  group_by(Family) %>%
  mutate(median=median(Abundance))

keep <- unique(ps.melt10_Bsal$Family[ps.melt10_Bsal$Abundance > 3.75]) 

ps.melt10_Bsal$Family[!(ps.melt10_Bsal$Family %in% keep)] <- "Other"

ps.melt20_Bsal <- ps.melt10_Bsal %>%
  group_by(Temperature, Day, Sample, Family) %>%
  dplyr::summarise(Abundance=sum(Abundance))

ps.melt20_Bsal <- ps.melt20_Bsal %>%
  arrange(Temperature, Day, Sample) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

ps.melt20_Bsal <- ps.melt20_Bsal %>%
  mutate(Sample_ordered = interaction(Temperature, Day, Sample, sep = "_")) %>%
  arrange(Temperature, Day) %>%
  mutate(Sample_ordered = factor(Sample_ordered, levels = unique(Sample_ordered)))


GTDB_Bsal <- ggplot(ps.melt20_Bsal, aes(x = Sample_ordered, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity", aes(fill=factor (Family))) + 
  labs(title= "B", x="", y="Relative Abundance (%)") +
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right", 
        legend.key.size = unit(1,'cm'),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())+
  scale_fill_manual(values= c ("Xanthomonadaceae" = "#FB9A99",
                               "Weeksellaceae" = "#FDBF6F",
                               "Sphingobacteriaceae" = "#FF7F00", 
                               "Pseudomonadaceae"= "#FFFF99",
                               "Other"= "gray",
                               "Mycobacteriaceae" = "#2E6F40", 
                               "Moraxellaceae" = "#B2DF8A",
                               "Microbacteriaceae" = "#1F78B4",
                               "Enterobacteriaceae" = "#6A3D9A", 
                               "Burkholderiaceae" = "#CAB2D6"), 
                    name="Family")+
  facet_grid(.~factor(Temperature, levels = c('15', '18',
                                              '20',
                                              '25', '29')), 
             scales = "free", switch = "x", space = "free_x")
print(GTDB_Bsal)


Relative_Abundance_GTDB<- GTDB_Bd + GTDB_Bsal + plot_layout(guides = "collect") &
  theme(legend.position = "right")



#### Community Composition Relative to Pathogen Optima ####

#### Creating the distance matrix from physeq object ####

####Bd

bray.dist.bd <- phyloseq::distance(physeq_filtered_Bd_Only, method="bray")

bray_df_Bd <- broom::tidy(bray.dist.bd) %>%
  rename(sample_a = item1, sample_b = item2, bray.dist = distance)

env.df <- as.data.frame(as.matrix(sample_data(ps_pruned)))


bray_combined_bd <- bray_df_Bd %>%
  left_join(env.df %>% rownames_to_column("sample_a"), by = "sample_a") %>%
  left_join(env.df %>% rownames_to_column("sample_b"), by = "sample_b", suffix = c("_a", "_b"))

baseline_Bd <- 20

bray_averages_bd <- bray_combined_bd %>%
  filter(Day_a == Day_b, Temperature_a == baseline_Bd | Temperature_b == baseline_Bd) %>%
  mutate(comparison_temp = if_else(Temperature_a == "20", Temperature_b, Temperature_a)) %>% 
  group_by(Day = Day_a, Temperature = comparison_temp) %>% # group data by the time x treatment
  summarize(                            
    avg_dist = mean(bray.dist, na.rm = TRUE),
    sd_dist  = sd(bray.dist),
    n        = n(),
    se_dist  = sd_dist / sqrt(n),
    .groups = "drop"
  )

bray_averages_bd_2 <- bray_averages_bd %>% filter(Temperature != "20")

Temp_Day_Eco_Bd<-ggplot(bray_averages_bd_2, aes(x = Day, y = avg_dist, color = Temperature, group = Temperature)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = avg_dist - se_dist, ymax = avg_dist + se_dist), width = 0.1) +
  labs(
    title = "C",
    subtitle = "Pathogen Optimum 20°C",
    x = "Time (days)",
    y = "Mean Bray-Curtis Dissimilarity (±SE)"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  theme_bw()+
  theme(plot.title = element_text(face = "bold"), panel.grid = element_blank(), text = element_text(size = 20), axis.text.y = element_text(size = 15))+
  scale_color_manual(values = c("15" = "#f2aa84", "18" = "#83e291", "25" = "#d86ecc", "29" = "#ffa1e7"))
print(Temp_Day_Eco_Bd)


####Bsal Eco Graph ####

bray.dist.bsal <- phyloseq::distance(physeq_filtered_Bsal_Only, method="bray")

bray_df_bsal <- broom::tidy(bray.dist.bsal) %>%
  rename(sample_a = item1, sample_b = item2, bray.dist = distance)

env.df <- as.data.frame(as.matrix(sample_data(ps_pruned)))

bray_combined_bsal <- bray_df_bsal %>%
  left_join(env.df %>% rownames_to_column("sample_a"), by = "sample_a") %>%
  left_join(env.df %>% rownames_to_column("sample_b"), by = "sample_b", suffix = c("_a", "_b"))

baseline <- 15

bray_averages_bsal <- bray_combined_bsal %>%
  filter(Day_a == Day_b, Temperature_a == baseline | Temperature_b == baseline) %>%
  mutate(comparison_temp = if_else(Temperature_a == "15", Temperature_b, Temperature_a)) %>% 
  group_by(Day = Day_a, Temperature = comparison_temp) %>% # group data by the time x treatment
  summarize(                            
    avg_dist = mean(bray.dist, na.rm = TRUE),
    sd_dist  = sd(bray.dist, na.rm= TRUE),
    n        = n(),
    se_dist  = sd_dist / sqrt(n),
    .groups = "drop"
  )

bray_averages_bsal_2 <- bray_averages_bsal %>% filter(Temperature != "15")

Temp_Day_Eco_Bsal<-ggplot(bray_averages_bsal_2, aes(x = Day, y = avg_dist, color = Temperature, group = Temperature)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = avg_dist - se_dist, ymax = avg_dist + se_dist), width = 0.1) +
  labs(
    title = "D",
    subtitle = "Pathogen Optimum: 15°C",
    x = "Time (days)",
    y = "Mean Bray-Curtis Dissimilarity (±SE)"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  theme_bw()+
  theme(plot.title = element_text(face = "bold"), panel.grid = element_blank(), text = element_text(size = 20), axis.text.y = element_text(size = 15))+
  scale_color_manual(values = c("20" = "#96dcf8", "18" = "#83e291", "25" = "#d86ecc", "29" = "#ffa1e7"))
print(Temp_Day_Eco_Bsal)


####Combine to Make Figure_3####
Figure3<- Relative_Abundance_GTDB / ((Temp_Day_Eco_Bd / Temp_Day_Eco_Bsal) )

ggsave("Figure3.png", plot= Figure3, height = 18, width = 14, units = "in")
