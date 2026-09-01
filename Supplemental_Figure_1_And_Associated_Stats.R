#### Supplemental_Figure_1####
###Created on: 11/11/25  
###Emma Thompson 
# Creating Supplemental Figure 1: Comparisons of alpha diversity between temperature treatments over time. Alpha diversity in SynComs containing Bd (A) or Bsal (B) is shown for all metrics: Shannon diversity, Chao1, and evenness indices. Day is on the x-axis and boxplot color corresponds to the trial temperature. Boxplots depict the median and interquartile range (IQR) with whiskers extending to 1.5X IQR, with outliers denoted as solid black points and all points as gray points. Letters represent significant interactions between temperature and day as determined by Tukey post hoc test.

#### Packages ####
library(ggplot2)
library(vegan)
library(dplyr)
library(phyloseq)
library(ape)
library(tidyr)
library(microbiome)
library(ANCOMBC)
library(ggsignif)
library(ggpubr)
library(patchwork)
library(stringr)
library("multcompView")
library(purrr)

#### Reading in Data####
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


####Removing Controls####

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

#### Alpha Diveristy Bd ####

AdivBd_Only <- suppressWarnings( data.frame(
  "Shannon" = phyloseq::estimate_richness(physeq_filtered_Bd_Only, measures = "Shannon")$Shannon,
  "Chao1" = phyloseq::estimate_richness(physeq_filtered_Bd_Only, measures = "Chao1")$Chao1,
  "Observed" = phyloseq::estimate_richness(physeq_filtered_Bd_Only, measures = "Observed")$Observed,
  "Day" = phyloseq::sample_data(physeq_filtered_Bd_Only)$Day,
  "Temperature" = phyloseq::sample_data(physeq_filtered_Bd_Only)$Temperature))

AdivBd_Only$Evenness <- AdivBd_Only$Shannon / log(AdivBd_Only$Observed)

AdivBd_Only <- AdivBd_Only[AdivBd_Only$Shannon != 0, ]


PalleteAdivGTDB <- c("#f2aa84","#83e291", "#96dcf8", "#d86ecc", "#ffa1e7")

Bd_Only_Adiv_plot <- AdivBd_Only %>%
  gather(key = metric, value = value, c("Shannon", "Chao1", "Evenness")) %>%
  mutate(
    metric = factor(metric, levels = c("Shannon", "Chao1", "Evenness")),
    Temperature = factor(Temperature, levels = c("15","18","20","25","29")), 
    Day = factor(Day, levels = c("0","6","12","18"))
  ) %>%
  ggplot(aes(
    x = Day,
    y = value,
    fill = Temperature
  )) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 1.5, alpha = 0.5, position = position_jitter(0.2)) +
  labs(title= "A", x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(1, 'cm'),
    axis.text.x = element_text()  # or element_blank() if you prefer hidden labels
  ) +
  scale_fill_manual(values = PalleteAdivGTDB, name = "Temperature") 

#### Bsal Alpha Diversity ####

AdivBsal_Only <- suppressWarnings( data.frame(
  "Shannon" = phyloseq::estimate_richness(physeq_filtered_Bsal_Only, measures = "Shannon")$Shannon,
  "Chao1" = phyloseq::estimate_richness(physeq_filtered_Bsal_Only, measures = "Chao1")$Chao1,
  "Observed" = phyloseq::estimate_richness(physeq_filtered_Bsal_Only, measures = "Observed")$Observed,
  "Day" = phyloseq::sample_data(physeq_filtered_Bsal_Only)$Day,
  "Temperature" = phyloseq::sample_data(physeq_filtered_Bsal_Only)$Temperature))

AdivBsal_Only$Evenness <- AdivBsal_Only$Shannon / log(AdivBsal_Only$Observed)

AdivBsal_Only <- AdivBsal_Only[AdivBsal_Only$Shannon != 0, ]


PalleteAdivGTDB <- c("#f2aa84","#83e291", "#96dcf8", "#d86ecc", "#ffa1e7")

Bsal_Only_Adiv_plot <- AdivBsal_Only %>%
  gather(key = metric, value = value, c("Shannon", "Chao1", "Evenness")) %>%
  mutate(
    metric = factor(metric, levels = c("Shannon", "Chao1", "Evenness")),
    Temperature = factor(Temperature, levels = c("15","18","20","25","29")), 
    Day = factor(Day, levels = c("0","6","12","18"))
  ) %>%
  ggplot(aes(
    x = Temperature,
    y = value,
    fill = Temperature
  )) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 1.5, alpha = 0.5, position = position_jitter(0.2)) +
  labs(title = "B", x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(1, 'cm'),
    axis.text.x = element_text()  # or element_blank() if you prefer hidden labels
  ) +
  scale_fill_manual(values = PalleteAdivGTDB, name = "Temperature") 


#### Adding the letters to the Adiv plots ####

##Bd
metrics_list <- c("Shannon", "Chao1", "Evenness")

adiv_long_Bd <- AdivBd_Only %>%
  gather(key = metric, value = value, all_of(metrics_list)) %>%
  mutate(metric = factor(metric, levels = metrics_list),
         Temperature = factor(Temperature, levels = c("15","18","20","25","29")),
         Day = factor(Day, levels = c("0","6","12","18")))

letters.df_Bd <- map_df(metrics_list, function(m) {
  
  sub_data_Bd <- adiv_long_Bd %>% filter(metric == m)
  fit_Bd <- aov(value ~ Temperature * Day, data = sub_data_Bd)
  tuk_res_Bd <- TukeyHSD(fit_Bd, "Temperature:Day")[["Temperature:Day"]][, "p adj"]
  letters_raw_Bd <- multcompLetters(tuk_res_Bd)$Letters
  data.frame(
    Group = names(letters_raw_Bd),
    Letter = as.character(letters_raw_Bd),
    metric = m,
    stringsAsFactors = FALSE)}) %>%
  separate(Group, into = c("Temperature", "Day"), sep = "[:\\-]") %>%
  mutate(metric = factor(metric, levels = metrics_list),
    Temperature = factor(Temperature, levels = c("15","18","20","25","29")),
    Day = factor(Day, levels = c("0","6","12","18")))
placement_Bd <- adiv_long_Bd %>%
  group_by(metric, Day, Temperature) %>%
  summarise(
    Placement.Value = max(value, na.rm = TRUE) * 1.05, .groups = "drop")

letters.df_Bd <- left_join(letters.df_Bd, placement_Bd, by = c("metric", "Temperature", "Day"))

Bd_Only_Adiv_plot <- adiv_long_Bd %>%
  ggplot(aes(x = Day, y = value, fill = Temperature)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    color = "black", 
    size = 1.5, 
    alpha = 0.5, 
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
  geom_text(
    data = letters.df_Bd,
    aes(x = Day, y = Placement.Value, label = Letter, group = Temperature),
    position = position_dodge(width = 1),
    vjust = -0.3,
    size = 5,
    fontface = "bold",
    inherit.aes = FALSE, 
    check_overlap = TRUE) +
  labs(title = "A", x = "", y = "") +
  facet_wrap(~ metric, ncol= 1, scales = "free") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(1, 'cm'),
    axis.text.x = element_text()) +
  scale_fill_manual(values = PalleteAdivGTDB, name = "Temperature")


##Bsal 
adiv_long_Bsal <- AdivBsal_Only %>%
  gather(key = metric, value = value, all_of(metrics_list)) %>%
  mutate(metric = factor(metric, levels = metrics_list),
         Temperature = factor(Temperature, levels = c("15","18","20","25","29")),
         Day = factor(Day, levels = c("0","6","12","18")))

letters.df_Bsal <- map_df(metrics_list, function(m) {
  sub_data_Bsal <- adiv_long_Bsal %>% filter(metric == m)
  fit_Bsal <- aov(value ~ Temperature * Day, data = sub_data_Bsal)
  tuk_res_Bsal <- TukeyHSD(fit_Bsal, "Temperature:Day")[["Temperature:Day"]][, "p adj"]
  letters_raw_Bsal <- multcompLetters(tuk_res_Bsal)$Letters
  data.frame(
    Group = names(letters_raw_Bsal),
    Letter = as.character(letters_raw_Bsal),
    metric = m,
    stringsAsFactors = FALSE)}) %>%
  separate(Group, into = c("Temperature", "Day"), sep = "[:\\-]") %>%
  mutate(metric = factor(metric, levels = metrics_list),
    Temperature = factor(Temperature, levels = c("15","18","20","25","29")),
    Day = factor(Day, levels = c("0","6","12","18")))
placement_Bsal <- adiv_long_Bsal %>%
  group_by(metric, Day, Temperature) %>%
  summarise(
    Placement.Value = max(value, na.rm = TRUE) * 1.05, .groups = "drop")

letters.df_Bsal <- left_join(letters.df_Bsal, placement_Bsal, by = c("metric", "Temperature", "Day"))

Bsal_Only_Adiv_plot <- adiv_long_Bsal %>%
  ggplot(aes(x = Day, y = value, fill = Temperature)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    color = "black", 
    size = 1.5, 
    alpha = 0.5, 
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
  geom_text(data = letters.df_Bsal,
    aes(x = Day, y = Placement.Value, label = Letter, group = Temperature),
    position = position_dodge(width = 1),
    vjust = -0.3,
    size = 5,
    fontface = "bold",
    inherit.aes = FALSE, 
    check_overlap = TRUE) +
  labs(title = "B", x = "", y = "") +
  facet_wrap(~ metric, ncol = 1, scales = "free") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(1, 'cm'),
    axis.text.x = element_text()) +
  scale_fill_manual(values = PalleteAdivGTDB, name = "Temperature")

Adiv_with_Letters<- Bd_Only_Adiv_plot + Bsal_Only_Adiv_plot + plot_layout(guides = "collect") & theme(legend.position = "right") 

Adiv_with_Letters<- Adiv_with_Letters & theme(text = element_text(size = 10))

#### Supplemental Table 2: Statistical results for two-way ANOVA on alpha diversity ####

##Bd 
Shannon_Bd_ANOVA_2<-aov(Shannon ~ Temperature * Day, data=AdivBd_Only)
summary(Shannon_Bd_ANOVA_2)

Chao1_Bd_ANOVA_2<- aov(Chao1 ~ Temperature * Day, data=AdivBd_Only)
summary(Chao1_Bd_ANOVA_2)

Evenness_Bd_ANOVA_2<-aov(Evenness ~ Temperature * Day, data=AdivBd_Only)
summary(Evenness_Bd_ANOVA_2)


## Bsal 

AdivBsal_Only$Day <- as.factor(AdivBsal_Only$Day)
AdivBsal_Only$Temperature <- as.factor(AdivBsal_Only$Temperature)

Shannon_Bsal_ANOVA_2<- aov(Shannon ~ Temperature * Day, data=AdivBsal_Only)
summary(Shannon_Bsal_ANOVA_2)

Chao1_Bsal_ANOVA_2<-aov(Chao1 ~ Temperature* Day, data=AdivBsal_Only)
summary(Chao1_Bsal_ANOVA_2)

Evenness_Bsal_ANOVA_2<-aov(Evenness ~ Temperature * Day, data=AdivBsal_Only)
summary(Evenness_Bsal_ANOVA_2)

####Supplemental Table 3: Results of post hoc testing of alpha diversity by day for SynComs containing Bsal ####


Shan_Tukey_Day_Bsal<-TukeyHSD(Shannon_Bsal_ANOVA_2, "Day")
print(Shan_Tukey_Day_Bsal)

Chao1_Tukey_Day_Bsal<-TukeyHSD(Chao1_Bsal_ANOVA_2, "Day")
print(Chao1_Tukey_Day_Bsal)

Even_Tukey_Day_Bsal<-TukeyHSD(Evenness_Bsal_ANOVA_2, "Day")
print(Even_Tukey_Day_Bsal)

