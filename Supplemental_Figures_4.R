#### Supplemental_Figures_4&6####
###Created on: 3/10/26   
###Emma Thompson 
# Creating Supplemental Figure 4:Log abundances of bacterial families over time and across temperature treatments. Samples are ordered by time within each temperature subset. A SynComs containing Bd, B SynComs containing Bsal.
#Figure 6: Comparison of bacterial relative abundance using either the Genome Taxonomy Database (GTDB) or SILVA taxonomy. A-B: GTDB, C-D: SILVA, with A & C representing SynComs containing Bd and B & D representing SynComs containing Bsal. Color depicts unique taxonomic families within the synthetic communities (all matched to inoculated strains). Notable differences include the abundance of Sphingobacterium and the other category (families with less than 3.75% relative abundance within the microcosm) in the SILVA-identified reads. Also, the SILVA-identified reads showed an additional family, Oxalobacteriaceae, which likely contains many reads that were identified as Burkholderiaceae by GTDB. With the exception of classifiers, all parameters in data processing remained the same. GTDB v2.6.0 and SILVA v138.2 were used for taxonomic classification. 

####Packages ####
library(ggplot2)
library(vegan)
library(dplyr)
library(phyloseq)
library(ape)
library(tidyr)
library(microbiome)
library(stringr)
library(patchwork)
library(agricolae)
library(ggpubr)
library(ggrepel)
library(tidyverse)
library(ANCOMBC)

####Figure 4 ####

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

physeq_filtered <- subset_samples(physeq16s, Contents != "Synth")

physeq_filtered1 <- subset_samples(physeq_filtered, Contents != "Bd")

physeq_filtered2 <- subset_samples(physeq_filtered1, Contents != "Bsal")

physeq_filtered3 <- subset_samples(physeq_filtered2, Contents != "Ext_cntrl")

physeq_filtered4 <- subset_samples(physeq_filtered3, Contents != "PCR_cntrl")

physeq_filtered5 <- subset_samples(physeq_filtered4, Contents != "Media")

####Pruning Low Abundance Taxa ####

taxa_to_keep <- taxa_sums(physeq_filtered5) > 500
ps_pruned <- prune_taxa(taxa_to_keep, physeq_filtered5)

#### Trimming ####

physeq.trim <- subset_taxa(ps_pruned, !is.na(Class) &
                             !Class %in% c("", "Unknown", "uncharacterized",
                                           "unidentified")) 
physeq.trim <- subset_taxa(physeq.trim, !is.na(Order) &
                             !Order %in% c("", "Unknown", "uncharacterized",
                                           "unidentified")) 

physeq.trim <- subset_taxa(physeq.trim, !is.na(Family) &
                             !Family %in% c("", "Unknown", "uncharacterized",
                                            "unidentified")) 


####Subset ####
T15 <- subset_samples(physeq.trim, Temperature %in% c("15", NA))
T18 <- subset_samples(physeq.trim, Temperature %in% c("18", NA))
T20 <- subset_samples(physeq.trim, Temperature %in% c("20", NA))
T25 <- subset_samples(physeq.trim, Temperature %in% c("25", NA))
T29 <- subset_samples(physeq.trim, Temperature %in% c("29", NA))

####Subest by Contents ####

physeq.trim_Bd <- subset_samples(physeq.trim, Contents != "Bsal_Synth")

metadata_df_Bd <- data.frame(sample_data(physeq.trim_Bd))

physeq.trim_Bsal <- subset_samples(physeq.trim, Contents != "Bd_Synth")

metadata_df_Bsal<- data.frame(sample_data(physeq.trim_Bsal))

####ANCOM_BC Temperature ####

sample_data(physeq.trim_Bsal)$Temperature <- as.factor(sample_data(physeq.trim_Bsal)$Temperature)

set.seed(123)
Temp_ANCOMBC_Bsal <- ancombc2(data = physeq.trim_Bsal, tax_level = "Family",
                              fix_formula = "Temperature", rand_formula = NULL,
                              p_adj_method = "holm", # authors of ANCOMBC recommend holm correction over BH
                              pseudo_sens = TRUE,  # pseudo_sens: whether to perform the sensitivity analysis; TRUE is default
                              prv_cut = 0.20, # Taxa with prevalence (the proportion of samples in which the taxon is present) less than prv_cut will be excluded; 0.10 is default
                              lib_cut = 1000, # Samples with library sizes less than lib_cut will be excluded in the analysis; 0 is default
                              group = "Temperature",
                              alpha = 0.001,  # set your alpha (p-value cut-off)
                              verbose = TRUE,
                              global = FALSE, # global: discern taxa that demonstrate differential abundance between a MINIMUM of two groups when analyzing three or more experimental groups
                              pairwise = TRUE, # pairwise: designed to identify taxa that exhibit differential abundance between ANY two groups within a set of three or more experimental groups
                              dunnet = FALSE, # Dunnet: multiple pairwise comparisons against a pre-specified group (e.g., control or reference group)
                              trend = FALSE # trend: when you expect results to align with specific patterns (e.g., dose-response; monotonically increasing, decreasing, or umbrella shaped)
)


# view output
#res_prim_Bsal = Temp_ANCOMBC_Bsal$res # if you only have two groups
res_pair_Bsal= Temp_ANCOMBC_Bsal$res_pair # if you have more than two groups and want the result for all pairwise comparisons
feature_Bsal_Temp= Temp_ANCOMBC_Bsal$feature_table

####Bsal Temp Heatmap ####

feature_Bsal_Temp_df <- as.data.frame(t(feature_Bsal_Temp))

feature_Bsal_Temp_df$sample_id <- rownames(feature_Bsal_Temp_df)

metadata_df_Bsal$sample_id <- rownames(metadata_df_Bsal)

merged_df_Bsal_Temp <- merge(metadata_df_Bsal, feature_Bsal_Temp_df, by = "sample_id")

long_df_Bsal_Temp <- merged_df_Bsal_Temp %>%
  rownames_to_column(var = "SampleID") %>%
  pivot_longer(
    cols = -(SampleID:Replicate),
    names_to = "Taxon",
    values_to = "Abundance"
  ) %>%
  filter(Taxon != "Microbacteriaceae")

taxon_order_Bsal_Time <- long_df_Bsal_Temp %>%
  group_by(Taxon) %>%
  summarise(Total_Abundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(Total_Abundance)) %>%
  pull(Taxon)

long_df_Bsal_Temp$Taxon <- factor(long_df_Bsal_Temp$Taxon, levels = taxon_order_Bsal_Time)

sample_info_Bsal_Temp <- long_df_Bsal_Temp %>%
  distinct(SampleID, Temperature, Day) %>%
  arrange(Temperature, Day, SampleID) %>%
  mutate(SampleOrder = row_number())


long_df_Bsal_Temp <- long_df_Bsal_Temp %>%
  left_join(sample_info_Bsal_Temp, by = c("SampleID", "Temperature", "Day"))

temp_labels_Bsal_Temp <- sample_info_Bsal_Temp %>%
  group_by(Temperature) %>%
  summarise(
    xmin = min(SampleOrder),
    xmax = max(SampleOrder),
    xmid = mean(c(xmin, xmax)))


HeatMap_Bsal_Temp<-ggplot(long_df_Bsal_Temp, aes(x = SampleOrder, y = Taxon, fill = log(Abundance+1))) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#00008B") +
  scale_x_continuous(
    breaks = temp_labels_Bsal_Temp$xmid,
    labels = paste(temp_labels_Bsal_Temp$Temperature)) +
  labs(
    x = "Temperature ",
    y = "Family",
    fill = "log(Abundance +1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 15)) +
  geom_vline(
    xintercept = temp_labels_Bsal_Temp$xmax + 0.5,
    color = "black",
    linewidth = 0.3)

####Temperature Bd ####
sample_data(physeq.trim_Bd)$Temperature <- relevel(
  as.factor(sample_data(physeq.trim_Bd)$Temperature),
  ref = "20"
)

set.seed(123)
Temp_ANCOMBC_Bd <- ancombc2(data = physeq.trim_Bd, tax_level = "Family",
                            fix_formula = "Temperature", rand_formula = NULL,
                            p_adj_method = "holm", # authors of ANCOMBC recommend holm correction over BH
                            pseudo_sens = TRUE,  # pseudo_sens: whether to perform the sensitivity analysis; TRUE is default
                            prv_cut = 0.20, # Taxa with prevalence (the proportion of samples in which the taxon is present) less than prv_cut will be excluded; 0.10 is default
                            lib_cut = 1000, # Samples with library sizes less than lib_cut will be excluded in the analysis; 0 is default
                            group = "Temperature",
                            alpha = 0.001,  # set your alpha (p-value cut-off)
                            verbose = TRUE,
                            global = FALSE, # global: discern taxa that demonstrate differential abundance between a MINIMUM of two groups when analyzing three or more experimental groups
                            pairwise = TRUE, # pairwise: designed to identify taxa that exhibit differential abundance between ANY two groups within a set of three or more experimental groups
                            dunnet = FALSE, # Dunnet: multiple pairwise comparisons against a pre-specified group (e.g., control or reference group)
                            trend = FALSE # trend: when you expect results to align with specific patterns (e.g., dose-response; monotonically increasing, decreasing, or umbrella shaped)
                            
)

res_prim_Bd= Temp_ANCOMBC_Bd$res # if you only have two groups
res_pair_Bd= Temp_ANCOMBC_Bd$res_pair # if you have more than two groups and want the result for all pairwise comparisons
Temp_ANCOMBC_Bd$feature_table

#### Bd Heatmap Temp ####

feature_Bd_Temp= Temp_ANCOMBC_Bd$feature_table

feature_Bd_Temp_df <- as.data.frame(t(feature_Bd_Temp))

feature_Bd_Temp_df$sample_id <- rownames(feature_Bd_Temp_df)

metadata_df_Bd$sample_id <- rownames(metadata_df_Bd)

merged_df_Bd_Temp <- merge(metadata_df_Bd, feature_Bd_Temp_df, by = "sample_id")

long_df_Bd_Temp <- merged_df_Bd_Temp %>%
  rownames_to_column(var = "SampleID") %>%
  pivot_longer(
    cols = -(SampleID:Replicate),   # keep metadata columns
    names_to = "Taxon",
    values_to = "Abundance"
  )%>%
  filter(Taxon != "Microbacteriaceae")

taxon_order_Bd_Time <- long_df_Bd_Temp %>%
  group_by(Taxon) %>%
  summarise(Total_Abundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(Total_Abundance)) %>%
  pull(Taxon)

long_df_Bd_Temp$Taxon <- factor(long_df_Bd_Temp$Taxon, levels = taxon_order_Bd_Time)

sample_info_Bd_Temp <- long_df_Bd_Temp %>%
  distinct(SampleID, Temperature, Day) %>%
  arrange(Temperature, Day, SampleID) %>%
  mutate(SampleOrder = row_number())


long_df_Bd_Temp <- long_df_Bd_Temp %>%
  left_join(sample_info_Bd_Temp, by = c("SampleID", "Temperature", "Day"))

temp_labels_Bd_Temp <- sample_info_Bd_Temp %>%
  group_by(Temperature) %>%
  summarise(
    xmin = min(SampleOrder),
    xmax = max(SampleOrder),
    xmid = mean(c(xmin, xmax)))


HeatMap_Bd_Temp<-ggplot(long_df_Bd_Temp, aes(x = SampleOrder, y = Taxon, fill = log(Abundance+1))) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#00008B") +
  scale_x_continuous(
    breaks = temp_labels_Bd_Temp$xmid,
    labels = paste(temp_labels_Bd_Temp$Temperature)) +
  labs(
    x = "Temperature ",
    y = "Family",
    fill = "log(Abundance +1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 15)) +
  geom_vline(
    xintercept = temp_labels_Bd_Temp$xmax + 0.5,
    color = "black",
    linewidth = 0.3)

#### TEMP PLOTS TOGETHER ####

ANCOMBC_Results<- HeatMap_Bd_Temp / HeatMap_Bsal_Temp + plot_layout(guides = "collect") &
  theme(legend.position = "right")
ANCOMBC_w_Tags<-ANCOMBC_Results + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = 'bold'))


