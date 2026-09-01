###Figure_4
###Created on: 4/05/26 
###Emma Thompson 
### Creating Figure 4: Differences in the metatranscriptomic profile of microcosms with pathogen and temperature. A) Venn diagram comparing the difference in the total number of unique transcripts identified within each experimental group. Genes expressed within each experimental group and shared between groups represent a percentage of the total number of genes expressed during the trial. PCoA of gene expression of microcosms containing Bd (B) and Bsal (C). No ellipse is shown for the 29°C trial with Bsal since n < 4 due to poor RNA quality of one replicate.

####Packages####
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(readxl)
library(vegan)
library(patchwork)
library(ggtext)
library(cowplot)
library(ggpubr)
library(ggVennDiagram)
library(ggvenn)

#### Reading in Bd Transcripts ####
Bd_annotations<-read.csv("Merged_No_CDS_Bd_KO.csv", header = TRUE)

Sample_Names_Bd<- read.csv("Sample_Names_Bd.csv", header = TRUE)

####Formating Bd Dataframe ####
Bd_annotations_Unique_KOs <- Bd_annotations %>%
  group_by(KO) %>%
  summarize(across(everything(), sum)) %>% #Across all columns sum counts of all the same KOs 
  ungroup()

Bd_annotations_Unique_KOs_df <- as.data.frame(Bd_annotations_Unique_KOs) 

rownames(Bd_annotations_Unique_KOs_df) <- Bd_annotations_Unique_KOs_df[, 1] 

Bd_annotations_Unique_KOs_df$KO <- NULL  

Sample_Names_Bd<- Sample_Names_Bd[, -c(4,5,6,7)] 
rownames(Sample_Names_Bd) <- Sample_Names_Bd[, 1] 

Sample_Names_Bd$Temperature <- as.factor(Sample_Names_Bd$Temperature) 

####Setting up DESeq2 For Bd Samples####

DeSeqDataset_Bd <- DESeqDataSetFromMatrix(
  countData = Bd_annotations_Unique_KOs_df,
  colData = Sample_Names_Bd,
  design = ~ Temperature)

DeSeqDataset_Bd2 <- DeSeqDataset_Bd[rowSums(counts(DeSeqDataset_Bd)) > 10, ] 

DeSeqDataset_Bd2 <- DESeq(DeSeqDataset_Bd2)

Results_Bd <- results(DeSeqDataset_Bd2, contrast = c("Temperature", "29", "20"))

Results_Bd_df <- as.data.frame(Results_Bd) 

Only_Sig_Results_Bd <- Results_Bd_df[which(Results_Bd_df$padj < 0.05 & abs(Results_Bd_df$log2FoldChange) > 1), ]

####Data Visualization Bd####

MA_Bd<-plotMA(Results_Bd, ylim = c(-5,5))


vsd_Bd<- vst(DeSeqDataset_Bd2, blind = FALSE)
Bd_PCA<-plotPCA(vsd_Bd, intgroup = "Temperature")

vsd_Bd_df <- assay(vsd_Bd) 
dist_matrix_vsd_Bd <- dist(t(vsd_Bd_df)) 
group_factor_from_Metadata_Bd <- colData(DeSeqDataset_Bd)$Temperature 
model_betdisp_Bd <- betadisper(dist_matrix_vsd_Bd, group_factor_from_Metadata_Bd) 
anova(model_betdisp_Bd)
BetaDisp_PCoA_Bd<-plot(model_betdisp_Bd) 

Bd_points <- data.frame(model_betdisp_Bd$vectors[,1:2])
Bd_points$group <- model_betdisp_Bd$group
Bd_centroids <- data.frame(model_betdisp_Bd$centroids[,1:2])
Bd_centroids$group <- rownames(Bd_centroids)
Bd_segments <- Bd_points %>%
  rename(x = PCoA1, y = PCoA2) %>%
  mutate(xend = Bd_centroids$PCoA1[match(group, Bd_centroids$group)],
         yend = Bd_centroids$PCoA2[match(group, Bd_centroids$group)])
Bd_points$Temperature <- Sample_Names_Bd$Temperature
Bd_centroids$Temperature <- c("15", "29")


PCoA_Bd <- ggplot() +
  geom_segment(data = Bd_segments,aes(x = x, y = y, xend = xend, yend = yend, color = group), alpha = 0.5) +
  geom_point(data = Bd_points, aes(x = PCoA1, y = PCoA2, color = group), size = 2) +
  geom_point(data = Bd_centroids,aes(x = PCoA1, y = PCoA2, fill = group), size = 3, shape = 21, color = "black") +
  stat_ellipse(data = Bd_points, aes(x = PCoA1, y = PCoA2, color = group)) +
  scale_color_manual(name = "Temperature",
                     values = c("29" = "#ffa1e7", "20" = "#96dcf8"),
                     labels = c("20°C", "29°C")) +
  scale_fill_manual(name = "Temperature",
                    values = c("29" = "#ffa1e7", "20" = "#96dcf8"),
                    labels = c("20°C", "29°C")) +
  labs( x = paste0("PCoA1 (", round(100 * model_betdisp_Bd$eig[1] /sum(model_betdisp_Bd$eig), 1), "%)"),
        y = paste0("PCoA2 (", round(100 * model_betdisp_Bd$eig[2] /sum(model_betdisp_Bd$eig), 1), "%)")) +
  theme_classic() +
  coord_fixed()

#### Reading in Bsal data ####

Bsal_annotations<-read.csv("Merged_No_CDS_Bsal_COUNTS.csv", header = TRUE)

Sample_Names_Bsal<- read.csv("Bsal_Sample_Names.csv", header = TRUE)

####Formating Dataframe Bsal####
Bsal_annotations<- Bsal_annotations[, -c(9)] 
Bsal_annotations <- Bsal_annotations %>% 
  filter(!is.na(KO) & KO!= "")

Bsal_annotations_Unique_KOs <- Bsal_annotations %>%
  group_by(KO) %>%
  summarize(across(everything(), sum)) %>% 
  ungroup()

Bsal_annotations_Unique_KOs_df <- as.data.frame(Bsal_annotations_Unique_KOs) 

rownames(Bsal_annotations_Unique_KOs_df) <- Bsal_annotations_Unique_KOs_df[, 1] 

Bsal_annotations_Unique_KOs_df$KO <- NULL 

Sample_Names_Bsal<- Sample_Names_Bsal[, -c(4,5,6,7)] 
rownames(Sample_Names_Bsal) <- Sample_Names_Bsal[, 1] 

Sample_Names_Bsal$Temperature <- as.factor(Sample_Names_Bsal$Temperature) 

####Setting up DESeq2 Bsal####

DeSeqDataset_Bsal <- DESeqDataSetFromMatrix(
  countData = Bsal_annotations_Unique_KOs_df,
  colData = Sample_Names_Bsal,
  design = ~ Temperature)

DeSeqDataset_Bsal2 <- DeSeqDataset_Bsal[rowSums(counts(DeSeqDataset_Bsal)) > 10, ] 

DeSeqDataset_Bsal2 <- DESeq(DeSeqDataset_Bsal2) 

Results_Bsal <- results(DeSeqDataset_Bsal2, contrast = c("Temperature", "29", "15"))

Results_Bsal_df <- as.data.frame(Results_Bsal) 

Only_Sig_Results_Bsal <- Results_Bsal_df[which(Results_Bsal_df$padj < 0.05 & abs(Results_Bsal_df$log2FoldChange) > 1), ]


####Data Visualization Bsal####

MA_Bsal<-plotMA(Results_Bsal, ylim = c(-5,5))


vsd_Bsal <- vst(DeSeqDataset_Bsal2, blind = FALSE)
Bsal_PCA<-plotPCA(vsd_Bsal, intgroup = "Temperature")

vsd_Bsal_df <- assay(vsd_Bsal) 
dist_matrix_vsd_Bsal <- dist(t(vsd_Bsal_df)) 
group_factor_from_Metadata <- colData(DeSeqDataset_Bsal2)$Temperature 
model_betdisp_Bsal <- betadisper(dist_matrix_vsd_Bsal, group_factor_from_Metadata) 
anova(model_betdisp_Bsal) 
BetaDisp_PCoA_Bsal<-plot(model_betdisp_Bsal) 

###Beta dispersion plot 
Bsal_points <- data.frame(model_betdisp_Bsal$vectors[,1:2])
Bsal_points$group <- model_betdisp_Bsal$group
Bsal_centroids <- data.frame(model_betdisp_Bsal$centroids[,1:2])
Bsal_centroids$group <- rownames(Bsal_centroids)
Bsal_segments <- Bsal_points %>%
  rename(x = PCoA1, y = PCoA2) %>%
  mutate(xend = Bsal_centroids$PCoA1[match(group, Bsal_centroids$group)],
         yend = Bsal_centroids$PCoA2[match(group, Bsal_centroids$group)])
Bsal_points$Temperature <- Sample_Names_Bsal$Temperature
Bsal_centroids$Temperature <- c("15", "29")


PCoA_Bsal <- ggplot() +
  geom_segment(data = Bsal_segments,aes(x = x, y = y, xend = xend, yend = yend, color = group), alpha = 0.5) +
  geom_point(data = Bsal_points, aes(x = PCoA1, y = PCoA2, color = group), size = 2) +
  geom_point(data = Bsal_centroids, aes(x = PCoA1, y = PCoA2, fill = group), size = 3, shape = 21, color = "black") +
  stat_ellipse(data = Bsal_points, aes(x = PCoA1, y = PCoA2, color = group)) +
  scale_color_manual(name = "Temperature",
                     values = c("29" = "#ffa1e7", "15" = "#f2aa84"),
                     labels = c("15°C", "29°C")) +
  scale_fill_manual(name = "Temperature",
                    values = c("29" = "#ffa1e7", "15" = "#f2aa84"),
                    labels = c("15°C", "29°C")) +
  labs(x = paste0("PCoA1 (", round(100 * model_betdisp_Bsal$eig[1] /sum(model_betdisp_Bsal$eig), 1), "%)"),
       y = paste0("PCoA2 (", round(100 * model_betdisp_Bsal$eig[2] / sum(model_betdisp_Bsal$eig), 1), "%)")) +
  theme_classic() +
  coord_fixed()

#### Venn Diagrams ####

genes_Bd   <- rownames(Only_Sig_Results_Bd)
genes_Bsal <- rownames(Only_Sig_Results_Bsal)

genes_15C <- rownames(Sig_normalized_counts_Bsal_df)[
  rowSums(Sig_normalized_counts_Bsal_df[, c("Bsal_1_15", "Bsal_2_15", "Bsal_3_15", "Bsal_4_15")]) > 0]

genes_29C_Bsal <- rownames(Sig_normalized_counts_Bsal_df)[
  rowSums(Sig_normalized_counts_Bsal_df[, c("Bsal_1_29", "Bsal_2_29", "Bsal_4_29")]) > 0]

four_way_gene_list <- list(
  "Bd 20°C"   = genes_20C_Bd,
  "Bd 29°C"   = genes_29C_Bd,
  "Bsal 15°C" = genes_15C,
  "Bsal 29°C" = genes_29C_Bsal)

Four_Way_Venn <- ggvenn(
  four_way_gene_list,
  fill_color = c("#7c5295", "#280051", "#4292c6", "#084594"), 
  stroke_size = 0.5,
  set_name_size = 4,   
  text_size = 3.5 ) +
  coord_cartesian(expand = FALSE) +
  theme(plot.margin = margin(0, 0, 0, 0))

Figure_4<- Four_Way_Venn / (PCoA_Bd + PCoA_Bsal) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold"))

####Supplemental Table 10: Statistical results for PERMANOVA on Gene Expression ####

##Bd
vsd_Bd<- vst(DeSeqDataset_Bd2, blind = FALSE)

vsd_mat_Bd <- t(assay(vsd_Bd))

sample_dist_Bd <- vegdist(vsd_mat_Bd, method = "euclidean")

sample_meta_Bd <- as.data.frame(colData(DeSeqDataset_Bd2))

permanova_res_Bd <- adonis2(sample_dist_Bd ~ Temperature, data = sample_meta_Bd, permutations = 999)

##Bsal 
vsd_Bsal<- vst(DeSeqDataset_Bsal2, blind = FALSE)

vsd_mat_Bsal <- t(assay(vsd_Bsal))

sample_dist_Bsal <- vegdist(vsd_mat_Bsal, method = "euclidean")

sample_meta_Bsal <- as.data.frame(colData(DeSeqDataset_Bsal2))

permanova_res_Bsal <- adonis2(sample_dist_Bsal ~ Temperature, data = sample_meta_Bsal, permutations = 999)



