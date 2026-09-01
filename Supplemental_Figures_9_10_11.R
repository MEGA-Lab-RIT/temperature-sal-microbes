###Supplemental Figures 9, 10 & 11
###Created on: 4/05/26 
###Emma Thompson 
### Creating Supplemental Figure 9 :  Venn diagrams comparing differentially expressed genes between treatments. Differentially expressed genes between 15°C and 29°C at padj < 0.05 and logFC > 1 were compared between trials. Genes expressed within each experimental group and shared between groups represent a percentage of the total number of genes expressed during the trial.
# Figure 10: Comparison of differentially expressed genes in the Bd microcosms at 20 vs. 29°C. Heat map of the top 50 differentially expressed KOs (i.e., padj < 0.05 and logFC > 1) between the 20°C (thermal optimum) and 29°C treatments ordered by their KEGG Orthology identifier. Columns represent metatranscriptomic replicates, heatmap color denotes DESeq2-normalized gene counts, and gray scale identifies the functional group of each KO.
# Figure 11: Comparison of differentially expressed genes in the Bsal microcosms at 15 vs. 29°C. Heat map of the top 50 differentially expressed genes (i.e., padj < 0.05 and logFC > 1) between the 15°C (thermal optimum) and 29°C treatments ordered by their KEGG Orthology identifier. Columns represent metatranscriptome replicates, heatmap color denotes DESeq2-normalized gene counts, and gray scale identifies the functional group of each KO.

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

#### Normalized Bsal ####
normalized_counts_Bd <- counts(DeSeqDataset_Bd2, normalized = TRUE)

Sig_normalized_counts_Bd <- normalized_counts_Bd[rownames(Only_Sig_Results_Bd), ] #extracting only the significant genes row names in one match the other 

Sig_normalized_counts_Bd_df <- as.data.frame(Sig_normalized_counts_Bd)

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

#### Normalized Bsal ####
normalized_counts_Bsal <- counts(DeSeqDataset_Bsal2, normalized = TRUE)

Sig_normalized_counts_Bsal <- normalized_counts_Bsal[rownames(Only_Sig_Results_Bsal), ] #extracting only the significant genes row names in one match the other 

Sig_normalized_counts_Bsal_df <- as.data.frame(Sig_normalized_counts_Bsal)

#### Figure 9 ####

genes_Bd   <- rownames(Only_Sig_Results_Bd)
genes_Bsal <- rownames(Only_Sig_Results_Bsal)

gene_lists <- list(Bd = genes_Bd, Bsal = genes_Bsal)

genes_15C <- rownames(Sig_normalized_counts_Bsal_df)[
  rowSums(Sig_normalized_counts_Bsal_df[, c("Bsal_1_15", "Bsal_2_15", "Bsal_3_15", "Bsal_4_15")]) > 0]

genes_29C_Bsal <- rownames(Sig_normalized_counts_Bsal_df)[
  rowSums(Sig_normalized_counts_Bsal_df[, c("Bsal_1_29", "Bsal_2_29", "Bsal_4_29")]) > 0]

genes_20C_Bd <- rownames(Sig_normalized_counts_Bd_df)[
  rowSums(Sig_normalized_counts_Bd_df[, c("Bd_1_20", "Bd_2_20", "Bd_4_20", "Bd_5_20")]) > 0]

genes_29C_Bd <- rownames(Sig_normalized_counts_Bd_df)[
  rowSums(Sig_normalized_counts_Bd_df[, c("Bd_1_29", "Bd_2_29", "Bd_4_29", "Bd_5_29")]) > 0]

temp_gene_lists_Bd <- list(
  "20°C" = genes_20C_Bd,
  "29°C" = genes_29C_Bd)

Venn_Temp_Plot_Bd <- ggvenn(
  temp_gene_lists_Bd, 
  fill_color = c("#7c5295", "#280051"), # Cool purple for 20C, deep dark purple for 29C
  stroke_size = 0.5, 
  set_name_size = 5,
  text_size = 4) +
  coord_cartesian(expand = FALSE) +
  theme(plot.margin = margin(0, 0, 0, 0))

genes_15C <- rownames(Sig_normalized_counts_Bsal_df)[
  rowSums(Sig_normalized_counts_Bsal_df[, c("Bsal_1_15", "Bsal_2_15", "Bsal_3_15", "Bsal_4_15")]) > 0]

genes_29C_Bsal <- rownames(Sig_normalized_counts_Bsal_df)[
  rowSums(Sig_normalized_counts_Bsal_df[, c("Bsal_1_29", "Bsal_2_29", "Bsal_4_29")]) > 0]

temp_gene_lists_Bsal <- list(
  "15°C" = genes_15C,
  "29°C" = genes_29C_Bsal)

Venn_Temp_Plot_Bsal <- ggvenn(
  temp_gene_lists_Bsal, 
  fill_color = c("#4292c6", "#084594"), # Cool purple for 20C, deep dark purple for 29C
  stroke_size = 0.5, 
  set_name_size = 5,
  text_size = 4) +
  coord_cartesian(expand = FALSE) +
  theme(plot.margin = margin(0, 0, 0, 0))

Groups <- list(
  "Pathogen Optimums" = c(genes_20C_Bd, genes_15C),
  "29°C Trials"      = c(genes_29C_Bd, genes_29C_Bsal))

Venn_Temp_Groups <- ggvenn(
  Groups, 
  fill_color = c( "#78be21", "#008000"),
  stroke_size = 0.5, 
  set_name_size = 4,
  text_size = 4) + 
  #labs (subtitle = "Comparison of Bd vs. Bsal at 29°C")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
  coord_cartesian(expand = FALSE) +
  theme(plot.margin = margin(0, 0, 0, 0))

Supplemental_Venns<- (Venn_Temp_Groups / (Venn_Temp_Plot_Bd + Venn_Temp_Plot_Bsal)+ plot_layout(heights = c(2, 1))) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold"))

####Figure 10 #### 

normalized_counts_Bd_df2 <- Sig_normalized_counts_Bd_df %>%
  rownames_to_column(var = "gene")

top_50_genes_bd <- normalized_counts_Bd_df2 %>%
  mutate(total_count = rowSums(across(-gene), na.rm = TRUE)) %>%
  arrange(desc(total_count)) %>%
  slice(1:50) %>%
  select(-total_count)

normalized_counts_Bd_df3 <- top_50_genes_bd %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "count"
  ) %>%
  mutate(countfactor = cut(
    count,
    breaks = c(-1, 0, 1000, 5000, 10000, 25000, 50000, 75000, 100000, max(count, na.rm = TRUE)),
    labels = c("0","1000","1000-5000","5000-10000","10000-25000","25000-50000","50000-75000","75000-100000",">100000"))) %>%
  mutate(countfactor = factor(
    as.character(countfactor),
    levels = rev(c("0","1000","1000-5000","5000-10000","10000-25000","25000-50000","50000-75000","75000-100000",">100000"))))

normalized_counts_Bd_df3<- top_50_genes_bd %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "count") %>%
  mutate(log_count = log10(count + 1)) %>%
  mutate(countfactor = cut(log_count,
                           breaks = c(-0.1, 0, 1, 2, 3, 4, 5, max(log_count, na.rm = TRUE)),
                           labels = c("0","1", "1-2","2-3","3-4", "4-5", ">5"))) %>%
  mutate(countfactor = factor(
    as.character(countfactor),
    levels = rev(c("0","1", "1-2","2-3","3-4", "4-5", ">5"))))


normalized_counts_Bd_df3$sample <- factor(
  normalized_counts_Bd_df3$sample,
  levels = c("Bd_1_20","Bd_2_20","Bd_4_20","Bd_5_20","Bd_1_29","Bd_2_29","Bd_4_29", "Bd_5_29"))

Plain_HM_palette_Bd <- c( "#000", "#7f4886", "#c69ecb","#f5ebfa")

Heatmap_Bd_Plain<- ggplot(normalized_counts_Bd_df3, aes(x = sample, y = gene, fill = countfactor)) +
  geom_tile(color = "white", size = 0.25) +
  geom_vline(xintercept = 4.5, color = "white", size = 2) +
  labs(x = "", y = "", fill = "Gene Count") +
  scale_x_discrete(breaks = c("Bd_5_20", "Bd_1_29"),
                   labels = c("20°C                                   ", "                                   29°C")) +
  scale_fill_manual(values = Plain_HM_palette_Bd) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = ggtext::element_markdown(size = 10),
        legend.position = "right")

All_KOs_Function <- read.csv("Regular_HM_KOs_Function.csv")

All_KOs_Function <- All_KOs_Function %>% 
  rename(gene = KO.Number,)

KO_anot_Bd <- left_join(normalized_counts_Bd_df3, All_KOs_Function, by = "gene", relationship = "many-to-many")%>% 
  arrange(Function, gene)

ordered_genes_Bd <- unique(KO_anot_Bd$gene)

KO_anot_Bd$gene <- factor(KO_anot_Bd$gene, levels = ordered_genes_Bd)
normalized_counts_Bd_df3$gene <- factor(normalized_counts_Bd_df3$gene, levels = ordered_genes_Bd)

annotation_pallet_reg<- c("#FFFFFF", "#F2F2F2","#E6E6E6", "#D9D9D9", "#CCCCCC","#BFBFBF" ,"#B3B3B3", "#A6a6a6", "#999999", "#8c8c8c", "#808080" , "#737373", "#666666" , "#595959" , "#4d4d4d" ,"#333333" , "#262626" , "#1a1a1a")

annotation_KO_Bd <- ggplot(KO_anot_Bd, aes(y = gene, x = 1, fill = Function)) +
  geom_tile() +
  scale_fill_manual(values = annotation_pallet_reg, name = "Function", labels = function(x) str_wrap(x, width = 60)) +
  theme_void() +
  guides(fill = guide_legend(reverse = TRUE))+
  theme(
    legend.position = "left",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12))

Bd_HM_w_Annotations <-annotation_KO_Bd + Heatmap_Bd_Plain + plot_layout(widths = c(0.1, 1))

####Figure 11 ####
normalized_counts_Bsal_df2 <- Sig_normalized_counts_Bsal_df %>%
  rownames_to_column(var = "gene")

top_50_genes <- normalized_counts_Bsal_df2 %>%
  mutate(total_count = rowSums(across(-gene), na.rm = TRUE)) %>%
  arrange(desc(total_count)) %>%
  slice(1:50) %>%
  select(-total_count)

normalized_counts_Bsal_df3 <- top_50_genes %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "count"
  ) %>%
  mutate(countfactor = cut(
    count,
    breaks = c(-1, 0, 1000, 5000, 10000, 25000, 50000, 75000, 100000, max(count, na.rm = TRUE)),
    labels = c("0","1000","1000-5000","5000-10000","10000-25000","25000-50000","50000-75000","75000-100000",">100000"))) %>%
  mutate(countfactor = factor(
    as.character(countfactor),
    levels = rev(c("0","1000","1000-5000","5000-10000","10000-25000","25000-50000","50000-75000","75000-100000",">100000"))))


normalized_counts_Bsal_df3$sample <- factor(
  normalized_counts_Bsal_df3$sample,
  levels = c("Bsal_1_15","Bsal_2_15","Bsal_3_15","Bsal_4_15","Bsal_1_29","Bsal_2_29","Bsal_4_29"))


normalized_counts_Bsal_df3_v2<- top_50_genes %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "count") %>%
  mutate(log_count = log10(count + 1)) %>%
  mutate(countfactor = cut(log_count,
                           breaks = c(-0.1, 0, 1, 2, 3, 4, 5, max(log_count, na.rm = TRUE)),
                           labels = c("0","1", "1-2","2-3","3-4", "4-5", ">5"))) %>%
  mutate(countfactor = factor(
    as.character(countfactor),
    levels = rev(c("0","1", "1-2","2-3","3-4", "4-5", ">5"))))


normalized_counts_Bsal_df3_v2$sample <- factor(
  normalized_counts_Bsal_df3_v2$sample,
  levels = c("Bsal_1_15","Bsal_2_15","Bsal_3_15","Bsal_4_15","Bsal_1_29","Bsal_2_29","Bsal_4_29"))

Plain_HM_palette_Bsal <- c( "#000", "#284060", "#6e8ebf",  "#ebf2fa")


Heatmap_Bsal_Plain<- ggplot(normalized_counts_Bsal_df3_v2, aes(x = sample, y = gene, fill = countfactor)) +
  geom_tile(color = "white", size = 0.25) +
  geom_vline(xintercept = 4.5, color = "white", size = 2) +
  labs(x = "", y = "", fill = "Gene Count") +
  scale_x_discrete(breaks = c("Bsal_2_15", "Bsal_2_29"),
                   labels = c("15", "29")) +
  scale_fill_manual(values = Plain_HM_palette_Bsal) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = ggtext::element_markdown(size = 10),
        legend.position = "right")

KO_anot_Bsal <- left_join(normalized_counts_Bsal_df3_v2, All_KOs_Function, by = "gene", relationship = "many-to-many")%>% 
  arrange(Function, gene)

ordered_genes_Bsal <- unique(KO_anot_Bsal$gene)

KO_anot_Bsal$gene <- factor(KO_anot_Bsal$gene, levels = ordered_genes_Bd)
normalized_counts_Bsal_df3_v2$gene <- factor(normalized_counts_Bsal_df3_v2$gene, levels = ordered_genes_Bd)

annotation_pallet_reg<- c("#FFFFFF", "#F2F2F2","#E6E6E6", "#D9D9D9", "#CCCCCC","#BFBFBF" ,"#B3B3B3", "#A6a6a6", "#999999", "#8c8c8c", "#808080" , "#737373", "#666666" , "#595959" , "#4d4d4d" ,"#333333" , "#262626" , "#1a1a1a")

annotation_KO_Bsal <- ggplot(KO_anot_Bsal, aes(y = gene, x = 1, fill = Function)) +
  geom_tile() +
  scale_fill_manual(values = annotation_pallet_reg, name = "Function", labels = function(x) str_wrap(x, width = 60)) +
  theme_void() +
  guides(fill = guide_legend(reverse = TRUE))+
  theme(
    legend.position = "left",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12))

Bsal_HM_w_Annotations <-annotation_KO_Bsal + Heatmap_Bsal_Plain + plot_layout(widths = c(0.1, 1))


