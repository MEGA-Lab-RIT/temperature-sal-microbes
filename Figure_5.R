###Figure_5
###Created on: 4/05/26 
###Emma Thompson 
### Creating Figure 5: Expression patterns of genes of interest identified via metatranscriptomic sequencing. Differentially expressed genes within microcosms with Bd (A) or Bsal (B) were filtered by KO number to represent only those that possess functions associated with fungal resistance, mucin utilization, or heat stress/heat shock. Gene counts are depicted as both log-transformed and DESeq2-normalized. 


####Packages####
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(readxl)
library(vegan)
library(patchwork)
library(ggtext)

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

####Extracting Counts Bd####
normalized_counts_Bd <- counts(DeSeqDataset_Bd2, normalized = TRUE) 

Sig_normalized_counts_Bd <- normalized_counts_Bd[rownames(Only_Sig_Results_Bd), ] 

Sig_normalized_counts_Bd_df <- as.data.frame(Sig_normalized_counts_Bd)

#### Visualization Bd####
normalized_counts_Bd_df2 <- Sig_normalized_counts_Bd_df %>%
  rownames_to_column(var = "gene")

KOs_OI_Function <- read.csv("KO_To_Function.csv")
kos_of_interest <- KOs_OI_Function$KO_Number

GOI_Bd<- normalized_counts_Bd_df2 %>%
  filter(gene %in% kos_of_interest)

GOI_Bd2<- GOI_Bd %>%
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

GOI_Bd2$sample <- factor(
  GOI_Bd2$sample,
  levels = c("Bd_1_20","Bd_2_20","Bd_4_20","Bd_5_20","Bd_1_29","Bd_2_29","Bd_4_29", "Bd_5_29"))

All_KOs_Function <- read.csv("KO_to_Function.csv")

All_KOs_Function <- rename(All_KOs_Function, gene = KO_Number)

KO_Bd_GOI_anot <- left_join(GOI_Bd2, All_KOs_Function, by = "gene", relationship = "many-to-many") %>% 
  arrange(Function, gene)

ordered_genes <- unique(KO_Bd_GOI_anot$gene)

KO_Bd_GOI_anot$gene <- factor(KO_Bd_GOI_anot$gene, levels = ordered_genes)
GOI_Bd2$gene <- factor(GOI_Bd2$gene, levels = ordered_genes)

palette_Bd <- c( "#000", "#402443", "#5f3665", "#7f4886","#b27db9", "#c69ecb","#f5ebfa")

Heatmap_Bd_GOI2<- ggplot(GOI_Bd2, aes(x = sample, y = gene, fill = countfactor)) +
  geom_tile(color = "white", size = 0.25) +
  geom_vline(xintercept = 4.5, color = "white", size = 2) +
  labs(x = "", y = "", fill = "Gene Count", title = "A") + 
  scale_x_discrete( breaks = c("Bd_5_20", "Bd_1_29"),
                    labels = c("20°C", "29°C")) + 
  scale_fill_manual(values = palette_Bd) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, vjust= 0.5, hjust = 0.5, size=15),
        axis.text.y = element_text(size = 8),
        legend.position = "right", 
        legend.text = element_text(size = 12))


subheader_pallete_gray<- c("#f3f3f3", "#D3D3D3","#B6b6b6", "#9E9E9E","#69686D","#4F4E52","#353437", "#000000")



annotation_Bd_KOI <- ggplot(KO_Bd_GOI_anot, aes(y = gene, x = 1, fill = Function)) +
  geom_tile() +
  scale_fill_manual(values = subheader_pallete_gray, name = "Function", labels = function(x) str_wrap(x, width = 60), guide = guide_legend(reverse = TRUE)) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12))


GOI_w_Strip_Bd <-annotation_Bd_KOI + Heatmap_Bd_GOI2 + plot_layout(widths = c(0.1, 1))

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

####Visualization Bsal ####
GOI_Bsal<- normalized_counts_Bsal_df2 %>%
  filter(gene %in% kos_of_interest)

GOI_Bsal2<- GOI_Bsal %>%
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


GOI_Bsal2$sample <- factor(
  GOI_Bsal2$sample,
  levels = c("Bsal_1_15","Bsal_2_15","Bsal_3_15","Bsal_4_15","Bsal_1_29","Bsal_2_29","Bsal_4_29"))

KO_Bsal_GOI_anot <- left_join(GOI_Bsal2, All_KOs_Function, by = "gene", relationship = "many-to-many")%>% 
  arrange(Function, gene)

ordered_genes <- unique(KO_Bsal_GOI_anot$gene)

KO_Bsal_GOI_anot$gene <- factor(KO_Bsal_GOI_anot$gene, levels = ordered_genes)
GOI_Bsal2$gene <- factor(GOI_Bsal2$gene, levels = ordered_genes)

palette2 <- c("#000000", "#1a2a40", "#284060", "#3a5d8c", "#6e8ebf", "#9cb4d9", "#ebf2fa", "#fff")

Heatmap_Bsal_GOI2<- ggplot(KO_Bsal_GOI_anot, aes(x = sample, y = gene, fill = countfactor)) +
  geom_tile(color = "white", size = 0.25) +
  geom_vline(xintercept = 4.5, color = "white", size = 2) +
  labs(x = "", y = "", fill = "Gene Count", title = "B") +
  scale_x_discrete( breaks = c("Bsal_2_15", "Bsal_2_29"),
                    labels = c("15°C", "29°C")) +
  scale_fill_manual(values = palette2) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), 
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size=15),
        axis.text.y = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 12))


annotation_Bsal_KOI <- ggplot(KO_Bsal_GOI_anot, aes(y = gene, x = 1, fill = Function)) +
  geom_tile() +
  scale_fill_manual(values = subheader_pallete_gray, name = "Function", labels = function(x) str_wrap(x, width = 60)) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12))


GOI_w_Strip <-annotation_Bsal_KOI + Heatmap_Bsal_GOI2 + plot_layout(widths = c(0.1, 1))


####Figure 5 ####

Min_annotation_Bd_KOI <- ggplot(KO_Bd_GOI_anot, aes(y = gene, x = 1, fill = Function)) +
  geom_tile() +
  scale_fill_manual(values = subheader_pallete_gray, name = "Function", labels = function(x) str_wrap(x, width = 60)) +
  theme_void() +
  theme(
    legend.position = "none")


Min_annotation_Bsal_KOI <- ggplot(KO_Bsal_GOI_anot, aes(y = gene, x = 1, fill = Function)) +
  geom_tile() +
  scale_fill_manual(values = subheader_pallete_gray, name = "Function", labels = function(x) str_wrap(x, width = 60)) +
  theme_void() +
  theme(
    legend.position = "none")

anotation_legend <- get_legend(annotation_Bd_KOI)

anot_leg_plot<-as_ggplot(anotation_legend) 

Figure_5<- plot_spacer() + anot_leg_plot + plot_spacer() + Min_annotation_Bd_KOI + Heatmap_Bd_GOI2 + Min_annotation_Bsal_KOI + Heatmap_Bsal_GOI2 + plot_layout(widths = c(4, 2, 4, 1, 5, 1, 5))

