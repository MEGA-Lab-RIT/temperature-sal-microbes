#### Supplemental_Figures_7_8_And_Stats####
###Created on: 3/10/26   
###Emma Thompson 
# Creating Supplemental Figure 7: Pathogen abundance on day 18 for each experimental temperature. Log abundance of zoospore genomic equivalents (ZGEs) of Bd (A) and Bsal (B) on day 18 measured by qPCR. Boxplots shows median and interquartile range with whiskers extending to 1.5X interquartile range, with outliers denoted as solid black points and all points as gray points.
# Figure 8: Pathogen abundance over time for each experimental temperature. Mean log zoospore genomic equivalent (ZGEs) and standard error for Bd (A) or Bsal (B) were measured every six days (n = 5) via qPCR. For reference: Bd’s thermal optimum is 20°C, and Bsal’s is 15°C.

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

#### Figure 7 ###

## Bd

Bd_Abun2_1 <- read.csv("R_BdAbundance.csv")

Bd_Abun3 <- Bd_Abun2_1 %>%
  filter(Day == 18)

Bd_Abun3$Day <- factor(Bd_Abun3$Day)


Bd_Abun3$Temperature <- factor(Bd_Abun3$Temperature,
                               levels = c("15", "18", "20", "25", "29"))

Bd_Abun3$logZGE <- as.numeric(as.character(Bd_Abun3$logZGE))

Bd_Abun3_no_zeros<- Bd_Abun3[Bd_Abun3$logZGE != 0, ]

Bd_Abundance <-ggplot( data = Bd_Abun3_no_zeros, aes(x = Temperature, y = logZGE, fill = Temperature)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 1.5, alpha = 0.5, position = position_jitter(0.2)) +
  labs(title = "", x = "", y = "Log (ZGE)") +
  theme_bw() +
  coord_cartesian(ylim = c(0, 16)) +
  theme(
    legend.position = "right",
    legend.key.size = unit(1, 'cm'),
    axis.text.x = element_text()) +
  scale_fill_manual(values = PalleteAdivGTDB, name = "Temperature") 


## Bsal 
Bsal_Abun2_1 <- read.csv("R_BsalAbundance.csv")

Bsal_Abun3 <- Bsal_Abun2_1 %>%
  filter(Day == 18)

Bsal_Abun3$Day <- factor(Bsal_Abun3$Day)


Bsal_Abun3$Temperature <- factor(Bsal_Abun3$Temperature,
                                 levels = c("15", "18", "20", "25", "29"))

Bsal_Abun3$logZGE <- as.numeric(as.character(Bsal_Abun3$logZGE))

Bsal_Abun3_no_zeros<- Bsal_Abun3[Bsal_Abun3$logZGE != 0, ]

Bsal_Abundance <-
  ggplot( data = Bsal_Abun3_no_zeros, aes(x = Temperature, y = logZGE, fill = Temperature)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 1.5, alpha = 0.5, position = position_jitter(0.2)) +
  labs(title= "", x = "", y = "Log (ZGE)") +
  theme_bw() +
  coord_cartesian(ylim = c(0, 16)) +
  theme(legend.position = "right",
    legend.key.size = unit(1, 'cm'),
    axis.text.x = element_text()) +
  scale_fill_manual(values = PalleteAdivGTDB, name = "Temperature") 

## Full Figure 
Abundance<- Bd_Abundance + Bsal_Abundance + plot_layout(guides = "collect") & theme(legend.position = "right")
Abundance__wTags <- Abundance + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = 'bold'))


#### Figure 8 ####

Bd_Abun2_1_no_zeros<- Bd_Abun2_1[Bd_Abun2_1$logZGE != 0, ]

avg_Bd_Abun <- Bd_Abun2_1_no_zeros %>%
  group_by(Temperature, Day) %>%
  summarise(
    mean_ZGE = mean(logZGE, na.rm = TRUE),
    sd_ZGE   = sd(logZGE, na.rm = TRUE),
    n        = sum(!is.na(ZGE)),
    se_ZGE   = sd_ZGE / sqrt(n)) %>%
  ungroup()

Bd_Abun_OverTime <-ggplot(avg_Bd_Abun, aes(x = Day, y = mean_ZGE, color = factor(Temperature))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_ZGE - se_ZGE, ymax = mean_ZGE + se_ZGE), width = 0.5) +
  labs(x = "Day",
       y = "Log(ZGE)",
       color = "Temperature") +
  ylim(0, 18) +
  scale_x_continuous(breaks = c(0, 6, 12, 18)) +
  scale_color_manual(values = c("#f2aa84","#83e291", "#96dcf8", "#d86ecc", "#ffa1e7"))+ 
  theme_minimal()+
  theme(axis.line = element_line(color = "black", linewidth = 0.5))

##Bsal 

Bsal_Abun2_1_no_zeros<- Bsal_Abun2_1[Bsal_Abun2_1$logZGE != 0, ]

avg_Bsal_Abun <- Bsal_Abun2_1_no_zeros %>%
  group_by(Temperature, Day) %>%
  summarise(
    mean_ZGE = mean(logZGE, na.rm = TRUE),
    sd_ZGE   = sd(logZGE, na.rm = TRUE),
    n        = sum(!is.na(ZGE)),
    se_ZGE   = sd_ZGE / sqrt(n)) %>%
  ungroup()

Bsal_Abun_OverTime <-ggplot(avg_Bsal_Abun, aes(x = Day, y = mean_ZGE, color = factor(Temperature))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_ZGE - se_ZGE, ymax = mean_ZGE + se_ZGE), width = 0.5) +
  labs(x = "Day",
       y = "Log(ZGE))",
       color = "Temperature") +
  ylim(0, 18) +
  scale_color_manual(values = c("#f2aa84","#83e291", "#96dcf8", "#d86ecc", "#ffa1e7"))+ 
  theme_minimal() +
  theme(axis.line = element_line(color = "black", linewidth = 0.5))+
  scale_x_continuous(breaks = seq(0, 18, by = 6))

##Full Figure 
Abundance_Over_Time <- (Bd_Abun_OverTime / Bsal_Abun_OverTime) + plot_layout(guides = "collect") &
  theme(legend.position = "right")
Abundance_Over_Time_wTags <- Abundance_Over_Time + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = 'bold'))
print(Abundance_Over_Time_wTags)

####Supplemental Table 8: Two-way ANOVA on pathogen abundance ####

## Bd
Bd_Abun2_1$Temperature <- as.factor(Bd_Abun2_1$Temperature)
Bd_Abun2_1$Day <- as.factor(Bd_Abun2_1$Day)

anova_model_Bd <- aov(ZGE ~ Temperature * Day, data = Bd_Abun2_1)
summary(anova_model_Bd)

## Bsal 
Bsal_Abun2_1$Temperature <- as.factor(Bsal_Abun2_1$Temperature)
Bsal_Abun2_1$Day <- as.factor(Bsal_Abun2_1$Day)

anova_model_Bsal <- aov(ZGE ~ Temperature * Day, data = Bsal_Abun2_1)
summary(anova_model_Bsal)

####Supplemental Table 9: Results of post hoc testing of Bsal abundance by day ####
anova_model_Bsal$Day <- as.factor(anova_model_Bsal$Day)
anova_model_Bsal$Temperature <- as.factor(anova_model_Bsal$Temperature)

Bsal_Day_PostHoc<-TukeyHSD(anova_model_Bsal, "Day")
print(Bsal_Day_PostHoc)

