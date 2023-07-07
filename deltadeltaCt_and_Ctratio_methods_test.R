# DeltaDeltaCt method computation #
# The script computes the deltadeltaCt method from qPCR Ct values and generates #
# individual plots for each gene.
# Input : dataframe containing list of target genes in columns and samples in rows #
# Víctor A. López-Agudelo 
# 06.07.2023

rm(list = ls())

library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(glue)

#setwd("~/Library/CloudStorage/OneDrive-InstitutfürKlinischeMolekularbiologie,Christian-Albrechts-UniversitätKiel/MedCluster/other_stuffs/qPCR_computation")



# Example Data to test the function
# Data taken from : https://toptipbio.com/delta-delta-ct-pcr/
# data <- data.frame(Treatment = c("Control","Control","Control","Treated1","Treated2","Treated3"),
#                   TargetGene = c(30.55,30.55,30.79,25.83,25.63,25.69),
#                   ReferenceGene = c(17.18,16.96,17.11,18.0,17.95,17.88))


# targetGene = "TargetGene"
# referenceGene = "ReferenceGene"
# controlGroup = "Control"
# treatmentColumn = "Treatment"

# Delta Delta Ct method in R #

computeDeltaDeltaCt <- function(data, targetGene, referenceGene, controlGroup, treatmentColumn) {
  
  # Subset the dataframe for target and reference genes
  targetCt <- data[, targetGene]
  referenceCt <- data[, referenceGene]
  
  # Calculate Ct ratio
  deltaCt <- targetCt - referenceCt
  Ctratio <- 2^(-targetCt)/2^ (-referenceCt)
  
  # Calculate average delta Ct for control and experimental groups
  controlDeltaCt <- mean(deltaCt[data[, treatmentColumn] == controlGroup])
  #experimentalDeltaCt <- ave(deltaCt, data[, treatmentColumn], FUN = mean)
  
  # Calculate delta delta Ct
  #deltaDeltaCt <- experimentalDeltaCt - controlDeltaCt
  deltaDeltaCt <- deltaCt - controlDeltaCt
  
  # Calculate fold change
  foldChange <- 2^(-deltaDeltaCt)
  
  # Return the results as a dataframe
  result <- data.frame(Treatment = data[, treatmentColumn],
                       DeltaDeltaCt = deltaDeltaCt,
                       FoldChange = foldChange,
                       Ctratio = Ctratio,
                       TargetGene = targetGene)
  
  return(result)
}

##::::::::::::::::::::::: Expression test data :::::::::::::::::::::::::::::::##

Ct_table <- read.table("stim_3t3_ct_gene_ct_housekeeping.txt", sep = "\t", header=T)

# Comparison WT vs KO for each "Treatment" #

targetgene_vector <- colnames(Ct_table)[-c(1:2)]
genotype_vector <- c("wt","ko")

delta_delta_ct_table <- NULL

for (z in 1:length(targetgene_vector)){

  for (x in 1:length(genotype_vector)){
  
  delta_delta_ct_output <- computeDeltaDeltaCt(Ct_table %>% filter(genotype == glue("{genotype_vector[x]}")),targetgene_vector[z],"actb","ctrl","sample") %>%
    mutate(Genotype = genotype_vector[x])
  
  delta_delta_ct_table <- rbind(delta_delta_ct_table,delta_delta_ct_output)
  
  
  }
  
}

delta_delta_ct_table <- delta_delta_ct_table %>%
  mutate(TreatmentGenotype = paste0(Treatment,"-",Genotype)) %>%
  mutate(sqrt_ct = sqrt(Ctratio))

# Order of Treatment and Genotype #
delta_delta_ct_table$TreatmentGenotype <- factor(delta_delta_ct_table$TreatmentGenotype ,
                                        levels = c("ctrl-wt","ctrl-ko","blm-wt","blm-ko","tgf_egf-wt","tgf_egf-ko","tgf-wt","tgf-ko"), 
                                        ordered = T)


# Create Individual Barplot Per Gene - FoldChange Approach #
dir.create("plots")
dir.create("plots/FoldChange_based_on_Ctrl_WT_or_Ctrl_KO")

for (j in 1:length(targetgene_vector)) {
  
ggbarplot(delta_delta_ct_table %>% filter(TargetGene == targetgene_vector[j]) , x = "TreatmentGenotype", y = "FoldChange", color = "Genotype",
          add = c("mean_se", "dotplot"),
          palette= c("wt"="darkgrey","ko"="darkred"),
          facet.by = c("TargetGene")
) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y = "FoldChange", x = "Treatment/Genotype") + stat_compare_means(method = "wilcox.test", label = "p.label",
                                                                                                                    comparisons = list(c("ctrl-wt","ctrl-ko"),
                                                                                                                                       c("blm-wt","blm-ko"),
                                                                                                                                       c("tgf_egf-wt","tgf_egf-ko"),
                                                                                                                                       c("tgf-wt","tgf-ko")))
ggsave(filename = glue("plots/FoldChange_based_on_Ctrl_WT_or_Ctrl_KO/barplot_FoldChange_Normalized-Ctrls-WTvsKO_{targetgene_vector[j]}.pdf"), width = 5,height=5)    

}


# Create Individual Barplot Per Gene - Relative Expression Approach #

dir.create("plots/RelativeExpressions")

for (j in 1:length(targetgene_vector)) {
  
  ggbarplot(delta_delta_ct_table %>% filter(TargetGene == targetgene_vector[j]) , x = "TreatmentGenotype", y = "Ctratio", color = "Genotype",
            add = c("mean_se", "dotplot"),
            palette= c("wt"="darkgrey","ko"="darkred"),
            facet.by = c("TargetGene")
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y = "Relative Expression", x = "Treatment/Genotype") + stat_compare_means(method = "wilcox.test", label = "p.label",
                                                                                                                                       comparisons = list(c("ctrl-wt","ctrl-ko"),
                                                                                                                                                          c("blm-wt","blm-ko"),
                                                                                                                                                          c("tgf_egf-wt","tgf_egf-ko"),
                                                                                                                                                          c("tgf-wt","tgf-ko")))
  ggsave(filename = glue("plots/RelativeExpressions/barplot_RelativeExpression-WTvsKO_{targetgene_vector[j]}.pdf"), width = 5,height=5)    
  
}
