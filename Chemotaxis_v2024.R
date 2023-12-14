# Investigating Chemotaxis gene expression and Cell infiltration into CL lesions - a focus on T cell signatures.



# Introduction ----
# In this analysis of deep exploration, I investigate the signatures of chemotaxis  in lesions from CL patients, and found
# that CCR5 is associated with delayed healing in patients.

# Libraries ----
library(ggthemes)
library(Hmisc)
library(gplots)
#library(naniar)
library(corrplot)
library(broom)
library(ggrepel)
#library(viridis)
#library(vegan)
#library(fcros)
library(ggraph)
library(corrplot)
library(corrr)
library(ggpubr)
library(tidygraph)
#library(reshape2)
#library(gt)
#library(patchwork)

library(tidyverse)

#'%ni%' <- Negate('%in%')

# Import files from 2nd Lesion dataset (Amorim et al. 2019) ----
# The files imported here were processed from the gene expression matrix GSE127831_Amorim_GEO_raw.txt available at
# GEO #GSE127831.

# targets:
targets <- read_csv("~/Dropbox/CamilaAnalysis/Lesion2ndRNAseqDataset2019/GSE127831_studydesign_mod.csv")

#targets_CL <- targets %>%
#  filter(group_RNAseq == "CL")

#targets <- targets %>% # 2nd
#  select(sample_RNAseq,
#         group_RNAseq,
#         treatment_outcome)

#samples_sbv <- targets %>%
  #filter(treatment_outcome %in% c("cure","failure"))
#  filter(treatment_outcome %in% c("Cure","Failure"))
#samples_sbv <- samples_sbv$sample_RNAseq

samples_cure <- targets %>%
  #filter(treatment_outcome == "cure")
  filter(treatment_outcome == "Cure")
samples_cure <- samples_cure$sample_RNAseq

samples_failure <- targets %>%
  # filter(treatment_outcome == "failure")
  filter(treatment_outcome == "Failure")
samples_failure <- samples_failure$sample_RNAseq

samples_CL <- targets %>%
  filter(group_RNAseq == "CL")
samples_CL <- samples_CL$sample_RNAseq

samples_HS <- targets %>%
  filter(group_RNAseq == "HS")
samples_HS <- samples_HS$sample_RNAseq

# myTopHits:
load("../Lesion2ndRNAseqDataset2019/Robjects/myTopHits")

# notfilt_norm_log2cpm:
load("../Lesion2ndRNAseqDataset2019/Robjects/notfilt_norm_log2cpm")
load("../Lesion2ndRNAseqDataset2019/Robjects/log2.cpm.filtered.norm")

# Cell estimations - res_cell:
# Cell estimations were performed with MCP-counter and xCell, using immunedeconv package

load("../Lesion2ndRNAseqDataset2019/Robjects/res_xcell")
load("../Lesion2ndRNAseqDataset2019/Robjects/res_mcp_counter")

## model mcp counter a bit:
res_mcp_counter
res_mcp_counter <- res_mcp_counter %>%
  filter(cell_type != "Cancer associated fibroblast") # these scores are not important here
res_mcp_counter <- as.matrix(column_to_rownames(res_mcp_counter, "cell_type"))
rownames(res_mcp_counter)[3] <- "CTL score"

## and then model xCell:
res_xcell <- res_xcell %>%
  filter(cell_type != "microenvironment score", # these scores are not important here
         cell_type != "stroma score",
         cell_type != "immune score",
         cell_type != "Cancer associated fibroblast",
         cell_type != "Granulocyte-monocyte progenitor") # these ones were not detected in the dataset, many zeros
res_xcell <- as.matrix(column_to_rownames(res_xcell, "cell_type"))

# Cell estimation method chosen:
res_cell <- res_mcp_counter
#res_cell <- res_xcell

# Import chemokine description file and visualization (This chunk was not used) ----
#description_pre <- read_csv("outputs_inputs//description.csv")
#description <- description_pre %>%
#  distinct(Gene, .keep_all = TRUE) %>%
#  drop_na(Gene) %>%
#  filter(Gene != "CCL4L1")

##CCL4L1 --> Annotation category: only annotated on alternate loci in reference assembly
##Note: In NCBI's Homo sapiens Annotation Release 106, this GeneID (GeneID:388372) was annotated as the CCL4L2 gene copy of this chemokine (C-C motif) ligand 4-like family member, while highly similar GeneID:9560 was annotated as the CCL4L1 gene copy. However, it has now been determined that these gene copies can be distinguished based on 3' terminal exon splice site usage, as described in PubMedID:15843566. Genomic sequence analysis at the annotated locations indicates that GeneID:388372 actually represents the CCL4L1 gene copy, while GeneID:9560 represents the CCL4L2 gene copy. The gene symbols have therefore been switched accordingly. The CCL4L1 gene copy is present only on the ALT_REF_LOCI_2 alternate haplotype of the GRCh38 human reference genome assembly, while the CCL4L2 gene copy is present on the primary chromosome 17 assembly (NC_000017.11) and on both the ALT_REF_LOCI_1 and ALT_REF_LOCI_2 alternate haplotypes. [17 Jun 2014]. https://www.ncbi.nlm.nih.gov/gene/388372
##CXCL4 --> PF4
##CXCL4VI --> variant of PF4: PF4V1
##CXCL7 --> PPBP
##CCR11 --> ACKR4
##CXCR7 --> ACKR3

#description %>% 
#  select(1,3,4,5) %>%
#  gt()

#write_tsv(as.data.frame(description %>% 
#                          select(1,3,4,5)), "outputs_inputs/Description File clean.txt")

# List of human chemokines and receptors:
#chemokines_list <- description %>% filter(Family == "chemokine")
#chemokines_list <- chemokines_list$Gene
#receptor_list <- description %>% filter(Family == "receptor")
#receptor_list <- receptor_list$Gene

#write_tsv(as.data.frame(chemokines_list), "outputs_inputs/chemokine_list.txt")
#write_tsv(as.data.frame(receptor_list), "outputs_inputs/receptor_list.txt")

# List of chemokines and receptors selected by Lais Sacramento in Feb, 2023 ----
chemokines_list <- c("CCL3","CCL4",
                     "CXCL9","CXCL10","CXCL11",
                     "CXCL16",
                     "CXCL13",
                     "CCL17","CCL21",
                     "CCL1","CCL8","CCL16","CCL18")
receptor_list <- c("CCR5","CXCR3","CXCR6","CXCR5","CCR4","CCR8")

  
# DGE analysis highlights chemokines and receptors ----
# Volcano:
myTopHits_1 <- myTopHits %>%
  mutate(siggenes1 = case_when(
    logFC > 0.59 & adj.P.Val < 0.01 ~ "over",
    logFC < -0.59 & adj.P.Val < 0.01 ~ "under",
    TRUE ~ "not_sig")) %>%
  mutate(siggenes2 = case_when(
    logFC > 0.59 & adj.P.Val < 0.01 | logFC < -0.59 & adj.P.Val < 0.01 ~ geneSymbol,
    TRUE ~ NA_character_)) %>%
  mutate(siggenes3 = case_when(
    geneSymbol %in% c("CCL3","CCL4") ~ "c_biomarker",
    geneSymbol %in% chemokines_list[3:13] ~ "a_chemokine",
    geneSymbol %in% receptor_list ~ "b_receptor",
    TRUE ~ "other")) %>%
  mutate(siggenes4 = case_when(
    geneSymbol %in% c(chemokines_list, receptor_list) ~ geneSymbol,
    TRUE ~ NA_character_))

write_tsv(myTopHits_1, "outputs_inputs/myTopHits_1_DGE.txt")
# ggplot:
myTopHits_1 %>%
  arrange(desc(siggenes3)) %>%
  ggplot(., aes(y=-log10(adj.P.Val), x=logFC,
                color = siggenes3, size = siggenes3, alpha = siggenes3)) +
  geom_hline(yintercept = -log10(0.01), color="dark gray") + #geom_vline(xintercept = 0.59) + 
  geom_jitter() +
  theme_stata() +
  scale_color_manual(values = c(Salvador[4],Salvador[1],Salvador[5],"dark gray")) +
  scale_alpha_manual(values = c(1,1,1,0.4)) +
  scale_size_manual(values = c(2.5,2.5,2.5,1)) +
  geom_text_repel(aes(label = siggenes4), size = 5, fontface=4, force = 90) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 17, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 17, angle = 0), axis.title = element_text(size = 17),
        legend.position="none", legend.title = element_text(size = 17)) +
  #ylim(0,3.35) +
  xlab(paste("logFC Lesion vs. HS"))
