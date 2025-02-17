library(CRISPRcleanR)
library(tidyverse)

source('src/custom_ccr.R')

data('KY_Library_v1.0')
data('EssGenes.ribosomalProteins')
data('EssGenes.DNA_REPLICATION_cons')
data('EssGenes.KEGG_rna_polymerase')
data('EssGenes.PROTEASOME_cons')
data('EssGenes.SPLICEOSOME_cons')

BAGEL_essential <- read_tsv('data/CEGv2.txt') %>% pull(GENE)
BAGEL_nonEssential <- read_tsv('data/NEGv1.txt') %>% pull(GENE)

## gene-to-sgRNA conversion for reference genes
BAGEL_essential_sgRNAs <- ccr.genes2sgRNAs(KY_Library_v1.0, BAGEL_essential)
BAGEL_nonEssential_sgRNAs <- ccr.genes2sgRNAs(KY_Library_v1.0, BAGEL_nonEssential)

# Gene signatures list
SIGNATURES <- list(Ribosomal_Proteins = EssGenes.ribosomalProteins,
                   DNA_Replication = EssGenes.DNA_REPLICATION_cons,
                   RNA_polymerase = EssGenes.KEGG_rna_polymerase,
                   Proteasome = EssGenes.PROTEASOME_cons,
                   Spliceosome = EssGenes.SPLICEOSOME_cons,
                   CFE = BAGEL_essential,
                   non_essential = BAGEL_nonEssential)

results_path <- 'results'
cell_lines <- list.dirs(results_path, full.names = FALSE, recursive = FALSE)

# Pool all results
sgRNA_level_PrRc_tot <- data.frame(AUC = numeric(), 
                                   Recall = numeric(), 
                                   sigthreshold = numeric(),
                                   cell_line = character())
gene_level_PrRc_tot <- data.frame(AUC = numeric(), 
                                   Recall = numeric(), 
                                   sigthreshold = numeric(),
                                   cell_line = character())
sgRNA_level_ROC_tot <- data.frame(AUC = numeric(), 
                                   Recall = numeric(), 
                                   sigthreshold = numeric(),
                                   cell_line = character())
gene_level_ROC_tot <- data.frame(AUC = numeric(), 
                                   Recall = numeric(), 
                                   sigthreshold = numeric(),
                                   cell_line = character())
gene_signatures_tot <- data.frame(Ribosomal_Proteins = numeric(), 
                                   DNA_Replication = numeric(), 
                                   RNA_polymerase = numeric(),
                                   Proteasome = numeric(),
                                   Spliceosome = numeric(),
                                   CFE = numeric(),
                                   non_essential = numeric(),
                                   cell_line = character())

for (i in cell_lines){
  ## logFC sgRNA corrected
  tmp <- read_tsv(paste0(results_path, '/', i, '/', i, '_corrected_logFCs.tsv'))

  ## Results visualization
  FCs <- tmp$correctedFC
  names(FCs) <- tmp$sgRNA_id
  geneFCs <- ccr.geneMeanFCs(FCs, KY_Library_v1.0)

  ### sgRNA-level
  sgRNA_level_ROC <- ccr.ROC_Curve(FCs, 
                                   BAGEL_essential_sgRNAs,
                                   BAGEL_nonEssential_sgRNAs,
                                   FDRth = 0.05,
                                   expName = label,
                                   display = FALSE)
  sgRNA_level_ROC_tot <- rbind(sgRNA_level_ROC_tot, 
                               data.frame(AUC = sgRNA_level_ROC$AUC,
                                          Recall = sgRNA_level_ROC$Recall,
                                          sigthreshold = sgRNA_level_ROC$sigthreshold,
                                          cell_line = i))
  
  sgRNA_level_PrRc <- ccr.PrRc_Curve(FCs,
                                     BAGEL_essential_sgRNAs,
                                     BAGEL_nonEssential_sgRNAs,
                                     FDRth = 0.05,
                                     expName = label,
                                     display = FALSE)
  sgRNA_level_PrRc_tot <- rbind(sgRNA_level_PrRc_tot, 
                               data.frame(AUC = sgRNA_level_PrRc$AUC,
                                          Recall = sgRNA_level_PrRc$Recall,
                                          sigthreshold = sgRNA_level_PrRc$sigthreshold,
                                          cell_line = i))
  
  ### gene-level
  gene_level_ROC <- ccr.ROC_Curve(geneFCs,
                                  BAGEL_essential,
                                  BAGEL_nonEssential,
                                  FDRth = 0.05,
                                  expName = label,
                                  display = FALSE)
  gene_level_ROC_tot <- rbind(gene_level_ROC_tot, 
                               data.frame(AUC = gene_level_ROC$AUC,
                                          Recall = gene_level_ROC$Recall,
                                          sigthreshold = gene_level_ROC$sigthreshold,
                                          cell_line = i))
  
  gene_level_PrRc <- ccr.PrRc_Curve(geneFCs,
                                    BAGEL_essential,
                                    BAGEL_nonEssential,
                                    FDRth = 0.05,
                                    expName = label,
                                    display = FALSE)
  gene_level_PrRc_tot <- rbind(gene_level_PrRc_tot,
                                    data.frame(AUC = gene_level_PrRc$AUC,
                                                Recall = gene_level_PrRc$Recall,
                                                sigthreshold = gene_level_PrRc$sigthreshold,
                                                cell_line = i))
  
  ### Gene essentiality profile
  Recall_scores <- ccr.custom_VisDepAndSig(FCsprofile = geneFCs,
                                            SIGNATURES = SIGNATURES,
                                            pIs = 6,
                                            nIs = 7,
                                            th = 0.05)
  gene_signatures_tot <- rbind(gene_signatures_tot,
                                    data.frame(Ribosomal_Proteins = Recall_scores[1],
                                                DNA_Replication = Recall_scores[2],
                                                RNA_polymerase = Recall_scores[3],
                                                Proteasome = Recall_scores[4],
                                                Spliceosome = Recall_scores[5],
                                                CFE = Recall_scores[6],
                                                non_essential = Recall_scores[7],
                                                cell_line = i))
}

# Plot summary
## sgRNA-level metrics
sgRNA_level_tot <- rbind(sgRNA_level_ROC_tot %>% mutate(metric = 'ROC'), 
                         sgRNA_level_PrRc_tot %>% mutate(metric = 'PrRc'))

p <- ggplot(sgRNA_level_tot, aes(x = cell_line, y = AUC, fill = metric)) +
  geom_bar(stat='identity', position = 'dodge', color = NA) +
  theme_minimal() +
  labs(x = 'Cell line', y = 'AUROC/AUPrR', fill = 'Metric', title = 'sgRNA-level metrics') +
    theme(plot.title = element_text(hjust = 0.5, size = 16), panel.border = element_blank(), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    aspect.ratio = 1,
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14))
ggsave('results/sgRNA_level_metrics.pdf', p, width = 10, height = 10, dpi = 300)
write_tsv(sgRNA_level_tot, 'results/sgRNA_level_metrics.tsv')

## sgRNA-level recall
p <- ggplot(sgRNA_level_ROC_tot, aes(x = cell_line, y = Recall)) +
  geom_bar(stat='identity', position = 'dodge', color = NA, fill = 'lightblue') +
  theme_minimal() +
  labs(x = 'Cell line', y = 'Recall at 5% FDR', title = 'sgRNA-level recall') +
    theme(plot.title = element_text(hjust = 0.5, size = 16), panel.border = element_blank(), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    aspect.ratio = 1,
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14))
ggsave('results/sgRNA_level_recall.pdf', p, width = 10, height = 10, dpi = 300)
write_tsv(sgRNA_level_ROC_tot, 'results/sgRNA_level_recall.tsv')

## Gene-level metrics
gene_level_tot <- rbind(gene_level_ROC_tot %>% mutate(metric = 'ROC'), 
                         gene_level_PrRc_tot %>% mutate(metric = 'PrRc'))

p <- ggplot(gene_level_tot, aes(x = cell_line, y = AUC, fill = metric)) +
  geom_bar(stat='identity', position = 'dodge', color = NA) +
  theme_minimal() +
  labs(x = 'Cell line', y = 'AUROC/AUPrR', fill = 'Metric', title = 'Gene-level metrics') +
    theme(plot.title = element_text(hjust = 0.5, size = 16), panel.border = element_blank(), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    aspect.ratio = 1,
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14))
ggsave('results/gene_level_metrics.pdf', p, width = 10, height = 10, dpi = 300)
write_tsv(gene_level_tot, 'results/gene_level_metrics.tsv')

## Gene-level recall
p <- ggplot(gene_level_ROC_tot, aes(x = cell_line, y = Recall)) +
  geom_bar(stat='identity', position = 'dodge', color = NA, fill = 'lightblue') +
  theme_minimal() +
  labs(x = 'Cell line', y = 'Recall at 5% FDR', title = 'Gene-level recall') +
    theme(plot.title = element_text(hjust = 0.5, size = 16), panel.border = element_blank(), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    aspect.ratio = 1,
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14))
ggsave('results/gene_level_recall.pdf', p, width = 10, height = 10, dpi = 300)
write_tsv(gene_level_ROC_tot, 'results/gene_level_recall.tsv')

## Gene essentiality profile
### Pivot table wide to long
gene_signatures_tot <- gene_signatures_tot %>%
  pivot_longer(cols = c(Ribosomal_Proteins, DNA_Replication, RNA_polymerase, Proteasome, Spliceosome, CFE, non_essential),
               names_to = 'metric',
               values_to = 'Recall')

p <- ggplot(gene_signatures_tot, aes(x = metric, y = Recall, fill = cell_line)) +
  geom_bar(stat='identity', position = 'dodge', color = NA) +
  theme_minimal() +
  labs(x = 'Gene signature', y = 'Recall at 5% FDR', fill = 'Cell line', title = 'Gene essentiality profile') +
    theme(plot.title = element_text(hjust = 0.5, size = 16), panel.border = element_blank(), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    aspect.ratio = 1,
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)) +
    scale_fill_brewer(palette = 'Paired')
ggsave('results/gene_signatures.pdf', p, width = 10, height = 10, dpi = 300)
write_tsv(gene_signatures_tot, 'results/gene_signatures.tsv')
