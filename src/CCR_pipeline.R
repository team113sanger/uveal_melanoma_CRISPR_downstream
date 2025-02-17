library(CRISPRcleanR)
library(tidyverse)

source('src/custom_ccr.R')

data('KY_Library_v1.0')
data('EssGenes.ribosomalProteins')
data('EssGenes.DNA_REPLICATION_cons')
data('EssGenes.KEGG_rna_polymerase')
data('EssGenes.PROTEASOME_cons')
data('EssGenes.SPLICEOSOME_cons')

# Load data
crispr_screens <- read_tsv('data/01_raw_count_matrix.tsv')
crispr_screens <- crispr_screens %>% select(Plasmid, everything())

BAGEL_essential <- read_tsv('data/CEGv2.txt') %>% pull(GENE)
BAGEL_nonEssential <- read_tsv('data/NEGv1.txt') %>% pull(GENE)

## gene-tosgRNA conversion for reference genes
BAGEL_essential_sgRNAs <- ccr.genes2sgRNAs(KY_Library_v1.0, BAGEL_essential)
BAGEL_nonEssential_sgRNAs <- ccr.genes2sgRNAs(KY_Library_v1.0, BAGEL_nonEssential)

## Gene signatures list
SIGNATURES <- list(Ribosomal_Proteins = EssGenes.ribosomalProteins,
                   DNA_Replication = EssGenes.DNA_REPLICATION_cons,
                   RNA_polymerase = EssGenes.KEGG_rna_polymerase,
                   Proteasome = EssGenes.PROTEASOME_cons,
                   Spliceosome = EssGenes.SPLICEOSOME_cons,
                   CFE = BAGEL_essential,
                   non_essential = BAGEL_nonEssential)

# CCR pipeline
ncontrols <- 1

cell_lines <- colnames(crispr_screens)[7:ncol(crispr_screens)]
cell_lines <- split(cell_lines, f = rep(seq_len(10), each = 3))

for (i in seq_along(cell_lines)){
  replicates <- cell_lines[[i]]
  print(paste0('Analyzing replicates:  ', paste0(replicates, collapse = '  ')))
  
  
  ## Subset dataset
  fn <- cbind(crispr_screens[, c(3,5,1)], crispr_screens[, replicates])
  fn <- as.data.frame(fn)
  colnames(fn)[1:2] <- c('sgRNA', 'gene')
  
  
  ## Sample label
  label <- replicates[1]
  label <- str_split(label, '_')[[1]][1]
  
  outdir <- paste0('results/', label, '/')
  
  
  ## Create directory if it doesn't exist
  if (!file.exists(outdir)){
    dir.create(outdir)
  }
  

  ## Save raw counts
  write_tsv(fn, file = paste0(outdir, label, '_raw_counts.tsv'))

  
  ## Normalize data
  print(paste0('Normalizing data for ', label))
  normANDfcs <- ccr.NormfoldChanges(Dframe = fn,
                                    min_reads = 30,
                                    EXPname = label,
                                    ncontrols = ncontrols,
                                    libraryAnnotation = KY_Library_v1.0,
                                    display = FALSE,
                                    saveToFig = TRUE,
                                    outdir = outdir)
  write_tsv(normANDfcs$norm_counts, file = paste0(outdir, label, '_count_norm.tsv'))
  write_tsv(normANDfcs$logFCs, file = paste0(outdir, label, '_logFCs.tsv'))
  
  
  ## Genome sorting and correction
  print(paste0('Genome sorting and correction for ', label))
  gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs$logFCs, KY_Library_v1.0)
  correctedFCs <- ccr.GWclean(gwSortedFCs, 
                              label = 'CCR_correction',
                              return.segments.adj = TRUE,
                              display = FALSE,
                              saveTO = outdir,
                              verbose = -1)
  
  
  ## BAGEL
  ### Prepare corrected logFCs for BAGEL
  print(paste0('Run BAGEL on ', label))
  bf_correctedFCs <- correctedFCs$corrected_logFCs %>% 
    rownames_to_column(var = 'sgRNA_id') %>%
    select(sgRNA_id, genes, correctedFC)
  
  write_tsv(bf_correctedFCs, file = paste0(outdir, label, '_input_bagel_corrected_logFCs.tsv'))
  write_tsv(correctedFCs$corrected_logFCs %>% rownames_to_column(var = 'sgRNA_id'), file = paste0(outdir, label, '_corrected_logFCs.tsv'))
  write_tsv(correctedFCs$segments, file = paste0(outdir, label, '_segments.tsv'))
  write_tsv(correctedFCs$segments_adj, file = paste0(outdir, label, '_corrected_segments.tsv'))
  
  ### run BAGEL
  system(paste0('bash src/run_BAGEL.sh ', label, '/', label, '_input_bagel_corrected_logFCs.tsv', ' ', label, '/', label, '_bf.tsv'), wait = TRUE)
  
  ### scale BFs using 5% FDR
  BFs <- read_tsv(paste0('results/', label, '/', label, '_bf.tsv'))
  BFs <- BFs %>% pull(BF, name = GENE)
  
  sigthreshold <- -ccr.ROC_Curve(-BFs,
                                  BAGEL_essential,
                                  BAGEL_nonEssential,
                                  FDRth = 0.05,
                                  display = FALSE)$sigthreshold
  
  BFs_scaled <- data.frame(GENE = names(BFs), BF = BFs-sigthreshold)
  write_tsv(BFs_scaled, file = paste0(outdir, label, '_bf_scaled.tsv'))
  
  
  ## Corrected counts
  print(paste0('Compute corrected counts for ', label))
  correctedCounts <- ccr.correctCounts(label, 
                                       normANDfcs$norm_counts,
                                       correctedFCs,
                                       KY_Library_v1.0,
                                       minTargetedGenes = 3,
                                       ncontrols = ncontrols,
                                       OutDir = outdir,
                                       verbose = -1)
  write_tsv(correctedCounts, file = paste0(outdir, label, '_corrected_counts.tsv'))
  
  
  ## MAGeCK
  ### Prepare corrected logFCs for MAGeCK
  uncorrected_fn <- normANDfcs$norm_counts
  Cnames <- colnames(uncorrected_fn)[3:(2 + ncontrols)]
  Tnames <- colnames(uncorrected_fn)[(3 + ncontrols):ncol(uncorrected_fn)]
  
  ### Implement design matrix
  assign(label, c(0, 1, 1, 1))
  
  design_matrix <- data.frame(Samples = c(Cnames, Tnames), baseline = 1, label = get(label))
  colnames(design_matrix)[3] <- label
  
  write_tsv(design_matrix, file = paste0(outdir, 'design_matrix.tsv'))
  
  ### run MAGeCK on uncorrected counts
  print(paste0('Run MAGeCK on uncorrected counts for ', label))
  uncorrected_fn_path <- paste0(outdir, label, '_count_norm.tsv')
  design_matrix_path <- paste0(outdir, 'design_matrix.tsv')
  system(paste0('bash src/run_MAGeCK.sh ', uncorrected_fn_path, ' ', design_matrix_path, ' ', label, ' ', 'mageck_uncorrected'), wait = TRUE)
  
  write_tsv(design_matrix, file = paste0(outdir, 'design_matrix.tsv'))
  
  ### run MAGeCK on corrected counts
  print(paste0('Run MAGeCK on corrected counts for ', label))
  corrected_fn_path <- paste0(outdir, label, '_corrected_counts.tsv')
  design_matrix_path <- paste0(outdir, 'design_matrix.tsv')
  system(paste0('bash src/run_MAGeCK.sh ', uncorrected_fn_path, ' ', design_matrix_path,' ', label, ' ', 'mageck_corrected'), wait = TRUE)
  
  ### Assess the impact on phenotype
  uncorrected_gs_fn <- paste0('results/', label, '/mageck_uncorrected/', label, '.gene_summary.txt')
  corrected_gs_fn <- paste0('results/', label, '/mageck_corrected/', label, '.gene_summary.txt')
  
  pdf(paste0('results/', label, '/Impact_on_phenotype.pdf'), width = 10, height = 10)
  impacted_phenotype_fdr <- ccr.custom_impactOnPhenotype(MO_uncorrectedFile = uncorrected_gs_fn,
                                                         MO_correctedFile = corrected_gs_fn,
                                                         expName = label,
                                                         sigFDR = 0.05,
                                                         display = TRUE)
  dev.off()
  
  
  ## Results visualization
  FCs <- correctedFCs$corrected_logFCs$correctedFC
  names(FCs) <- rownames(correctedFCs$corrected_logFCs)
  geneFCs <- ccr.geneMeanFCs(FCs, KY_Library_v1.0)
  
  ### sgRNA-level
  pdf(paste0('results/', label, '/sgRNA_ROC.pdf'))
  sgRNA_level_ROC <- ccr.ROC_Curve(FCs, 
                                   BAGEL_essential_sgRNAs,
                                   BAGEL_nonEssential_sgRNAs,
                                   FDRth = 0.05,
                                   expName = label,
                                   display = TRUE)
  dev.off()
  
  pdf(paste0('results/', label, '/sgRNA_PrRc.pdf'))
  sgRNA_level_PrRc <- ccr.PrRc_Curve(FCs,
                                     BAGEL_essential_sgRNAs,
                                     BAGEL_nonEssential_sgRNAs,
                                     FDRth = 0.05,
                                     expName = label,
                                     display = TRUE)
  dev.off()
  
  ### gene-level
  pdf(paste0('results/', label, '/gene_ROC.pdf'))
  gene_level_ROC <- ccr.ROC_Curve(geneFCs,
                                  BAGEL_essential,
                                  BAGEL_nonEssential,
                                  FDRth = 0.05,
                                  expName = label,
                                  display = TRUE)
  dev.off()
  
  pdf(paste0('results/', label, '/gene_PrRc.pdf'))
  gene_level_PrRc <- ccr.PrRc_Curve(geneFCs,
                                    BAGEL_essential,
                                    BAGEL_nonEssential,
                                    FDRth = 0.05,
                                    expName = label,
                                    display = TRUE)
  dev.off()
  
  ### Gene essentiality profile
  pdf(paste0('results/', label, '/gene_signatures.pdf'))
  Recall_scores <- ccr.VisDepAndSig(FCsprofile = geneFCs,
                                    SIGNATURES = SIGNATURES,
                                    TITLE = label,
                                    pIs = 6,
                                    nIs = 7,
                                    th = 0.05,
                                    plotFCprofile = TRUE)
  dev.off()
}
