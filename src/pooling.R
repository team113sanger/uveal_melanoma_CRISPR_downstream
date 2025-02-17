library(tidyverse)

results_path <- 'results'
cell_lines <- list.dirs(results_path, full.names = FALSE, recursive = FALSE)

# Pool all samples
logFC_sgRNA_corrected <- NULL
logFC_gene_corrected <- NULL
BAGEL_unscaled <- NULL
BAGEL_scaled <- NULL
MAGeCK_gene_uncorrected_beta <- NULL
MAGeCK_gene_corrected_beta <- NULL
MAGeCK_gene_uncorrected_fdr <- NULL
MAGeCK_gene_corrected_fdr <- NULL

for (i in cell_lines){
    ## logFC sgRNA corrected
    tmp <- read_tsv(paste0(results_path, '/', i, '/', i, '_corrected_logFCs.tsv')) %>% 
        select(sgRNA_id, genes, correctedFC)
    tmp_genes <- tmp %>% 
        group_by(genes) %>% 
        mutate(correctedFC = mean(correctedFC, na.rm = TRUE)) %>%
        ungroup() %>%
        select(genes, correctedFC) %>%
        distinct()
    colnames(tmp)[ncol(tmp)] <- i
    colnames(tmp_genes)[ncol(tmp_genes)] <- i
    
    if (is.null(logFC_sgRNA_corrected)){
        logFC_sgRNA_corrected <- tmp
    } else {
        logFC_sgRNA_corrected <- logFC_sgRNA_corrected %>% full_join(tmp, by = c('sgRNA_id', 'genes'))
    }

    ## logFC gene corrected    
    if (is.null(logFC_gene_corrected)){
        logFC_gene_corrected <- tmp_genes
    } else {
        logFC_gene_corrected <- logFC_gene_corrected %>% full_join(tmp_genes, by = c('genes'))
    }

    ## BAGEL unscaled
    tmp <- read_tsv(paste0(results_path, '/', i, '/', i, '_bf.tsv'))
    colnames(tmp)[ncol(tmp)] <- i

    if (is.null(BAGEL_unscaled)){
        BAGEL_unscaled <- tmp
    } else {
        BAGEL_unscaled <- BAGEL_unscaled %>% full_join(tmp, by = 'GENE')
    }

    ## BAGEL scaled
    tmp <- read_tsv(paste0(results_path, '/', i, '/', i, '_bf_scaled.tsv'))
    colnames(tmp)[ncol(tmp)] <- i

    if (is.null(BAGEL_scaled)){
        BAGEL_scaled <- tmp
    } else {
        BAGEL_scaled <- BAGEL_scaled %>% full_join(tmp, by = 'GENE')
    }

    ## MAGeCK gene uncorrected beta
    item <- paste0(i, '|beta')
    tmp <- read_tsv(paste0(results_path, '/', i, '/mageck_uncorrected/', i, '.gene_summary.txt')) %>% 
        select(Gene, item)
    colnames(tmp)[ncol(tmp)] <- i

    if (is.null(MAGeCK_gene_uncorrected_beta)){
        MAGeCK_gene_uncorrected_beta <- tmp
    } else {
        MAGeCK_gene_uncorrected_beta <- MAGeCK_gene_uncorrected_beta %>% full_join(tmp, by = 'Gene')
    }

    ## MAGeCK gene uncorrected fdr
    item <- paste0(i, '|fdr')
    tmp <- read_tsv(paste0(results_path, '/', i, '/mageck_uncorrected/', i, '.gene_summary.txt')) %>% 
        select(Gene, item)
    colnames(tmp)[ncol(tmp)] <- i

    if (is.null(MAGeCK_gene_uncorrected_fdr)){
        MAGeCK_gene_uncorrected_fdr <- tmp
    } else {
        MAGeCK_gene_uncorrected_fdr <- MAGeCK_gene_uncorrected_fdr %>% full_join(tmp, by = 'Gene')
    }

    ## MAGeCK gene corrected beta
    item <- paste0(i, '|beta')
    tmp <- read_tsv(paste0(results_path, '/', i, '/mageck_corrected/', i, '.gene_summary.txt')) %>% 
        select(Gene, item)
    colnames(tmp)[ncol(tmp)] <- i

    if (is.null(MAGeCK_gene_corrected_beta)){
        MAGeCK_gene_corrected_beta <- tmp
    } else {
        MAGeCK_gene_corrected_beta <- MAGeCK_gene_corrected_beta %>% full_join(tmp, by = 'Gene')
    }

    ## MAGeCK gene corrected fdr
    item <- paste0(i, '|fdr')
    tmp <- read_tsv(paste0(results_path, '/', i, '/mageck_corrected/', i, '.gene_summary.txt')) %>% 
        select(Gene, item)
    colnames(tmp)[ncol(tmp)] <- i

    if (is.null(MAGeCK_gene_corrected_fdr)){
        MAGeCK_gene_corrected_fdr <- tmp
    } else {
        MAGeCK_gene_corrected_fdr <- MAGeCK_gene_corrected_fdr %>% full_join(tmp, by = 'Gene')
    }
}

# save pooled data
write_tsv(logFC_sgRNA_corrected, file = paste0(results_path, '/logFC_sgRNA_corrected.tsv'))
write_tsv(logFC_gene_corrected, file = paste0(results_path, '/logFC_gene_corrected.tsv'))
write_tsv(BAGEL_unscaled, file = paste0(results_path, '/BAGEL_unscaled.tsv'))
write_tsv(BAGEL_scaled, file = paste0(results_path, '/BAGEL_scaled.tsv'))
write_tsv(MAGeCK_gene_uncorrected_beta, file = paste0(results_path, '/MAGeCK_gene_uncorrected_beta.tsv'))
write_tsv(MAGeCK_gene_corrected_beta, file = paste0(results_path, '/MAGeCK_gene_corrected_beta.tsv'))
write_tsv(MAGeCK_gene_uncorrected_fdr, file = paste0(results_path, '/MAGeCK_gene_uncorrected_fdr.tsv'))
write_tsv(MAGeCK_gene_corrected_fdr, file = paste0(results_path, '/MAGeCK_gene_corrected_fdr.tsv'))
