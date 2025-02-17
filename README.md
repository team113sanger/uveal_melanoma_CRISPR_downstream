## CRISPR - Single guide RNA downstream analysis
This report describes the methodology used in the analysis of 10 genome-wide knockout CRISPR-Cas9 screens for *Targeting the CDS1/2 axis as a therapeutic strategy in uveal melanoma and pan-cancer*. All the findings and results can be found in the `results` subfolder. 

IMPORTANT: We used the single-guide RNA (sgRNA) guides in common with Human CRISPR Library v1.0 only, as version 1.1 of the same library contains additional guides targeting essential genes with 10 guides per gene. This would result in an unbalanced targeting of these essential genes, therefore we decided to remove them. Nevertheless, the Human CRISPR Library v1.1 has superior chemistry.

### CRISPRcleanR
Single-guide RNA count pre-processing and bias correction
We performed all computational analyses considering only the sgRNAs in common with the Human CRISPR Library v1.0. 

The initial matrix contained the raw counts of 10 cell lines screened in 3 technical replicates and 1 plasmid. For each cell line, considering the technical replicates and the plasmid, we applied the following pipeline to preprocess and correct the sgRNA counts using CRISPRcleanR v3.0.1:
-	sgRNAs with less than 30 read counts in the plasmid were first removed
-	sgRNA raw counts were normalised by their total number in the replicate
-	Log fold changes (logFCs) for individual sgRNAs were quantified between post-library-transduction read counts and library plasmid read counts at the individual replicate level.
-	We corrected the gene-independent response to CRISPR-Cas9 targeting using the default parameters.

In addition, to assess the impact on phenotype by CRISPRcleanR correction, we computed corrected normalised sgRNAs’ treatment counts back from the corrected logFCs by applying the inverse transformation described in (Iorio et al. 2018). The data in this format can be used as input by MAGeCK.


### Quality control assessments
We performed several checks to control the quality of the correction in each sample:
-	We tested the genome-wide profile of sgRNAs’ logFCs (or gene-level logFCs, averaged across same-gene targeting sgRNAs) as a classifier of reference sets of core-fitness essential (CFE) and non-essential genes (NEG). These two reference sets were downloaded from the BAGEL GitHub repository (https://github.com/hart-lab/bagel). We computed the Area Under the ROC (AUROC) and the Area under the Precision-recall Curve (AUPrC).
-	We inspected enrichments of predefined sets of core-fitness essential genes near the top of the genome-wide essentiality profiles (composed of sgRNA or gene-level logFCs ranked in increasing order) and computed their classification recall at a 5% FDR determined as detailed in (Dempster et al. 2019). In addition to the CFE genes, we considered sets of genes involved in housekeeping cellular processes (i.e. ribosomal proteins, DNA replication, RNA polymerase, proteosome, spliceosome) assembled from the Molecular Signature Database (MSigDB).
-	To evaluate the effect of the CRISPRcleanR correction on the genes showing a significant loss/gain-of-fitness effect (fitness genes) in the uncorrected data, we performed a comparison of fitness gene sets (computed with MAGeCK) before (normalised uncorrected sgRNAs’ treatment counts) and after (normalised corrected sgRNAs’ treatment counts) CRISPRcleanR correction. In both cases, genes with a significant loss-of-fitness (LoF) effect were defined as having a beta score < 0 and an FDR < 5%, while genes with a significant gain-of-fitness (GoF) effect were defined as having a beta score > 0 and an FDR < 5%.



### Execution of BAGEL
We run BAGEL2 (https://github.com/hart-lab/bagel) on the corrected sgRNAs logFC profile for each cell line. We used CFEv2 and NEGv1 as reference genes (available in the same repository) and used the default parameters (10-fold cross-validation).
The resulting gene-level Bayes Factors (BF) were scaled at 5% FDR using the same procedure mentioned in the previous section.




### Execution of MAGeCK MLE
We ran MAGeCK MLE v0.5.9.5 (https://sourceforge.net/projects/mageck/) on the normalised uncorrected sgRNAs’ treatment counts and on the normalised corrected sgRNAs’ treatment counts (after applying the inverse transformation on the corrected logFCs as mentioned in the previous section) using the default parameters.



## Contents
Results repository structure (label is a boilerplate for the cell line name):

-	label = repository containing all the results for a cell line.
-	BAGEL_unscaled.tsv = pooled matrix of gene-level BFs (each cell line is a column).
-	BAGEL_scaled.tsv = pooled matrix of gene-level BFs (each cell line is a column). The BFs are scaled at 5% FDR using the reference essential and non-essential genes.
-	gene_level_metrics.pdf = AUROC/AUPrR curves at the gene level.
-	gene_level_recall.pdf = recall at 5% FDR at the gene level.
-	gene_signatures.pdf = recall at 5% FDR of gene signatures from MSigDB at the gene level.
-	logFC_sgRNA_corrected.tsv = pooled matrix of corrected sgRNAs’ logFCs (each cell line is a column).
-	logFC_gene_corrected.tsv = pooled matrix of corrected sgRNAs’ logFCs (each cell line is a column).
-	MAGeCK_gene_uncorrected_beta.tsv = pooled matrix of MAGeCK beta scores from normalised uncorrected sgRNAs’ treatment counts (each cell line is a column).
-	MAGeCK_gene_corrected_beta.tsv = pooled matrix of MAGeCK beta scores from  normalised corrected sgRNAs’ treatment counts (each cell line is a column).
-	MAGeCK_gene_uncorrected_fdr.tsv = pooled matrix of MAGeCK FDR scores from  normalised uncorrected sgRNAs’ treatment counts (each cell line is a column).
-	MAGeCK_gene_corrected_fdr.tsv = pooled matrix of MAGeCK FDR scores from  normalised corrected sgRNAs’ treatment counts (each cell line is a column).
-	sgRNA_level_metrics.pdf = AUROC/AUPrR curves at the sgRNA level.
-	sgRNA_level_recall.pdf = recall at 5% FDR at the sgRNA level.

In addition, Mel-202 and Mel-285 cell lines were included in the integrated Sanger + Broad dataset (Pacini et al. 2021) on the Project score. We compared the gene-level corrected logFCs we obtained with those in the integrated dataset. To obtain this integrated dataset, the authors applied batch correction with ComBat on the Sanger and Broad datasets, quantile-normalised them, and removed the first two principal components. Nevertheless, there is a high agreement between their corresponding gene essentiality profiles.
Results are summarised in these two scatter plots: Mel-202.pdf (Pearson's R = 0.7915) and Mel-285.pdf (0.7086).

For instance, considering the cell line Mel202, inside the folder you’ll find this structure:
```
├── CCR_correction
│   ├── 1.pdf
│   ├── 10.pdf
│   ├── 11.pdf
│   ├── 12.pdf
│   ├── 13.pdf
│   ├── 14.pdf
│   ├── 15.pdf
│   ├── 16.pdf
│   ├── 17.pdf
│   ├── 18.pdf
│   ├── 19.pdf
│   ├── 2.pdf
│   ├── 20.pdf
│   ├── 21.pdf
│   ├── 22.pdf
│   ├── 23.pdf
│   ├── 24.pdf
│   ├── 3.pdf
│   ├── 4.pdf
│   ├── 5.pdf
│   ├── 6.pdf
│   ├── 7.pdf
│   ├── 8.pdf
│   └── 9.pdf
├── Impact_on_phenotype.R
├── Impact_on_phenotype.pdf
├── Mel202_bf.tsv
├── Mel202_bf_scaled.tsv
├── Mel202_correctedCounts.RData
├── Mel202_corrected_counts.tsv
├── Mel202_corrected_logFCs.tsv
├── Mel202_corrected_segments.tsv
├── Mel202_count_norm.tsv
├── Mel202_fcs.pdf
├── Mel202_foldChanges.RData
├── Mel202_logFCs.tsv
├── Mel202_normCounts.RData
├── Mel202_normCounts.pdf
├── Mel202_raw_counts.tsv
├── Mel202_segments.tsv
├── gene_PrRc.pdf
├── gene_ROC.pdf
├── gene_signatures.pdf
├── mageck_corrected
│   ├── Mel202.gene_summary.txt
│   ├── Mel202.log
│   └── Mel202.sgrna_summary.txt
├── mageck_uncorrected
│   ├── Mel202.gene_summary.txt
│   ├── Mel202.log
│   └── Mel202.sgrna_summary.txt
├── sgRNA_PrRc.pdf
└── sgRNA_ROC.pdf
```

-	CCR_correction is a folder containing one plot per chromosome, with segments of sgRNAs’ equal log fold-change before and after the correction.
-	Impact_on_phenotype.pdf is a two-page pdf showing:
-	4 pie charts showing the CCR impact (GoF genes classified as LoF and vice versa) and distortion (original GoF and LoF genes classified as null phenotype).
-	1 barplot summarising the findings.
-	Impact_on_phenotype.R is the corresponding R object.
-	Mel202_bf.tsv is the BFs matrix outputted from BAGEL (used to derive the pooled BFs unscaled matrix).
-	Mel202_bf_scaled.tsv is the scaled BFs matrix at 5% FDR (used to derive the pooled BFs scaled matrix).
-	Mel202_corrected_counts.tsv is the matrix of normalised corrected sgRNAs’ treatment counts. 
-	Mel202_correctedCounts.RData is the corresponding R object.
-	Mel202_corrected_logFCs.tsv is the matrix of corrected logFCs (used to derive the pooled sgRNA- and gene-level logFCs matrices).
-	Mel202_corrected_segments.tsv is a matrix of genomic segments of equal sgRNAs’ logFCs (post-CRISPRcleanR correction).
-	Mel202_count_norm.tsv = normalised sgRNAs counts, scaled by the total number of reads per sample.
-	Mel202_fcs.pdf shows boxplots of sgRNAs logFCs before CRISPRcleanR correction.
-	Mel202_foldChanges.RData is the corresponding R object.
-	Mel202_logFCs.tsv shows sgRNAs logFCs derived from the normalised sgRNA counts (prior CRISPRcleanR correction).
-	Mel202_normCounts.pdf shows boxplots of sgRNA raw (left) and normalised (right) counts.
-	Mel202_normCounts.RData is the corresponding R object.
-	Mel202_raw_counts.tsv is the original matrix of sgRNAs raw read counts.
-	Mel202_segments.tsv is a matrix of genomic segments of sgRNAs’ equal logFCs identified by CRISPRcleanR with an indication of the type of correction to be applied (more details in (Iorio et al. 2018)).
-	gene_PrRc.pdf shows the precision-recall curve at the gene-level logFCs using reference essential and nonessential genes from (https://github.com/hart-lab/bagel) as positive and negative controls.
-	gene_ROC.pdf shows the ROC curve at the gene-level logFCs using reference essential and nonessential genes from (https://github.com/hart-lab/bagel) as positive and negative controls.
-	gene_signatures.pdf is the visualisation of a gene essentiality profile derived from applying CRISPRcleanR to the data, with superimposed reference positive control gene sets. The dashed red line indicates the depletion rank threshold at a fixed false discovery rate (i.e. 5%, computed as described above).
-	mageck_uncorrected is a folder containing the output of running MAGeCK MLE on the uncorrected normalised sgRNAs’ treatment counts.
-	mageck_corrected is a folder containing the output of running MAGeCK MLE on the corrected normalised sgRNAs’ treatment counts.
-	sgRNA_PrRc.pdf shows the precision-recall curve at the sgRNA-level logFCs using reference essential and nonessential genes from (https://github.com/hart-lab/bagel) as positive and negative controls.
-	sgRNA_ROC.pdf shows the ROC curve at the sgRNA-level logFCs using reference essential and nonessential genes from (https://github.com/hart-lab/bagel) as positive and negative controls.



## NOTE
The summary plot sgRNA_level_recall.pdf shows the recall at 5% FDR from sgRNA logFC profiles. In this plot, the cell line MP46 has a very low score (close to 0). This is due to few guides targeting non-essential genes (negative controls) with strong negative logFCs, which impairs the initial precision before recovering. The situation improves when looking at the gene logFC profiles. Nevertheless, this cell line has the lowest metrics across the quality assessments we performed, implying a poorer separation between core-fitness essential and non-essential genes.
