Overview
This repository contains code for genome-wide fitness and enrichment analysis using sequence-barcoded T-DNA insertions (RB-TDNAseq) for functional genomics in fungi as described in Functional genomics of lipid metabolism in the oleaginous yeast Rhodosporidium toruloides (https://www.biorxiv.org/content/early/2017/09/19/190546).

This analysis extends previous tools used to estimate gene fitness in bacteria using sequence-barcoded transposon insertion libraries (RB-TnSeq) as described in Rapid quantification of mutant fitness in diverse bacteria by sequencing randomly bar-coded transposons (Wetmore et al, PubMed ID 25968644).  Several steps in the following analysis use scripts developed for that study, or slightly modified versions.  Un-modified scripts are included in the directory Code/RB-TnSeq_scripts/ for convenience.  We used RB-TnSeq code as was available in early 2016.  An actively maintained codebase for that analysis is available at https://bitbucket.org/berkeleylab/feba.

There are two stages RB-TDNAseq, described in detail below: mapping sequence-barcoded T-DNA insertions in genes with a Tn-Seq like protocol and counting barcode abundance.

Prerequisites:
Various steps in this analysis use scripts written in Perl, Python and R.  They were developed on a system with Perl ver 5.10.1, Python 2.7.13 (within the Anaconda 4.4.0 scientific computing package), and R version 3.0.2.  

Example Data
Due to space constrains on BitBucket, example high-throughput sequence data for TnSeq and BarSeq are stored externally in the NCBI Short Read Archive (SRA) at https://www.ncbi.nlm.nih.gov/sra, under the accession numbers SRP116146 and SRP116193, respectively.


Mapping T-DNA Insertions (TnSeq)
The TnSeq pipeline consists of the following steps:
1) MapTnSeq_Trimmed.pl: Extract T-DNA/genome junctions and map to the genome.
2) RandomPoolConcatamers.py: Find useful barcodes (unique, anamibuously mapping barcodes).
3) InsertionLocationJGI.py: Annotate insertion locations relative to gene models.
4) InsertBias.py: Produce summary stats on the distribution of useful T-DNA insertions.


An example batch script (TnSeq.batch) to run these steps on the example sequence data is included in the TnSeq_Results directory.

Fitness Estimation (BarSeq)
The BarSeq pipeline consists of the following steps:
1) MultiCodes_Variable_Length.pl: Count up barcodes in individual BarSeq samples
2) combineBarSeq.pl: Combine barcode counts across samples 
3) Wilcoxon.py: Perform wilcoxon singed rank test between specified conditions
4) BarSeqR_modified.pl: Calculate fitness scores and T-like test statistics (T-stats)
5) AverageReplicates.py: Combine fitness scores and T-stats across biological replicates
6) ResultsSummary.py: Produce a summary file with fitness scores and T-stats for selected conditions, as well as p-values (from wilcoxon signed rank tests) and T-stats for inter-condition comparisons.

An example batch script (BarSeq.batch) to run these steps on the example sequence data is included in the BarSeq_Results directory.