1. Overview
This repository contains code for genome-wide fitness and enrichment analysis using sequence-barcoded T-DNA insertions (RB-TDNAseq) for functional genomics in fungi as described in Functional genomics of lipid metabolism in the oleaginous yeast Rhodosporidium toruloides (https://www.biorxiv.org/content/early/2017/09/19/190546).

This analysis extends previous tools used to estimate gene fitness in bacteria using sequence-barcoded transposon insertion libraries (RB-TnSeq) as described in Rapid quantification of mutant fitness in diverse bacteria by sequencing randomly bar-coded transposons (Wetmore et al, PubMed ID 25968644).  Several steps in the following analysis use scripts developed for that study, or slightly modified versions.  Un-modified scripts are included in the directory Code/RB-TnSeq_scripts/ for convenience.  We used RB-TnSeq code as was available in early 2016.  An actively maintained codebase for that analysis is available at https://bitbucket.org/berkeleylab/feba.

There are two stages RB-TDNAseq, described in detail below: mapping sequence-barcoded T-DNA insertions in genes with a Tn-Seq like protocol and counting barcode abundance.

2. Prerequisites:
Various steps in this analysis use scripts written in Perl, Python and R.  They were developed on a system with Perl ver 5.10.1, Python 2.7.13 (within the Anaconda 4.4.0 scientific computing package), and R version 3.0.2.  The TnSeq mapping step requires command line access to BLAT.

3. Example Data
Due to space constrains on BitBucket, example high-throughput sequence data for TnSeq and BarSeq are stored externally in the NCBI Short Read Archive (SRA) at https://www.ncbi.nlm.nih.gov/sra, under the accession numbers SRP116146 and SRP116193, respectively.


4. Mapping T-DNA Insertions (TnSeq)
The TnSeq pipeline consists of the following steps:
1) MapTnSeq_Trimmed.pl: Extract T-DNA/genome junctions and map to the genome.
2) RandomPoolConcatamers.py: Find useful barcodes (unique, anamibuously mapping barcodes).
3) InsertionLocationJGI.py: Annotate insertion locations relative to gene models.
4) InsertBias.py: Produce summary stats on the distribution of useful T-DNA insertions.


An example batch script (TnSeq.batch) to run these steps on the example sequence data is included in the TnSeq_Results directory.

4.1 MapTnSeq_Trimmed.pl
Inputs:
-model [T-DNA_model.txt]: A text file specifying the sequence of TDNA from the TnSeq priming site to the junction with the genome on the first line, and the rest of the sequence for the binary plasmid used to introduce the T-DNA on the subsequent line.  In the first line the presence of a random sequence barcode is denoted by a sequence of N's of the appropriate length.  (see Reference_Files/TDNAmodels/model_pDP11BC for an example).  See the header of MapTnSeq_Trimmed.pl script for additional variables to adjust mapping behavior, such as the number of bases on either side of the barcode to match to the T-DNA model, ect.
-trim [integer]: The number of bases to "trim" from the end T-DNA model's first line before looking for T-DNA/genome junctions.  Because T-DNA insertion is often not precise at the right border sequence, many insertions will be missed if the search algorithm requires the full sequence be present up to the genome.
-first [TnSeq_reads.fastq]: Fastq file with reads from the TnSeq protocol 
-genome [genome.fasta]: Fasta file with all scaffolds

Outputs:
The script will output to STDout by default.  Redirect to a text file for use in later steps.  For each read the script output a line with the following tab-delimited fields:
*read id
*barcode sequence
*scaffold to which genomic sequence was mapped
*start position of match on the scaffold
*strand of the match 
*unique match flag: if not 1, then the genomic sequence had multiple possible matches in the genome
*position in the read (after removing T-DNA sequence) where the match starts.  If this position is greater than 1, the read may be a chimeric product. Alternatively, there may be some random filler sequence added at the T-DNA insertion site, as is known to happen with ATMT. In some implementations of TnSeq analysis these reads could be discarded as probable chimeric artifacts.
*bit score of the match
*percent identity of the match (should be near 100)

Example call:
MapTnSeq_Trimmed.pl -trim 10 -model T-DNA_model.txt -first TnSeq_reads.fastq -genome genome.fasta > mapped_reads.tab


4.2 RandomPoolConcatamers.py
Inputs:
--input [mapped_reads1.tab,mapped_reads2.tab]: comma separated list of output files from MapTnSeq_Trimmeed.pl.
--output [output_pool_name] : base name for the output files
--offbytwo : flag to filter barcodes within 10 base pairs of a more abundant barcode that differs by only 2-4 bases/indels (likely artifactual barcodes produced by sequencing errors) 
Outputs:
*output_pool_name : a tab-delimited table of useful barcodes in the mutant pool (barcodes that only map to one genomic location for which BLAT hits to the genome are unambiguous)
*output_pool_name_multilocus : a tab-delimited table of barcodes that mapped to two or more distinct genomic locations
*output_pool_name_ambiguous : a tab-delimited table of barcodes associated with reads for which BLAT returned multiple significant alignments to the genome and thus location of insertion is ambiguous.
*output_pool_name_offByOne : a tab-delimited table of barcodes for which a much more abundant barcode exists that only differs by one base. These are almost certainly sequencing errors. 
*output_pool_name_offbytwo : a tab-delimited table of barcodes for which a much more abundant barcode was mapped at the same location or within 10 basepairs and that more abundant barcode had a Levenshtien edit distance of less than 5 from the listed barcode.  Most likely these barcodes are artifacts of sequencing errors.
*output_pool_name_pastEnd : barcodes associated with reads that map only to sequence past the end of the T-DNA insertion, or for which the vast majority of the reads map "pastEnd".  That is the sequence matches another part of the T-DNA sequence or the binary vector.  These barcodes may be located in concatameric insertions for which the genomic junction is difficult to amplify in the TnSeq protocol, or insertions for which the T-DNA was aberrantly processed, or they could represent plasmids in un-cured agrobacterium cells.

The format of the above tables is:
*barcode : barcode sequence
*rcbarcode : reverse compliment of the barcode
*nTot : total reads in which the barcode was present
*n : reads in which the barcode was present and the match to the genome started at the first base past T-DNA sequence
*scaffold : scaffold the barcode mapped to
*strand : strand of the match
*pos : starting position on the scaffold
*type : type of insertion inferred from all the reads with this barcode.  'Single' insertions had only reads mapping to the genome.  'Concat' insertions had a mix of reads mapping to the genome and to the T-DNA insertion suggesting a concatameric insertion. 'InCat' insertions are concatameric insertions in which the individual T-DNA's are inserted in opposing directions.
*nMainLocation : number of reads mapping to the most common genomic location
*nTot : total reads 
*nPastEnd : number of reads mapping to T-DNA
*All genomic mappings : a list of genomic locations, strand, and number of reads at that location. Can be useful for parsing true multi-locus insertions from inverted concatamers and artifactual small changes in mapped locations due to sequencing errors.
*All pastEnd mappings : a list of location in the past-end sequence from the T-DNA model file and number of reads at that location. Can be useful to understand the nature of the insertion (concatamer, inverted concatamer, etc)

Example call:
RandomPoolConcatamers.py --input [mapped_reads1.tab,mapped_reads2.tab]  --output [output_pool_name] --offbytwo


4.3 InsertionLocationJGI.py
Inputs:
--pool_file [output_pool_name] : pool of useful barcodes from RandomPoolConcatamers.py
--output [output_pool_name] : base name for the output file 
--gff [genes.gff3] : gene model file. This script was written specifically for the gff3 files produced by the JGI gene annotation pipeline.  Because gff3 is a loose standard with many variants, this script may need to be modified to correctly interpret gffs from other sources.

Output:
*output_pool_name_genes.txt : a modified version of the input pool file with the following columns added:
**nearest_gene : proteinID of nearest gene to the insertion
**location : 'intergenic_promoter', '5_prime_UTR', 'coding_region', '3_prime_UTR', or 'intergenic_terminator'. Note for that 'coding_region' in referee to the entire region between the start and stop codons, including introns.
**CDS_fraction : -1 if outside the coding region otherwise 1-100, the percent of bases from start to stop
**Basepair_from_CDS : base pairs before the start codon (negative numbers) or after the stop codon (positive numbers)
**exon : 'exon' or 'intron' denoting location of the insertion in an exon or intron.  If the insertion mapped outside the coding region; 'intergenic', '5_prime_UTR', or '3_prime_UTR'.
   
Example call:
InsertionLocationJGI.py --pool_file output_pool_name_genes.txt --output output_pool_name --gff [genes.gff3]


4.4 InsertBias.py
Inputs:
--input [output_pool_name_genes.txt] : output from InsertionLocationJGI.py
--fasta [genome.fasta] : Fasta file with all scaffolds
--genes [genes.gff3] : gene model file. This script was written specifically for the gff3 files produced by the JGI gene annotation pipeline.  Because gff3 is a loose standard with many variants, this script may need to be modified to correctly interpret gffs from other sources.

Outputs:
*first_scaffold_density.pdf : plot of observed insertion density (inserts per kb) over the length of the first (presumably biggest) scaffold and simulated density from random insertions weighted towards intergenic/promoter/intron/exon locations as observed from biases in the data.
*GChistogram.pdf : histogram of local GC content for insertion sites vs the whole genome
*AllChrom_N.pdf : plot of number of inserts versus scaffold length

Example call:
InsertBias.py --input "$resultsDir$poolBaseName"_genes.txt --fasta "$genomeSequence" --genes "$genesGFF"


5. Fitness Estimation (BarSeq)
The BarSeq pipeline consists of the following steps:
1) MultiCodes_Variable_Length.pl: Count up barcodes in individual BarSeq samples
2) combineBarSeq.pl: Combine barcode counts across samples 
3) Wilcoxon.py: Perform wilcoxon singed rank test between specified conditions
4) BarSeqR_modified.pl: Calculate fitness scores and T-like test statistics (T-stats)
5) AverageReplicates.py: Combine fitness scores and T-stats across biological replicates
6) ResultsSummary.py: Produce a summary file with fitness scores and T-stats for selected conditions, as well as p-values (from wilcoxon signed rank tests) and T-stats for inter-condition comparisons.

An example batch script (BarSeq.batch) to run these steps on the example sequence data is included in the BarSeq_Results directory.

5.1 MultiCodes_Variable_Length.pl
Inputs:
-index [BarSeq_Sample] : A unique name for the BarSeq sample.  This string will be used as the header row in the output files here and from there in the 'poolcount' file produced in the next step.  It should match the "index" field in the metadata file used in step 5.4 as input into the R code from Wetmore at al.
-nPreExpected 15 : an integer specifying the number of bases in each read expected before the defined sequence preceding the barcode
-preseq AGCGTACG : a defined sequence expected before the random barcode
-postseq AGAG : a defined sequence expected after the random barcode
-out [BarSeq_Sample] : basename for the output files
< [BarSeq_Sample.fastq] : The BarSeq sequence file to be analyzed

Outputs:
*BarSeq_Sample.codes : number of times each barcode was seen
*BarSeq_Sample.counts : summary distribution of barcode abundance
*BarSeq_Sample.close : list of barcode pairs off by one nucleotide


Example call:
MultiCodes_Variable_Length.pl -index [BarSeq_Sample] -nPreExpected 15 -preseq AGCGTACG -postseq AGAG -out [BarSeq_Sample] < [BarSeq_Sample.fastq]


5.2 combineBarSeq.pl
Inputs:
[poolcount_prefix] : prefix for the the output file 
[poolfile] : Output of InsertionLocationJGI.py or RandomPoolConcatamers.py tab-delimited file of mapped barcodes
[BarSeq_Sample1.codes BarSeq_Sample2.codes etc]: list of .codes outputs from MultiCodes_Varialbe_Length.pl, separated with spaces


Outputs:
*poolcount_prefix.poolcount : table of mapped barcodes and counts across samples

Example call:
combineBarSeq.pl [poolcount_prefix] [poolfile] [BarSeq_Sample1.codes BarSeq_Sample2.codes etc]

5.3 Wilcoxon.py
Inputs:
--poolcount [poolcount_prefix.poolcount] : poolcount file from combineBarSeq.pl
--barcodeLocations [poolfile] : Output of InsertionLocationJGI.py, tab-delimited file of mapped barcodes with information on nearest genes
--output [output_prefix] : name for output file
--conditions "Condition1:[C1_A,C1_B,C1_C],Condition2:[C2_A,C2_B,C2_C],Condition3:[C3_A,C3_B,C3_C]" : names of conditions and samples from independent biological replicates for each condition. For statistical comparisons to be valid each biological replicate must be explicitly paired between replicates. For example if a one starter culture is used to seed C1_A, C2_A, and C3_A while a different starter culture is used to seed C1_B, C2_B, and C3_B and a third culture is used to seed C1_C, C2_C, and C3_C.
--tests "[[Condition1,Condition2],[Condition1,Condition3]]" : String formatted as a nested list of conditions to be compared.  In this example Condition1 will be compared to both Condition2 and Condition3.

Outputs:
*output_prefix_wilcoxon.txt : table of genes and p-values between conditions from the Wilcoxon signed rank test.


Example call:
Wilcoxon.py --poolcount [poolcount_prefix.poolcount] --barcodeLocations [poolfile] --output [wilcoxon_output] --conditions "Condition1:[C1_A,C1_B,C1_C],Condition2:[C2_A,C2_B,C2_C],Condition3:[C3_A,C3_B,C3_C]" --tests "[[Condition1,Condition2],[Condition1,Condition3]]"


5.4 BarSeqR_modified.pl
Inputs:
-org [experiment_name] : short descriptive name for the experiment
-indir [inputdirectory] : directory to search for input files
-exps [experiment_metadata] : tab-delimited file with experimental metadata
-outdir [replicate_results] : a directory to write the results
-poolfile [poolfile] : Output of InsertionLocationJGI.py or RandomPoolConcatamers.py tab-delimited file of mapped barcodes
-genesfile [genes.GC] : A tab-delimited file with gene names and locations

Outputs:
Directory [replicate_results] with the outputs described in Wetmore et al 2015.  Most important in this pipeline are the files fit_logratios.tab (fitness scores) and fit_t.tab (T-like test statistics).  

Example call:
BarSeqR_modified.pl -org [experiment_name] -indir [inputdirectory] -exps [experiment_metadata] -outdir [replicate_results] -poolfile [poolfile] -genesfile [genes.GC]

More information on the experimental metafile.
It is critical for this step that the metadata file be formatted as the R code expects.  First of all make sure to convert any non-unix newline characters to '\n'.  Handler code does not correctly interpret typical line breaks from excel, etc.  See the example files for the metadata layout.  The key columns are 'Date_pool_expt_started', 'Description', and 'Index'.  The 'Date_pool_expt_started' field is used to match experimental conditions with their 'Time0' starter cultures.  Because Wetmore et al designed their software for large-scale surveys of many conditions with many species, it does not explicitly handle experimental replicates.  All Time0 samples from the same day are combined and used as the reference samples for all other samples.  To ensure that all experimental samples are compared against the Time0 samples from which they originated, and thus ensure maximal statistical power make sure that all samples coming from a given Time0 sample and that Time0 sample have the same date.  Enter a different date for the next set of samples and their common Time0 sample, etc and etc.  The 'Description' field can be any short descriptive name for the individual conditions, but all time zero samples must be called Time0. The 'Index' fields must be the same short name as in the column headers in the poolcount file produced by combineBarSeq.pl 


5.5 AverageReplicates.py
Inputs:
--inputdir [replicate_results] : directory of results from individual replicate cultures (-outdir from BarSeqR)
--conditions "Condition1:[C1_A,C1_B,C1_C],Condition2:[C2_A,C2_B,C2_C],Condition3:[C3_A,C3_B,C3_C]" : names of conditions and samples from independent biological replicates for each condition. For statistical comparisons to be valid each biological replicate must be explicitly paired between replicates. For example if a one starter culture is used to seed C1_A, C2_A, and C3_A while a different starter culture is used to seed C1_B, C2_B, and C3_B and a third culture is used to seed C1_C, C2_C, and C3_C.
--outputdir [averaged_results] : directory to write results combined across replicates

Outputs:
Combines outputs from BarSeqR across replicates and writes them to the results directory. The key outputs are:
*fit_logratios.tab : all fitness scores
*fit_t.tab : T-like test statistics (versus T0) 
*strain_fit.tab : strain-level fitness scores

Example call:
AverageReplicates.py --inputdir [replicate_results] --conditions "Condition1:[C1_A,C1_B,C1_C],Condition2:[C2_A,C2_B,C2_C],Condition3:[C3_A,C3_B,C3_C]" --outputdir [averaged_results]


5.6 ResultsSummary.py
Inputs:
--wilcoxonfile [output_prefix_wilcoxon.txt] : output from Wilcoxon.py
--tests "[[Condition1,Condition2],[Condition1,Condition3]]" : String formatted as a nested list of conditions to be compared.  In this example Condition1 will be compared to both Condition2 and Condition3.
--output [results_summary.txt] : name of output file
--rename "set[A-Z]" : regular expression to strip out setA, setB, etc from condition headers.  These short set identifiers are added by the BarSeqR script and are not used in this pipeline.

Outputs:
results_summary.txt : A tab-delimited table of average fitness scores for each condition, T-stats and P-values from Wilcoxon signed rank tests for each condition versus Time0 samples as well as for comparisons between conditions given given with the --tests flag.

Example call:
ResultsSummary.py --febadir [averaged_results] --wilcoxonfile [output_prefix_wilcoxon.txt] --output [results_summary.txt] --tests "[[Condition1,Condition2],[Condition1,Condition3]]" --rename "set[A-Z]" 


