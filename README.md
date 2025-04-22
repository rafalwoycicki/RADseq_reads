<a href="https://doi.org/10.5281/zenodo.15111568"><img src="https://zenodo.org/badge/957550864.svg" alt="DOI"></a>

# ddRADseq_reads
pipeline for reads preprocessing from ddRADseq experiments done similar to one described in Schweyen et al. 2014 (DOI: 10.1086/BBLv227n2p146) - but can be freely modified for different needs.

## The pipelines ddradseq_pre.bash and ddradseq_dedup.bash were written to help preprocessing reads from ddRADseq experiments using sequencing of double digested of genomic DNA inserted sorrounded by inline barcode on the P5 adaptor read and DBR region on the P7 adaptor read, the procedure modified from Schweyen et al.
## The new solution gave finally up to 10x more Stacks with Coverage 2-6 times higher than when preprocessing reads with the original proposed approach with STACKS's process_radtags and clone_filter. Comparison results in ComparisonFinal.ods file and in the attached graphs.

Old dataset 130nt All		p5_cpa	p7_cpa	p5_rnp	p7_rnp	p5_dpa	p7_dpa		Old dataset 130nt All		Stacks	Mean_cov
process_radtags eb2 ea2	Average	1,11%	1,11%	7,86%	0,23%	0,49%	0,49%		process_radtags eb2 ea2	Average	120675,625	10,51
ddradseq_reads eb1	Average	7,11%	7,11%	0,46%	0,22%	2,56%	2,56%		ddradseq_reads eb1	Average	131449,75	62,66

pipeline uses cutadapt 5.0 (https://doi.org/10.14806/ej.17.1.200).

when using specific conda envinorment like in my environment for cutadapt 5.0 pipeline, one needs to run the script by: source <<file name>> to enable loading this specific envinorment. When not using CONDA you can comment out the line with:
conda activate cutadapt

## ddradseq_pre.bash
pipeline inspired by process_radtags script (https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php) from the STACKS pipeline (https://catchenlab.life.illinois.edu/stacks/), outputing quality filtered, demultiplexed and cut site filtered paired and unpaired reads to be used by next STACKS steps.

There are several variables to be set before running the script (will me soon modified to be set in command line.

#### iftest="-0" # "1" for testing one barcode, "-0" for normal analysis of whole dataset
#### qualada="20,20" # quality trimming adapters
#### qualfil="20,20" # quality trimming filtering
#### errbar=1 # number of allowed errors for demultiplexing
#### errfil=0 # number of allowed errors for CUT site filtering
#### errfilanc=0 # number of allowed error for anchored CUT site filtering
#### thr=40 # number of processor threads to use
#### lenada=140 # minimum length for adapter removal
#### lenfil=140 # minimum length for filtering cut sites
#### readsP5="1.fq.gz" # reads from P5 primer - forward
#### readsP7="2.fq.gz" # reads from P7 primer - reverse
#### barcodes="barcodes_P1_cutadapt.txt" # fasta file with barcodes
#### directory="/mnt/qnap/projects/RafalWoycicki/" # directory to save BIG data files, "" - if localy

Remember:

You need to input proper PATH TO READS TEMPLATE at "$reads" variable

You need to input proper BARCODES file name at "$barcodes" variable

You need to input proper DIRECTORY name at $directory variable


#### Sequencing adapters:
#### A_p5_5p="P5read5prim=TACACGACGCTCTTCCGATCT" # read P5 5prim sequencing adapter sequence
#### A_p5_3p="P5read3prim=AGATCGGAAGAGCACACGTCT" # read P5 3prim sequencing adapter sequence
#### A_p7_5p="P7read5prim=AGACGTGTGCTCTTCCGATCT" # read P7 5prim sequencing adapter sequence
#### A_p7_3p="P7read3prim=AGATCGGAAGAGCGTCGTGTA" # read P7 3prim sequencing adapter sequence

#### Cutsites adapters:
#### C_p5_5p="P5read5primSBF1=^TGCAGG" # read P5 5prim cut site for the SBF1 RE
#### C_p5_3p="P5read3primMSE1_DBR=TTAGCNNNNNNNN" # read P5 3prim cut site for the MSE1 RE including DBR region
#### C_p7_5p="P7read5primDBR_MSE1=^NNNNNNNNGCTAA" # read P7 5prim cut site for the MSE1 RE including DBR region
#### C_p7_3p="P7read3primSBF1=CCTGCA" # read P7 3 prim cut site for the SBF1 RE

At this moment pipeline runs in 4 sections:

### 1. removing sequencing adapters:
uses adaptor sequences for both reads from both sites: 5' "TACACGACGCTCTTCCGATCT" and 3' "AGATCGGAAGAGCACACGTCT" for P5 reads as well as 5' "AGACGTGTGCTCTTCCGATCT" and 3' "AGATCGGAAGAGCGTCGTGTA" for P7 reads. The adaptor sequences are not required but when found are trimmed away. Uses as input paired reads.

Output file schema: Adapters.?.fq.gz

### 2. demultiplexing:
uses barcode sequences anchored/required at the beginning of the 5' site of P5 read. Found barcodes are trimmed away. This step outputs paired reads sorted by the barcode found as well as untrimmed reads where barcode was not found. Uses as input the output of step 1.

Output file schema: Barcodes.?.{barcode}.fq.gz

### 3. filteringcutsites:
uses specific cut site sequences with DBR region. At this moment these are: 

For P5 reads: from the 5' site this is anchored/required at the beginning SBF1 cut site "^TGCAGG" and at the 3' end optional MSE1 cut site with DBR region "TTAGCNNNNNNNN".

For P7 reads: from the 5' site this is anchored/required at the beginning MSE1 site with DBR region "^NNNNNNNNGCTAA" and at the 3' end optional SBF1 cut site "CCTGCA".

The output of this step consists of paired end reads where both the anchored 5' cut site was found in reads P5 and anchored 5' cut site with DBR region was found in read P7.
Untrimmed reads are also placed in separate output files.

The cut sites and DBR regions are left intact and not trimmed away.

Output file schema: Filtered.?.{barcode}.fq.gz ( and Filtered.?.{barcode}.untrimmed.fq.gz )

### 4. rescuing untrimmed reads:
this step tries to rescue P5 and P7 reads which in the previous step were in the output of untrimmed reads.

The script rescues P5 reads which does not have paired correct P7 read (therefore untrimmed by the previous step), but contain anchiored at the 5' of P5 read SBF1 cut site "^TGCAGG".

The script rescues P7 reads which does not have paired correct P5 read (therefore untrimmed by the previous step), but contain anchiored at the 5' of P7 read MSE1 cut site with DBR region "^NNNNNNNNGCTAA".

The cut sites and DBR regions are left intact and not trimmed away.

Output file schema: Rescued.?.{barcode}.fq.gz

### Statistics.
The script outputs the statistics in the form of the number of the reads left for the analysis from each step in a file named "Counts.stat"

## ddradseq_dedup.bash

pipeline inspired by clone filter script (https://catchenlab.life.illinois.edu/stacks/comp/clone_filter.php) from STACKS package for deduplication of reads based on DBR region in P7 reads.

There are several variables to be set before running the script (will me soon modified to be set in command line.

#### iftest="-0" # "1" if test "-0" if not test
#### directory="/mnt/qnap/projects/RafalWoycicki/" # "" - if localy
#### p5len=130 #length of p5 read with SBF1 cut site beginning with ^TGCA
#### p7len=140 #length of p7 read with MSE1 cut site and DBR.
#### len=130 # lengthof the final reads
#### thr=40 # nymber of processor therds to use for cutadapt
#### err=0 # numner of mismatches allowed adapter sequence for cutadapt
#### p7_5p_seq="P7read5primDBR_MSE1=^NNNNNNNNGC" # name and sequence of the nucleotides filtered by the final reads shortening cutadapt script (for details check the cutadapt manual)

This script takes as input the paired sequences after the cut sites filtering step of the ddradseq_pre.bash script in the form of "Filtered.?.{barcode}.fq.gz" files.

Algorithm:
#### 1st: first only the P5 reads of the $p5len and P7 reads of the $p7len are considered and both the reads are trimmed to these lengths respectively
#### 2nd: the DBR sequence of the p7 reads is combined with the p5 read sequence and this combined sequences of P5 reads: DBR_P5read are sorted and only the unique sequences are left after comparison (one sequence from the multiplicated sequences is taken randomly). The unique genome sequnces are the ones existing only once with the specific DBR. So the additional "reads" being combination of DBR and P5read are considered PCR duplicates and not purged.
#### 3rd: As in practice the paired P7 reads of the P5 read are not all unique (DBR_P7read), this step filters the p7 reads to be the same as p5 reads and to be unique.
#### 4th: shortening all the reads to the final $len length as required by ustacks (https://catchenlab.life.illinois.edu/stacks/comp/ustacks.php )part of STACKS package. 

The output file schema of the files ready to be used by "ustacks" is: Short.Filtered.{barcode}.p5p7dedupl.?.fq

### Statistics.
The script outputs the statistics in the form of the number of the reads left for the analysis from each step in a file named "Counts_dedupl.stat"





