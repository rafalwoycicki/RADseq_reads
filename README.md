<a href="https://doi.org/10.5281/zenodo.15111568"><img src="https://zenodo.org/badge/957550864.svg" alt="DOI"></a>

# ddRADseq_reads
pipeline for reads preprocessing from ddRADseq experiments done similar to one described in Schweyen et al. 2014 (DOI: 10.1086/BBLv227n2p146) - but can be freely modified for different needs.

pipeline inspired by process_radtags script from the STACKS pipeline (https://catchenlab.life.illinois.edu/stacks/), outputing quality filtered, demultiplexed and cut site filtered paired and unpaired reads to be used by next STACKS steps.

pipeline uses cutadapt 5.0 (https://doi.org/10.14806/ej.17.1.200).

when using specific conda envinorment like in my environment for cutadapt 5.0 pipeline, one needs to run the script by: source <<file name>> to enable loading this specific envinorment. When not using CONDA you can comment out the line with:
conda activate cutadapt

There are several variables to be set before running the script (will me soon modified to be set in command line.

iftest="-0" # "1" for testing one barcode, "-0" for normal analysis of whole dataset

qualada="20,20" # quality trimming adapters

qualfil="20,20" # quality trimming filtering

errbar=1 # number of allowed errors for demultiplexing

errfil=0 # number of allowed errors for CUT site filtering

errfilanc=0 # number of allowed error for anchored CUT site filtering

thr=40 # number of processor threads to use

lenada=140 # minimum length for adapter removal

lenfil=140 # minimum length for filtering cut sites

readsP5="1.fq.gz" # reads from P5 primer - forward

readsP7="2.fq.gz" # reads from P7 primer - reverse

barcodes="barcodes_P1_cutadapt.txt" # fasta file with barcodes

directory="/mnt/qnap/projects/RafalWoycicki/" # directory to save BIG data files, "" - if localy

Remember:

You need to input proper PATH TO READS TEMPLATE at "$reads" variable

You need to input proper BARCODES file name at "$barcodes" variable

You need to input proper DIRECTORY name at $directory variable

At this moment pipeline runs in 4 sections:

# 1. removing sequencing adapters:
uses adaptor sequences for both reads from both sites: 5' "TACACGACGCTCTTCCGATCT" and 3' "AGATCGGAAGAGCACACGTCT" for P5 reads as well as 5' "AGACGTGTGCTCTTCCGATCT" and 3' "AGATCGGAAGAGCGTCGTGTA" for P7 reads. The adaptor sequences are not required but when found are trimmed away. Uses as input paired reads.

# 2. demultiplexing:
uses barcode sequences anchored/required at the beginning of the 5' site of P5 read. Found barcodes are trimmed away. This step outputs paired reads sorted by the barcode found as well as untrimmed reads where barcode was not found. Uses as input the output of step 1.

# 3. filteringcutsites:
uses specific cut site sequences with DBR region. At this moment these are: 

For P5 reads: from the 5' site this is anchored/required at the beginning SBF1 cut site "^TGCAGG" and at the 3' end optional MSE1 cut site with DBR region "TTAGCNNNNNNNN".

For P7 reads: from the 5' site this is anchored/required at the beginning MSE1 site with DBR region "^NNNNNNNNGCTAA" and at the 3' end optional SBF1 cut site "CCTGCA".

The output of this step consists of paired end reads where both the anchored 5' cut site was found in reads P5 and anchored 5' cut site with DBR region was found in read P7.
Untrimmed reads are also placed in separate output files.

The cut sites and DBR regions are left intact and not trimmed away.

# 4. rescuing untrimmed reads:
this step tries to rescue P5 and P7 reads which in the previous step were in the output of untrimmed reads.

The script rescues P5 reads which does not have paired correct P7 read (therefore untrimmed by the previous step), but contain anchiored at the 5' of P5 read SBF1 cut site "^TGCAGG".

The script rescues P7 reads which does not have paired correct P5 read (therefore untrimmed by the previous step), but contain anchiored at the 5' of P7 read MSE1 cut site with DBR region "^NNNNNNNNGCTAA".

The cut sites and DBR regions are left intact and not trimmed away.

# Statistics.
The script outputs the statistics in the form of the number of the reads left for the analysis from each step in a file named "Counts.stat"







