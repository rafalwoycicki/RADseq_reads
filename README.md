<a href="https://doi.org/10.5281/zenodo.15111568"><img src="https://zenodo.org/badge/957550864.svg" alt="DOI"></a>

# ddRADseq_reads
## The set of pipelines for reads preprocessing and deduplicating from ddRADseq experiments, helping getting more reads to be used by main STACKS (https://catchenlab.life.illinois.edu/stacks/) pipeline steps.

### The pipelines `ddradseq_pre.bash` and `ddradseq_dedup.bash` and `cleanup_radseq.sh` were written to help preprocessing reads from ddRADseq experiments using sequencing of double digested genomic DNA insert sorrounded by inline barcode on the P5 forward adaptor read and DBR region on the P7 reverse adaptor read, the procedure modified from Schweyen et al. (`DOI: 10.1086/BBLv227n2p146`), as shown on this hand made schema (for other specific constructs one needs to modify the sequences of adapters):

<img src="https://github.com/rafalwoycicki/ddRADseq_reads/blob/main/construct_schema.jpg" width="500" />

#### The new solution gave finally up to 10x more Stacks with Coverage 2-6 times higher than when preprocessing reads with STACKS's process_radtags and clone_filter. Comparison results in ComparisonFinal.ods (https://github.com/rafalwoycicki/ddRADseq_reads/blob/main/ComparisonsFinal.ods) file and in the Comparisons.md (https://github.com/rafalwoycicki/ddRADseq_reads/blob/main/Comparisons.md) file.

The 3 separate datasets to test the software together will all supplementary infrormation needed (experiment setup, sequencing adaptors, RE, barcodes and DBR used) were kindly delivered by Prof. Maciej Pabijan (Institute of Zoology and Biomedical Research, Faculty of Biology, Jagiellonian University, Kraków, Poland).

---
### Dependencies

Pipeline uses:
- `cutadapt 5.0` (https://doi.org/10.14806/ej.17.1.200)
- `seqtk` (https://github.com/lh3/seqtk) - for ddradseq_dedup.bash
- `awk` scripting language

When using a conda environment, run the script using:
```bash
source <file name>
```
to enable the specific environment. If not using conda, comment out the line with `conda activate environment`.

Temporary files can be removed by running the script `cleanup_radseq.sh` (run with --help for instruction)

---
# ddradseq_pre.bash
pipeline inspired by process_radtags script (https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php) from the STACKS pipeline (https://catchenlab.life.illinois.edu/stacks/), outputing quality filtered, demultiplexed and cut site filtered paired and unpaired reads to be used by next STACKS steps.

There are several variables to be set on command line before running the script (the '--help' option prints them all):

```bash
iftest="1" # 1 for testing one barcode, "-0" for not testing
qualada="20,20" # quality trimming adapters
qualfil="20,20" # quality trimming filtering
errbar=0 # number of allowed errors for demultiplexing
errfil=0 # number of allowed errors for CUT site filtering
thr=10 # number of processor threads to use
lenada=140 # minimum length for adapter removal
lenfil=140 # minimum length for filtering cut sites
readsP5="1.fq.gz" # reads from P5 primer - forward
readsP7="2.fq.gz" # reads from P7 primer - reverse
barcodes="" # fasta file with barcodes
directory="" # output directory path (needs to end with "/"), set "" - if localy
cutsites="cutsites.txt" # REQUIRED: File with cutsite sequences
```
Remember:
- You need to input proper PATH TO READS TEMPLATE at "$reads" variable
- You need to input proper BARCODES file name at "$barcodes" variable
- You need to input proper DIRECTORY name at $directory variable

### Sequencing adapters. These are usually the same for most Illumina sequencing runs.
```bash
A_p5_5p="P5read5prim=TACACGACGCTCTTCCGATCT" # read P5 5prim sequencing adapter sequence
A_p5_3p="P5read3prim=AGATCGGAAGAGCACACGTCT" # read P5 3prim sequencing adapter sequence
A_p7_5p="P7read5prim=AGACGTGTGCTCTTCCGATCT" # read P7 5prim sequencing adapter sequence
A_p7_3p="P7read3prim=AGATCGGAAGAGCGTCGTGTA" # read P7 3prim sequencing adapter sequence
```

### Cutsites adapters in the <cutsites.txt> file.
These need to be adjusted to your specific RAD-seq experiment. Here filtering the citsites was done using linked cutadapt adapters option, so ^C_p5_5p...C_p5_3p for P5 read and ^C_p7_5p...C_p7_3p for P7 read, where both 5p sequences were required and 3p sequences optional. When linked adapters are used both required sequences needs to be found for pair to be retained.
```bash
  # read P5 5prim cut site for the SBF1 RE
C_p5_5p=P5read5primSBF1=^TGCAGG

# read P5 3prim cut site for the MSE1 RE including DBR region
C_p5_3p=P5read3primMSE1_DBR=TTAGCNNNNNNNN

# read P7 5prim cut site for the MSE1 RE including DBR region
C_p7_5p=P7read5primDBR_MSE1=^NNNNNNNNGCTAA

# read P7 3 prim cut site for the SBF1 RE
C_p7_3p=P7read3primSBF1=CCTGCA
```

### At this moment pipeline runs in 4 sections:

#### 1. Removing sequencing adapters:
Uses adaptor sequences for both reads from both sites: 5' "TACACGACGCTCTTCCGATCT" and 3' "AGATCGGAAGAGCACACGTCT" for P5 reads as well as 5' "AGACGTGTGCTCTTCCGATCT" and 3' "AGATCGGAAGAGCGTCGTGTA" for P7 reads. The adaptor sequences are not required but when found are trimmed away. Uses as input paired reads.

- Output file schema: Adapters.?.fq.gz

#### 2. Demultiplexing:
Uses barcode sequences anchored/required at the beginning of the 5' site of P5 read. Found barcodes are trimmed away. This step outputs paired reads sorted by the barcode found as well as untrimmed reads where barcode was not found. Uses as input the output of step 1.

- Output file schema: Barcodes.?.{barcode}.fq.gz

#### 3. Filteringcutsites:
Uses specific cut site sequences with DBR region. At this moment these are: 

- For P5 reads: from the 5' site this is anchored/required at the beginning SBF1 cut site "^TGCAGG" and at the 3' end optional MSE1 cut site with DBR region "TTAGCNNNNNNNN".
- For P7 reads: from the 5' site this is anchored/required at the beginning MSE1 site with DBR region "^NNNNNNNNGCTAA" and at the 3' end optional SBF1 cut site "CCTGCA".
- The output of this step consists of paired end reads where both the anchored 5' cut site was found in reads P5 and anchored 5' cut site with DBR region was found in read P7.
- Untrimmed reads are also placed in separate output files.
- The cut sites and DBR regions are left intact and not trimmed away.
- Output file schema: Filtered.?.{barcode}.fq.gz ( and Filtered.?.{barcode}.untrimmed.fq.gz )

#### 4. Rescuing untrimmed reads:
This step tries to rescue P5 and P7 reads which in the previous step were in the output of untrimmed reads.

- The script rescues P5 reads which does not have paired correct P7 read (therefore untrimmed by the previous step), but contain anchored SBF1 cut site "^TGCAGG" at the 5' of P5 read .
- The script rescues P7 reads which does not have paired correct P5 read (therefore untrimmed by the previous step), but contain anchored MSE1 cut site with DBR region "^NNNNNNNNGCTAA" at the 5' of P7 read .
- The cut sites and DBR regions are left intact and not trimmed away.
- Output file schema: Rescued.?.{barcode}.fq.gz

### Statistics.
The script outputs the statistics in the form of the number of the reads left for the analysis from each step in a file named "Counts.stat"

---

# ddradseq_dedup.bash

Pipeline inspired by `clone_filter` script (https://catchenlab.life.illinois.edu/stacks/comp/clone_filter.php) from STACKS package for deduplication of reads based on DBR region in P7 reads.

There are several variables to be set before running the script in the command line:
```bash
iftest="1" # 1 if test "-0" if not test
directory="" # output directory path (needs to end with "/"), set "" - if localy
p5len=130 #length of p5 read with SBF1 cut site beginning with ^TGCA
p7len=140 #length of p7 read with MSE1 cut site and DBR
len=130 # length of the final reads
thr=10 # number of processor therds to use for cutadapt
err=0 # numner of mismatches allowed adapter sequence for cutadapt
barcodes="" # fasta file with barcodes
dbr_sequence="NNNNNNNN" # DBR sequence in nucleotides
dbr_pattern="[AGCT]{8}" # the DBR regex pattern
motif_cut_adapter="GC" # the outer part of the adapter adjacent to the CUT site and DBR site from P7 adapter
motif_cut_rerest="TAA" # the inner part of the CUT site from P7 adapter site which stays with insert
```

- This script takes as input the paired sequences after the cut sites filtering step of the ddradseq_pre.bash script in the form of "Filtered.?.{barcode}.fq.gz" files.
- Additionally the script uses Rescued.2.{barcode}.fq.gz files from the previous step to deduplicate rescued reverse reads (as these contain DBR sequence).

### Algorithm:
1. First only the P5 reads of the $p5len and P7 reads of the $p7len are considered and both the reads are trimmed to these lengths respectively.
2. The DBR sequences of the p7 reads are combined with the p5 reads sequences and these combined sequences of P5 reads (DBR_P5reads) are sorted.
   - After sorting only the unique sequences are left after comparison (one sequence from the multiplicated sequences is retained randomly).
   - The unique genome sequnces are the ones existing only once with the specific DBR. So the additional "DBR_P5reads" reads are considered PCR duplicates and purged.
3. As in practice the paired P7 reads of the P5 reads are not all unique, this step double check that the p7 reads are unique (by using DBR_P7reads combination) and paired with the unique p5 reads.
4. Shortening all the reads to the final $len length as required by ustacks (https://catchenlab.life.illinois.edu/stacks/comp/ustacks.php) part of STACKS package. 

- The same procedure is run for rescued reverse (P7) reads in Rescued.2.{barcode}.fq.gz files

- The output file schema of the files ready to be used by "ustacks" are: Short.Filtered.{barcode}.p5p7dedupl.?.fq and Short.Rescued.2.{barcode}.fq.gz

### Statistics.
The script outputs the statistics in the form of the number of the reads left for the analysis from each step in a file named "Counts_dedupl.stat"





