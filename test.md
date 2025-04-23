
<a href="https://doi.org/10.5281/zenodo.15111568"><img src="https://zenodo.org/badge/957550864.svg" alt="DOI"></a>

# ddRADseq_reads
## The set of pipelines for reads preprocessing and deduplicating from ddRADseq experiments, helping getting more reads to be used by main STACKS (https://catchenlab.life.illinois.edu/stacks/) pipeline steps.

### The pipelines ddradseq_pre.bash and ddradseq_dedup.bash were written to help preprocessing reads from ddRADseq experiments using sequencing of double digested of genomic DNA inserted sorrounded by inline barcode on the P5 adaptor read and DBR region on the P7 adaptor read, the procedure modified from Schweyen et al.(DOI: 10.1086/BBLv227n2p146)
#### The new solution gave finally up to 10x more Stacks with Coverage 2-6 times higher than when preprocessing reads with the original proposed approach with STACKS's process_radtags and clone_filter. Comparison results in ComparisonFinal.ods (https://github.com/rafalwoycicki/ddRADseq_reads/blob/main/ComparisonsFinal.ods) file and in the Comparisons.md (https://github.com/rafalwoycicki/ddRADseq_reads/blob/main/Comparisons.md) file.

The 3 separate datasets to test the software together will all supplementary infrormation needed (experiment setup, sequencing adaptors, RE, barcodes and DBR used) were kindly delivered by Prof. Maciej Pabijan (Institute of Zoology and Biomedical Research, Faculty of Biology, Jagiellonian University, Kraków, Poland).

---

### Dependencies
Pipeline uses cutadapt 5.0 (https://doi.org/10.14806/ej.17.1.200), seqtk (https://github.com/lh3/seqtk), and awk scripting language. At this moment it is also using 'zcat' Linux command, so could be not working on MacOS. When using specific conda environment like in my environment for cutadapt 5.0 pipeline, one needs to run the script by: `source <file name>` to enable loading this specific environment. When not using conda environment one can comment out the line with `conda activate environment`.

---

# ddradseq_pre.bash
Pipeline inspired by process_radtags script (https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php) from the STACKS pipeline (https://catchenlab.life.illinois.edu/stacks/), outputting quality filtered, demultiplexed and cut site filtered paired and unpaired reads to be used by next STACKS steps.

There are several variables to be set on command line before running the script (the --help option prints them all):

- `iftest="-0"` — "1" for testing one barcode, "-0" for normal analysis of whole dataset
- `qualada="20,20"` — quality trimming adapters
- `qualfil="20,20"` — quality trimming filtering
- `errbar=1` — number of allowed errors for demultiplexing
- `errfil=0` — number of allowed errors for CUT site filtering
- `errfilanc=0` — number of allowed errors for anchored CUT site filtering
- `thr=40` — number of processor threads to use
- `lenada=140` — minimum length for adapter removal
- `lenfil=140` — minimum length for filtering cut sites
- `readsP5="1.fq.gz"` — reads from P5 primer - forward
- `readsP7="2.fq.gz"` — reads from P7 primer - reverse
- `barcodes="barcodes_P1_cutadapt.txt"` — fasta file with barcodes
- `directory="/mnt/qnap/projects/RafalWoycicki/"` — directory to save BIG data files, "" if locally

Remember:
- Set proper `PATH TO READS TEMPLATE` at `reads` variable
- Set proper `BARCODES` file name at `barcodes` variable
- Set proper `DIRECTORY` name at `directory` variable

#### Sequencing adapters:
- `A_p5_5p="P5read5prim=TACACGACGCTCTTCCGATCT"`
- `A_p5_3p="P5read3prim=AGATCGGAAGAGCACACGTCT"`
- `A_p7_5p="P7read5prim=AGACGTGTGCTCTTCCGATCT"`
- `A_p7_3p="P7read3prim=AGATCGGAAGAGCGTCGTGTA"`

#### Cutsite adapters:
- `C_p5_5p="P5read5primSBF1=^TGCAGG"`
- `C_p5_3p="P5read3primMSE1_DBR=TTAGCNNNNNNNN"`
- `C_p7_5p="P7read5primDBR_MSE1=^NNNNNNNNGCTAA"`
- `C_p7_3p="P7read3primSBF1=CCTGCA"`

The pipeline runs in 4 sections:

### 1. Removing sequencing adapters
Adapter sequences are trimmed when found in both reads.  
Input: paired reads  
Output: `Adapters.?.fq.gz`

### 2. Demultiplexing
Barcodes at the start of P5 reads are matched and trimmed.  
Output: `Barcodes.?.{barcode}.fq.gz`

### 3. Filtering cut sites
Reads are filtered based on cut site and DBR sequences.  
Output: `Filtered.?.{barcode}.fq.gz` and `Filtered.?.{barcode}.untrimmed.fq.gz`

### 4. Rescuing untrimmed reads
Attempts to rescue reads with valid cut sites not found in pairs.  
Output: `Rescued.?.{barcode}.fq.gz`

### Statistics
Outputs counts per step to `Counts.stat`.

---

# ddradseq_dedup.bash

Pipeline inspired by clone_filter (https://catchenlab.life.illinois.edu/stacks/comp/clone_filter.php) from STACKS package, for deduplication of reads using DBR in P7 reads.

Variables to set (soon to be command-line):
- `iftest="-0"`
- `directory="/mnt/qnap/projects/RafalWoycicki/"`
- `p5len=130`
- `p7len=140`
- `len=130`
- `thr=40`
- `err=0`
- `p7_5p_seq="P7read5primDBR_MSE1=^NNNNNNNNGC"`

### Steps:
1. Select reads of defined length
2. Form and sort DBR_P5read pairs, filter for uniqueness
3. Filter P7 reads by DBR_P7read uniqueness
4. Trim reads to final `len`

Output:  
`Short.Filtered.{barcode}.p5p7dedupl.?.fq`

### Statistics
Counts saved in `Counts_dedupl.stat`
