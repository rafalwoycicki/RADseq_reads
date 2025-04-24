#!/bin/bash

### !!! need to run script using: source <<file name>>
conda activate cutadapt

# pipeline for trimming and demultiplexing illumina reads after RADseq
# using cutadapt 5.0

#if test variable
iftest="1" # 1 for testing one barcode, "-0" for not testing
# general variables
qualada="20,20" # quality trimming adapters
qualfil="20,20" # quality trimming filtering
errbar=1 # number of allowed errors for demultiplexing
errfil=0 # number of allowed errors for CUT site filtering
thr=10 # number of processor threads to use
lenada=140 # minimum length for adapter removal
lenfil=140 # minimum length for filtering cut sites
#crucial variables
readsP5="1.fq.gz" # reads from P5 primer - forward
readsP7="2.fq.gz" # reads from P7 primer - reverse
barcodes="" # fasta file with barcodes
directory="" # "" - if localy
cutsites="" # REQUIRED: File with cutsite sequences
#adapters
A_p5_5p="P5read5prim=TACACGACGCTCTTCCGATCT" # read P5 5prim sequencing adapter sequence
A_p5_3p="P5read3prim=AGATCGGAAGAGCACACGTCT" # read P5 3prim sequencing adapter sequence
A_p7_5p="P7read5prim=AGACGTGTGCTCTTCCGATCT" # read P7 5prim sequencing adapter sequence
A_p7_3p="P7read3prim=AGATCGGAAGAGCGTCGTGTA" # read P7 3prim sequencing adapter sequence

# Help message
usage() {
  echo "Usage: source ddradseq_pre.bash [options]"
  echo "Options:"
  echo "  --iftest <value>         Set test mode: 1 for testing one barcode (default), -0 for not testing"
  echo "  --qualada <value>        Quality for adapter trimming (default: 20,20)"
  echo "  --qualfil <value>        Quality for filtering (default: 20,20)"
  echo "  --errbar <value>         Allowed errors demultiplexing (default: 1)"
  echo "  --errfil <value>         Allowed errors cutsite filtering (default: 0)"
  echo "  --thr <value>            Number of threads (default: 10)"
  echo "  --lenada <value>         Min length after adapter removal (default: 140)"
  echo "  --lenfil <value>         Min length after cutsite filtering (default: 140)"
  echo "  --readsP5 <file>         Forward reads file (default: 1.fq.gz)"
  echo "  --readsP7 <file>         Reverse reads file (default: 2.fq.gz)"
  echo "  --barcodes <file>        Barcodes file (fasta file with barcodes) (default not set)"
  echo "  --directory <path>       Project output directory set to blank/nothing if locally (default)"
  echo "  --A_p5_5p <value>        P5 5' adapter sequence (default: P5read5prim=TACACGACGCTCTTCCGATCT)"
  echo "  --A_p5_3p <value>        P5 3' adapter sequence (default: P5read3prim=AGATCGGAAGAGCACACGTCT)"
  echo "  --A_p7_5p <value>        P7 5' adapter sequence (default: P7read5prim=AGACGTGTGCTCTTCCGATCT)"
  echo "  --A_p7_3p <value>        P7 3' adapter sequence (default: P7read3prim=AGATCGGAAGAGCGTCGTGTA)"
  echo "  --cutsites <file>        REQUIRED: Line by line list of 4 cutsite sequences (C_p5_5p=, C_p5_3p=, C_p7_5p=, C_p7_3p=) to be used as linked adapters by CUTADAPT"
  echo "  --help                   Show this help message"
}

# PROTECTION: Require at least one argument
if [ "$#" -eq 0 ]; then
  echo "ERROR: No options provided."
  usage
  return 1
fi

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --iftest) iftest="$2"; shift; shift ;;
    --qualada) qualada="$2"; shift; shift ;;
    --qualfil) qualfil="$2"; shift; shift ;;
    --errbar) errbar="$2"; shift; shift ;;
    --errfil) errfil="$2"; shift; shift ;;
    --thr) thr="$2"; shift; shift ;;
    --lenada) lenada="$2"; shift; shift ;;
    --lenfil) lenfil="$2"; shift; shift ;;
    --readsP5) readsP5="$2"; shift; shift ;;
    --readsP7) readsP7="$2"; shift; shift ;;
    --barcodes) barcodes="$2"; shift; shift ;;
    --directory) directory="$2"; shift; shift ;;
    --A_p5_5p) A_p5_5p="$2"; shift; shift ;;
    --A_p5_3p) A_p5_3p="$2"; shift; shift ;;
    --A_p7_5p) A_p7_5p="$2"; shift; shift ;;
    --A_p7_3p) A_p7_3p="$2"; shift; shift ;;
    --cutsites) cutsites="$2"; shift; shift;;
    --help) usage; return 0 ;;
    *) echo "Unknown option: $1"; usage; return 1 ;;
  esac
done

# Load cutsite sequence from file
if [ -z "$cutsites" ]; then
  echo "ERROR: --cutsites file is required."
  return 1
fi

if [ ! -f "$cutsites" ]; then
  echo "ERROR: Cut-site file '$cutsites' not found."
  return 1
fi

while IFS='=' read -r key value; do
  if [[ "$key" =~ ^C_ && -n "$value" ]]; then
    declare "$key=$value"
  fi
done < "$cutsites"


# Print all variables for verification
echo "Variables set:"
echo "iftest=$iftest"
echo "qualada=$qualada"
echo "qualfil=$qualfil"
echo "errbar=$errbar"
echo "errfil=$errfil"
echo "thr=$thr"
echo "lenada=$lenada"
echo "lenfil=$lenfil"
echo "readsP5=$readsP5"
echo "readsP7=$readsP7"
echo "barcodes=$barcodes"
echo "directory=$directory"
echo "A_p5_5p=$A_p5_5p"
echo "A_p5_3p=$A_p5_3p"
echo "A_p7_5p=$A_p7_5p"
echo "A_p7_3p=$A_p7_3p"
echo "C_p5_5p=$C_p5_5p"
echo "C_p5_3p=$C_p5_3p"
echo "C_p7_5p=$C_p7_5p"
echo "C_p7_3p=$C_p7_3p"

# Start

echo "# reading barcodes"

cat $barcodes | grep '>' | sed 's/>//' > barcodes.list

echo "# counting input reads"

gzip -cd $readsP5 | wc -l | awk '{print $1/4}' > Start.1.count
gzip -cd $readsP7 | wc -l | awk '{print $1/4}' > Start.2.count
echo "Original" > first.column
paste first.column Start.1.count Start.2.count > Start.counts

echo "# 1. removing sequencing adapters:"

cutadapt -j "$thr" -g "$A_p5_5p" -a "$A_p5_3p" -G "$A_p7_5p" -A "$A_p7_3p" -n 2 -m "$lenada" --max-n 0 -q "$qualada" -o "$directory"Adapters.1.fq.gz -p "$directory"Adapters.2.fq.gz $readsP5 $readsP7 > Adapters.txt

gzip -cd "$directory"Adapters.1.fq.gz | wc -l | awk '{print $1/4}' > Adapters.1.count
gzip -cd "$directory"Adapters.2.fq.gz | wc -l | awk '{print $1/4}' > Adapters.2.count
echo "After_adapters_removal" > first.column
paste first.column Adapters.1.count Adapters.2.count > Adapters.counts

echo "# 2. demultiplexing:"

cutadapt --no-indels -j "$thr" -e "$errbar" -g ^file:"$barcodes" --action=trim --untrimmed-output "$directory"Barcodes.1.untrimmed.fq.gz --untrimmed-paired-output "$directory"Barcodes.2.untrimmed.fq.gz -o "$directory"Barcodes.1.{name}.fq.gz -p "$directory"Barcodes.2.{name}.fq.gz "$directory"Adapters.1.fq.gz "$directory"Adapters.2.fq.gz > Barcodes.txt

cat barcodes.list | head -n "$iftest" | while read barcode; do # normal lociation of start of the loop

gzip -cd "$directory"Barcodes.1."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Barcodes.1."$barcode".count
gzip -cd "$directory"Barcodes.2."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Barcodes.2."$barcode".count
echo "$barcode"."demux" > first.column
paste first.column Barcodes.1."$barcode".count Barcodes.2."$barcode".count > Barcodes."$barcode".counts

echo "# 3. filteringcutsites:"

echo "$barcode"

cutadapt --no-indels -j "$thr" -e "$errfil" -a "$C_p5_5p...$C_p5_3p" -A "$C_p7_5p...$C_p7_3p" -m "$lenfil" --max-n 0 -q "$qualfil" --action=retain --untrimmed-output "$directory"Filtered.1."$barcode".untrimmed.fq.gz --untrimmed-paired-output "$directory"Filtered.2."$barcode".untrimmed.fq.gz -o "$directory"Filtered.1."$barcode".fq.gz -p "$directory"Filtered.2."$barcode".fq.gz "$directory"Barcodes.1."$barcode".fq.gz "$directory"Barcodes.2."$barcode".fq.gz > Filtered."$barcode".txt

gzip -cd "$directory"Filtered.1."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Filtered.1."$barcode".count
gzip -cd "$directory"Filtered.2."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Filtered.2."$barcode".count
echo "$barcode"."filtered" > first.column
paste first.column Filtered.1."$barcode".count Filtered.2."$barcode".count > Filtered."$barcode".counts

echo "# 4. rescuing untrimmed reads:"

# P5 forward

cutadapt --no-indels -j "$thr" -e "$errfil" -g "$C_p5_5p" -m "$lenfil" --max-n 0 -q "$qualfil" --action=retain --discard-untrimmed -o "$directory"Rescued.1."$barcode".fq.gz "$directory"Filtered.1."$barcode".untrimmed.fq.gz > Rescued.1."$barcode".txt

gzip -cd "$directory"Rescued.1."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Rescued.1."$barcode".count

# P7 reverse

cutadapt --no-indels -j "$thr" -e "$errfil" -g "$C_p7_5p" -m "$lenfil" --max-n 0 -q "$qualfil" --action=retain --discard-untrimmed -o "$directory"Rescued.2."$barcode".fq.gz "$directory"Filtered.2."$barcode".untrimmed.fq.gz > Rescued.2."$barcode".txt

gzip -cd "$directory"Rescued.2."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Rescued.2."$barcode".count

echo "$barcode"."rescued" > first.column
paste first.column Rescued.1."$barcode".count Rescued.2."$barcode".count > Rescued."$barcode".counts

done

cat Start.counts Adapters.counts Barcodes.*.counts Filtered.*.counts Rescued.*.counts | sed 's/ /\t/g'  > Counts.stat

rm *.count

rm first.column

echo "Done!"
