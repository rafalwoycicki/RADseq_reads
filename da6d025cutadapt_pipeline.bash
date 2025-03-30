#!/bin/bash

### !!! need to run script using: source <<file name>>
conda activate cutadapt

# pipeline for trimming and demultiplexing illumina reads after RADseq
# using cutadapt 5.0
qualada="20,20" # quality trimming adapters
qualfil="20,20" # quality trimming filtering
errbar=1 # number of allowed errors for demultiplexing
errfil=0 # number of allowed errors for CUT site filtering
errfilanc=1 # number of allowed error for anchored CUT site filtering
thr=6 # number of processor threads to use
lenada=130 # minimum length for adapter removal
lenfil=100 # minimum length for filtering cut sites
reads="/mnt/qnap/projects/RADseq/MPabijan_RAD/X204SC24124326-Z01-F001/01.RawData/P1/P1_MKDL240011297-1A_22LMFJLT4_L2_" # template for reads 
readsP5="$reads""1.fq.gz" # reads from P5 primer - forward
readsP7="$reads""2.fq.gz" # reads from P7 primer - reverse

barcodes=barcodes_P1_cutadapt.txt # fasta file with barcodes
#barcodesplus=barcodes_P1_cutadapt_plusCUT.txt #fasta file with barcodes and P5 cut site

echo "$qualada $qualfil $errbar $errfil $errfilanc $thr $lenada $lenfil $reads $readsP5 $readsP7 $barcodes"

cat $barcodes | grep '>' | sed 's/>//' > barcodes.list

zcat $readsP5 | wc -l | awk '{print $1/4}' > Start.1.count
zcat $readsP7 | wc -l | awk '{print $1/4}' > Start.2.count
echo "Original" > first.column
paste first.column Start.1.count Start.2.count > Start.counts

# 1. removing sequencing adapters:

cutadapt -j "$thr" -g P5read5prim=TACACGACGCTCTTCCGATCT -a P5read3prim=AGATCGGAAGAGCACACGTCT -G P7read5prim=AGACGTGTGCTCTTCCGATCT -A P7read3prim=AGATCGGAAGAGCGTCGTGTA -n 2 -m "$lenada" --max-n 0 -q "$qualada" -o Adapters.1.fq.gz -p Adapters.2.fq.gz $readsP5 $readsP7 > Adapters.txt

zcat Adapters.1.fq.gz | wc -l | awk '{print $1/4}' > Adapters.1.count
zcat Adapters.2.fq.gz | wc -l | awk '{print $1/4}' > Adapters.2.count
echo "After_adapters_removal" > first.column
paste first.column Adapters.1.count Adapters.2.count > Adapters.counts

# 2. demultiplexing

cutadapt --no-indels -j "$thr" -e "$errbar" -g ^file:"$barcodes" --action=trim --untrimmed-output Barcodes.1.untrimmed.fq.gz --untrimmed-paired-output Barcodes.2.untrimmed.fq.gz -o Barcodes.1.{name}.fq.gz -p Barcodes.2.{name}.fq.gz Adapters.1.fq.gz Adapters.2.fq.gz > Barcodes.txt

cat barcodes.list | while read barcode; do

zcat Barcodes.1."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Barcodes.1."$barcode".count
zcat Barcodes.2."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Barcodes.2."$barcode".count
echo "$barcode"."demux" > first.column
paste first.column Barcodes.1."$barcode".count Barcodes.2."$barcode".count > Barcodes."$barcode".counts

# 3. filteringcutsites:

cutadapt --no-indels -j "$thr" -e "$errfil" -a "P5read5primSBF1=^TGCAGG;e=$errfilanc...P5read3primMSE1_DBR=TTAGCNNNNNNNN" -A "P7read5primDBR_MSE1=^NNNNNNNNGCTAA;e=$errfilanc...P7read3primSBF1=CCTGCA" -m "$lenfil" --max-n 0 -q "$qualfil" --action=retain --untrimmed-output Filtered.1."$barcode".untrimmed.fq.gz --untrimmed-paired-output Filtered.2."$barcode".untrimmed.fq.gz -o Filtered.1."$barcode".fq.gz -p Filtered.2."$barcode".fq.gz Barcodes.1."$barcode".fq.gz Barcodes.2."$barcode".fq.gz > Filtered."$barcode".txt

zcat Filtered.1."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Filtered.1."$barcode".count
zcat Filtered.2."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Filtered.2."$barcode".count
echo "$barcode"."filtered" > first.column
paste first.column Filtered.1."$barcode".count Filtered.2."$barcode".count > Filtered."$barcode".counts

# 4. rescuing untrimmed read:

# P5 forward

cutadapt --no-indels -j "$thr" -e "$errfil" -g "P5read5primSBF1=^TGCAGG;e=$errfilanc" -m "$lenfil" --max-n 0 -q "$qualfil" --action=retain --discard-untrimmed -o Rescued.1."$barcode".fq.gz Filtered.1."$barcode".untrimmed.fq.gz > Rescued.1."$barcode".txt

zcat Rescued.1."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Rescued.1."$barcode".count

# P7 reverse

cutadapt --no-indels -j "$thr" -e "$errfil" -g "P7read5primDBR_MSE1=^NNNNNNNNGCTAA;e=$errfilanc" -m "$lenfil" --max-n 0 -q "$qualfil" --action=retain --discard-untrimmed -o Rescued.2."$barcode".fq.gz Filtered.2."$barcode".untrimmed.fq.gz > Rescued.2."$barcode".txt

zcat Rescued.2."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Rescued.2."$barcode".count

echo "$barcode"."rescued" > first.column
paste first.column Rescued.1."$barcode".count Rescued.2."$barcode".count > Rescued."$barcode".counts

done

cat Start.counts Adapters.counts Barcodes.*.counts Filtered.*.counts Rescued.*.counts | sed 's/ /\t/g'  > Counts.stat

rm *.count

rm first.column


