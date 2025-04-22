#!/bin/bash

### !!! need to run script using: source <<file name>>
conda activate cutadapt

# pipeline for trimming and demultiplexing illumina reads after RADseq
# using cutadapt 5.0
iftest="-0" # 1 for testing one barcode, "-0" for not testing
qualada="20,20" # quality trimming adapters
qualfil="20,20" # quality trimming filtering
errbar=1 # number of allowed errors for demultiplexing
errfil=0 # number of allowed errors for CUT site filtering
thr=40 # number of processor threads to use
lenada=140 # minimum length for adapter removal
lenfil=140 # minimum length for filtering cut sites
readsP5="1.fq.gz" # reads from P5 primer - forward
readsP7="2.fq.gz" # reads from P7 primer - reverse
barcodes="barcodes_P1_cutadapt.txt" # fasta file with barcodes
directory="/mnt/qnap/projects/RafalWoycicki/" # "" - if localy
#adapters
A_p5_5p="P5read5prim=TACACGACGCTCTTCCGATCT" # read P5 5prim sequencing adapter sequence
A_p5_3p="P5read3prim=AGATCGGAAGAGCACACGTCT" # read P5 3prim sequencing adapter sequence
A_p7_5p="P7read5prim=AGACGTGTGCTCTTCCGATCT" # read P7 5prim sequencing adapter sequence
A_p7_3p="P7read3prim=AGATCGGAAGAGCGTCGTGTA" # read P7 3prim sequencing adapter sequence
#cutsites
C_p5_5p="P5read5primSBF1=^TGCAGG" # read P5 5prim cut site for the SBF1 RE
C_p5_3p="P5read3primMSE1_DBR=TTAGCNNNNNNNN" # read P5 3prim cut site for the MSE1 RE including DBR region
C_p7_5p="P7read5primDBR_MSE1=^NNNNNNNNGCTAA" # read P7 5prim cut site for the MSE1 RE including DBR region
C_p7_3p="P7read3primSBF1=CCTGCA" # read P7 3 prim cut site for the SBF1 RE

echo "variable: $iftest $qualada $qualfil $errbar $errfil $errfilanc $thr $lenada $lenfil $readsP5 $readsP7 $barcodes $directory $A_p5_5p $A_p5_3p $A_p7_5p $A_p7_3p $C_p5_5p $C_p5_3p $C_p7_5p $C_p7_3p"

cat $barcodes | grep '>' | sed 's/>//' > barcodes.list

zcat $readsP5 | wc -l | awk '{print $1/4}' > Start.1.count
zcat $readsP7 | wc -l | awk '{print $1/4}' > Start.2.count
echo "Original" > first.column
paste first.column Start.1.count Start.2.count > Start.counts

echo "# 1. removing sequencing adapters:"

cutadapt -j "$thr" -g "$A_p5_5p" -a "$A_p5_3p" -G "$A_p7_5p" -A "$A_p7_3p" -n 2 -m "$lenada" --max-n 0 -q "$qualada" -o "$directory"Adapters.1.fq.gz -p "$directory"Adapters.2.fq.gz $readsP5 $readsP7 > Adapters.txt

zcat "$directory"Adapters.1.fq.gz | wc -l | awk '{print $1/4}' > Adapters.1.count
zcat "$directory"Adapters.2.fq.gz | wc -l | awk '{print $1/4}' > Adapters.2.count
echo "After_adapters_removal" > first.column
paste first.column Adapters.1.count Adapters.2.count > Adapters.counts

echo "# 2. demultiplexing:"

cutadapt --no-indels -j "$thr" -e "$errbar" -g ^file:"$barcodes" --action=trim --untrimmed-output "$directory"Barcodes.1.untrimmed.fq.gz --untrimmed-paired-output "$directory"Barcodes.2.untrimmed.fq.gz -o "$directory"Barcodes.1.{name}.fq.gz -p "$directory"Barcodes.2.{name}.fq.gz "$directory"Adapters.1.fq.gz "$directory"Adapters.2.fq.gz > Barcodes.txt

cat barcodes.list | head -n "$iftest" | while read barcode; do # normal lociation of start of the loop

zcat "$directory"Barcodes.1."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Barcodes.1."$barcode".count
zcat "$directory"Barcodes.2."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Barcodes.2."$barcode".count
echo "$barcode"."demux" > first.column
paste first.column Barcodes.1."$barcode".count Barcodes.2."$barcode".count > Barcodes."$barcode".counts

echo "# 3. filteringcutsites:"

echo "$barcode"

cutadapt --no-indels -j "$thr" -e "$errfil" -a "$C_p5_5p...$C_p5_3p" -A "$C_p7_5p...$C_p7_3p" -m "$lenfil" --max-n 0 -q "$qualfil" --action=retain --untrimmed-output "$directory"Filtered.1."$barcode".untrimmed.fq.gz --untrimmed-paired-output "$directory"Filtered.2."$barcode".untrimmed.fq.gz -o "$directory"Filtered.1."$barcode".fq.gz -p "$directory"Filtered.2."$barcode".fq.gz "$directory"Barcodes.1."$barcode".fq.gz "$directory"Barcodes.2."$barcode".fq.gz > Filtered."$barcode".txt

zcat "$directory"Filtered.1."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Filtered.1."$barcode".count
zcat "$directory"Filtered.2."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Filtered.2."$barcode".count
echo "$barcode"."filtered" > first.column
paste first.column Filtered.1."$barcode".count Filtered.2."$barcode".count > Filtered."$barcode".counts

echo "# 4. rescuing untrimmed reads:"

# P5 forward

cutadapt --no-indels -j "$thr" -e "$errfil" -g "$C_p5_5p" -m "$lenfil" --max-n 0 -q "$qualfil" --action=retain --discard-untrimmed -o "$directory"Rescued.1."$barcode".fq.gz "$directory"Filtered.1."$barcode".untrimmed.fq.gz > Rescued.1."$barcode".txt

zcat "$directory"Rescued.1."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Rescued.1."$barcode".count

# P7 reverse

cutadapt --no-indels -j "$thr" -e "$errfil" -g "$C_p7_5p" -m "$lenfil" --max-n 0 -q "$qualfil" --action=retain --discard-untrimmed -o "$directory"Rescued.2."$barcode".fq.gz "$directory"Filtered.2."$barcode".untrimmed.fq.gz > Rescued.2."$barcode".txt

zcat "$directory"Rescued.2."$barcode".fq.gz | wc -l | awk '{print $1/4}' > Rescued.2."$barcode".count

echo "$barcode"."rescued" > first.column
paste first.column Rescued.1."$barcode".count Rescued.2."$barcode".count > Rescued."$barcode".counts

done

cat Start.counts Adapters.counts Barcodes.*.counts Filtered.*.counts Rescued.*.counts | sed 's/ /\t/g'  > Counts.stat

rm *.count

rm first.column

echo "Done!"
