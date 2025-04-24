#!/bin/bash

# Must be sourced: source ddradseq_dedup.bash [options]

conda activate cutadapt

iftest="1" # 1 if test "-0" if not test
directory="" # set the directory to save the files, set "" - if localy
p5len=130 #length of p5 read with SBF1 cut site beginning with ^TGCA
p7len=140 #length of p7 read with MSE1 cut site and DBR.
len=130 # length of the final reads
thr=10 # nymber of processor therds to use for cutadapt
err=0 # numner of mismatches allowed adapter sequence for cutadapt
p7_5p_seq="P7read5primDBR_MSE1=^NNNNNNNNGC" # name and sequence of the nucleotides filtered by the final reads shortening cutadapt script (for details check the cutadapt manual)

# Help message
usage() {
  echo "Usage: source ddradseq_dedup.bash [options]"
  echo "Options:"
  echo "  --iftest <value>         Set test mode: 1 for testing one barcode (default), -0 for not testing"
  echo "  --directory <path>       Project output directory set to blank/nothing if locally (default: locally)"
  echo "  --p5len <value>          length of p5 read with SBF1 cut site beginning with ^TGCA (default: 130)"
  echo "  --p7len <value>          length of p7 read with MSE1 cut site and DBR (default: 140)"
  echo "  --len <value>            length of the final reads (default: 130)"
  echo "  --thr <value>            Number of threads (default: 10)"
  echo "  --err <value>            Number of mismatches allowed in adapter sequence for cutadapt (default: 0)"
  echo "  --p7_5p_seq <value>      P7 5' outer cutsite sequence (default: P7read5primDBR_MSE1=^NNNNNNNNGC)"
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
    --directory) directory="$2"; shift; shift ;;
    --p5len) p5len="$2"; shift; shift ;;
    --p7len) p7len="$2"; shift; shift ;;
    --len) len="$2"; shift; shift ;;
    --thr) thr="$2"; shift; shift ;;
    --err) err="$2"; shift; shift ;;
    --p7_5p_seq) p7_5p_seq="$2"; shift; shift ;;
    --help) usage; return 0 ;;
    *) echo "Unknown option: $1"; usage; return 1 ;;
  esac
done

# Print all variables for verification
echo "Variables set:"
echo "iftest=$iftest"
echo "directory=$directory"
echo "p5len=$p5len"
echo "p7len=$p7len"
echo "len=$len"
echo "thr=$thr"
echo "err=$err"
echo "p7_5p_seq=$p7_5p_seq"

# Start

cat barcodes.list | head -n "$iftest" | while read barcode; do # normal lociation of start of the loop

# deduplicating in Filtered in both reads

echo "$barcode"

echo "1"
seqtk seq -A "$directory"Filtered.1."$barcode".fq.gz > "$directory"Filtered.1."$barcode".fasta
seqtk seq -A "$directory"Filtered.2."$barcode".fq.gz > "$directory"Filtered.2."$barcode".fasta

echo "2"
cat "$directory"Filtered.1."$barcode".fasta | awk '/^>/{if (seq) print header "\t" seq; split($1,h," "); header=substr(h[1],2); seq=""} !/^>/ {seq = seq $0} END {if (seq) print header "\t" seq}' > "$directory"Filtered.1."$barcode".table
cat "$directory"Filtered.2."$barcode".fasta | awk '/^>/{if (seq) print header "\t" seq; split($1,h," "); header=substr(h[1],2); seq=""} !/^>/ {seq = seq $0} END {if (seq) print header "\t" seq}' > "$directory"Filtered.2."$barcode".table

echo "3"
paste "$directory"Filtered.1."$barcode".table "$directory"Filtered.2."$barcode".table | awk -v p5len="$p5len" -v p7len="$p7len" '$1==$3 && length($2)>=p5len && length($4)>=p7len { $2=substr($2,1,p5len); $4=substr($4,1,p7len); if (match($4, /[AGCT]{8}GCTAA/)) { sub(/([AGCT]{8})GCTAA/, substr($4, RSTART, RLENGTH-5)"GC\tTAA", $4); } print $1"\t"$2"\t"$4; }' > "$directory"Filtered."$barcode".all.table

echo "4"
cat "$directory"Filtered."$barcode".all.table | awk '{ print $1"\t"$3"_"$2 }' | sort -k2,2 > "$directory"Filtered."$barcode".p5sorted.table
cat "$directory"Filtered."$barcode".all.table | awk '{ print $1"\t"$3"_"$4 }' | sort -k2,2 > "$directory"Filtered."$barcode".p7sorted.table

echo "5"
cat "$directory"Filtered."$barcode".p5sorted.table | uniq -c -f1 | awk '{sub(/^ +/, "", $0); sub(/ +/, "\t", $0); print $0}' | awk '$1>0' | cut -f2 | sort > Filtered."$barcode".p5dedupl.list

echo "6"
awk 'NR==FNR {dedup[$1]; next} $1 in dedup' Filtered."$barcode".p5dedupl.list "$directory"Filtered."$barcode".p7sorted.table > "$directory"Filtered."$barcode".p7sorted.p5dedupl.table

echo "7"
cat "$directory"Filtered."$barcode".p7sorted.p5dedupl.table | sort -k2,2 | uniq -c -f1 | awk '{sub(/^ +/, "", $0); sub(/ +/, "\t", $0); print $0}' | awk '$1==1' | cut -f2 | sort > Filtered."$barcode".p7sorted.p5dedupl.list

echo "8"
comm -12 Filtered."$barcode".p7sorted.p5dedupl.list Filtered."$barcode".p5dedupl.list > Filtered."$barcode".p5p7dedupl.list

echo "9"
seqtk subseq "$directory"Filtered.1."$barcode".fq.gz Filtered."$barcode".p5p7dedupl.list > "$directory"Filtered."$barcode".p5p7dedupl.1.fq
seqtk subseq "$directory"Filtered.2."$barcode".fq.gz Filtered."$barcode".p5p7dedupl.list > "$directory"Filtered."$barcode".p5p7dedupl.2.fq

echo "10"
#shortening the reads so they are equal in length
cutadapt -j "$thr" --length "$len" -o "$directory"Short.Filtered."$barcode".p5p7dedupl.1.fq "$directory"Filtered."$barcode".p5p7dedupl.1.fq > FilterShort.1.txt
cutadapt -j "$thr" --length "$len" -e "$err" -g "$p7_5p_seq" -o "$directory"Short.Filtered."$barcode".p5p7dedupl.2.fq "$directory"Filtered."$barcode".p5p7dedupl.2.fq > FilterShort.2.txt #removing DBR with additional "GC" bases from P7 reads

echo "11"
cat "$directory"Short.Filtered."$barcode".p5p7dedupl.1.fq | wc -l | awk '{print $1/4}' > Short.Filtered."$barcode".p5p7dedupl.1.count
cat "$directory"Short.Filtered."$barcode".p5p7dedupl.2.fq | wc -l | awk '{print $1/4}' > Short.Filtered."$barcode".p5p7dedupl.2.count
echo "$barcode"."dedupl_filtered" > first.column
paste first.column Short.Filtered."$barcode".p5p7dedupl.1.count Short.Filtered."$barcode".p5p7dedupl.2.count > FilDedup."$barcode".p5p7dedupl.counts

echo "12"
# deduplicating in Rescued in P7 read only

seqtk seq -A "$directory"Rescued.2."$barcode".fq.gz > "$directory"Rescued.2."$barcode".fasta

cat "$directory"Rescued.2."$barcode".fasta | tr '\n' '\t' | sed 's/\t>/\n/g' | sed 's/>//' | sed 's/\t$/\n/' | sed 's/ /\t/' | cut -f1,3 > "$directory"Rescued.2."$barcode".table

cat "$directory"Rescued.2."$barcode".table | sed "s/\t/_fc_/" | sed "s/_fc_/\t_fc_\t/" | awk -v p7len="$p7len" 'length($3)>=p7len' | awk -v p7len="$p7len" -F"\t" -vOFS="\t" '{ $3=substr($3,1,p7len) };1' | sed "s/\t_fc_\t\([AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT]\)GCTAA/\t\1GC\tTAA/" > "$directory"Rescued.2."$barcode".all.table

cat "$directory"Rescued.2."$barcode".all.table | awk '{ print $1"\t"$2"_"$3 }' | sort -k2,2 > "$directory"Rescued.2."$barcode".p7sorted.table

cat "$directory"Rescued.2."$barcode".p7sorted.table | uniq -c -f1 | sed 's/^    //g' | sed 's/^  //g' | sed 's/^ //g' | sed 's/ /\t/' | awk '$1>0' | cut -f2 | sort > Rescued.2."$barcode".p7dedupl.list

seqtk subseq "$directory"Rescued.2."$barcode".fq.gz Rescued.2."$barcode".p7dedupl.list > "$directory"Rescued.2."$barcode".p7dedupl.fq

cutadapt -j "$thr" --length "$len" -e "$err" -g "$p7_5p_seq" -o "$directory"Short.Rescued.2."$barcode".p7dedupl.fq "$directory"Rescued.2."$barcode".p7dedupl.fq > ResShort.2.txt

cat "$directory"Short.Rescued.2."$barcode".p7dedupl.fq | wc -l | awk '{print $1/4}' > Short.Rescued.2."$barcode".p7dedupl.count
echo "$barcode"."dedupl_rescued2" > first.column
paste first.column Short.Rescued.2."$barcode".p7dedupl.count > ResDedup.2."$barcode".p7dedupl.counts

done

cat Start.counts Adapters.counts Barcodes.*.counts Filtered.*.counts Rescued.*.counts FilDedup.*.counts ResDedup.*.counts | sed 's/ /\t/g'  > Counts_dedupl.stat

rm *.count

rm -rf "$directory"*.fasta
rm -rf "$directory"*.table

rm first.column
