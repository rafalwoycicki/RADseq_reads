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
barcodes="" # fasta file with barcodes
dbr_sequence="NNNNNNNN" # DBR sequence in nucleotides
dbr_pattern="[AGCT]{8}"
motif_cut_adapter="GC" # the outer part of the adapter adjacent to the CUT site and DBR site from P7 adapter
motif_cut_rerest="TAA"

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
  echo "  --barcodes <file>        Barcodes file (fasta file with barcodes) (default not set)"
  echo "  --dbr_sequence <value>	P7 DBR sequence (default 8 Ns: NNNNNNNN)"
  echo "  --dbr_pattern <regex>      P7 DBR regex pattern before motif (default: [AGCT]{8})"
  echo "  --motif_cut_adapter <value>  The outer part of the adapter adjacent to the CUT site and DBR site from P7 adapter (default: GC)"
  echo "  --motif_cut_rerest <value>   The inner part of the CUT site from P7 adapter site which stays with insert (default: TAA)"
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
    --barcodes) barcodes="$2"; shift; shift ;;
    --dbr_sequence) dbr_sequence="$2"; shift; shift ;;
    --dbr_pattern) dbr_pattern="$2"; shift; shift;;
    --motif_cut_adapter) motif_cut_adapter="$2"; shift; shift;;
    --motif_cut_rerest) motif_cut_rerest="$2"; shift; shift;;
    --help) usage; return 0 ;;
    *) echo "Unknown option: $1"; usage; return 1 ;;
  esac
done

# mofitf_suffix
motif_suffix="${motif_cut_adapter}${motif_cut_rerest}"
# p7_5p_seq
p7_5p_seq="^""${dbr_sequence}${motif_cut_adapter}"

# Print all variables for verification
echo "Variables set:"
echo "iftest=$iftest"
echo "directory=$directory"
echo "p5len=$p5len"
echo "p7len=$p7len"
echo "len=$len"
echo "thr=$thr"
echo "err=$err"
echo "barcodes=$barcodes"
echo "dbr_sequence=$dbr_sequence"
echo "dbr_pattern=$dbr_pattern"
echo "p7_5p_seq=$p7_5p_seq (auto-generated)"
echo "motif_suffix=$motif_suffix (auto-generated)"
echo "motif_cut_adapter=$motif_cut_adapter"
echo "motif_cut_rerest=$motif_cut_rerest"

# Start

cat $barcodes | grep '>' | sed 's/>//' > barcodes.list

cat barcodes.list | head -n "$iftest" | while read barcode; do # normal lociation of start of the loop

echo "$barcode"

echo "# deduplicating paired reads"

gzip -cd "$directory"Filtered.1."$barcode".fq.gz | awk 'NF && NR%4==1{sub(/^@/,">"); print} NF && NR%4==2{print}' | awk '/^>/{if (seq) print header "\t" seq; split($1,h," "); header=substr(h[1],2); seq=""} !/^>/ {seq = seq $0} END {if (seq) print header "\t" seq}' > "$directory"Filtered.1."$barcode".table
gzip -cd "$directory"Filtered.2."$barcode".fq.gz | awk 'NF && NR%4==1{sub(/^@/,">"); print} NF && NR%4==2{print}' | awk '/^>/{if (seq) print header "\t" seq; split($1,h," "); header=substr(h[1],2); seq=""} !/^>/ {seq = seq $0} END {if (seq) print header "\t" seq}' > "$directory"Filtered.2."$barcode".table

paste "$directory"Filtered.1."$barcode".table "$directory"Filtered.2."$barcode".table | awk -v p5len="$p5len" -v p7len="$p7len" -v dbr="$dbr_pattern" -v suffix="$motif_suffix" -v adapter="$motif_cut_adapter" -v rerest="$motif_cut_rerest" 'BEGIN{OFS="\t"; pattern="(" dbr ")" suffix} $1==$3 && length($2)>=p5len && length($4)>=p7len { $2=substr($2,1,p5len); $4=substr($4,1,p7len); if (match($4, pattern)) { sub(pattern, substr($4, RSTART, RLENGTH - length(suffix)) adapter "\t" rerest, $4) } print $1, $2, $4 }' > "$directory"Filtered."$barcode".all.table

cat "$directory"Filtered."$barcode".all.table | awk '{ print $1"\t"$3"_"$2 }' | sort -k2,2 > "$directory"Filtered."$barcode".p5sorted.table
cat "$directory"Filtered."$barcode".all.table | awk '{ print $1"\t"$3"_"$4 }' | sort -k2,2 > "$directory"Filtered."$barcode".p7sorted.table

cat "$directory"Filtered."$barcode".p5sorted.table | uniq -c -f1 | awk '{sub(/^ +/, "", $0); sub(/ +/, "\t", $0); print $0}' | awk '$1>0' | cut -f2 | sort > Filtered."$barcode".p5dedupl.list

awk 'NR==FNR {dedup[$1]; next} $1 in dedup' Filtered."$barcode".p5dedupl.list "$directory"Filtered."$barcode".p7sorted.table > "$directory"Filtered."$barcode".p7sorted.p5dedupl.table

cat "$directory"Filtered."$barcode".p7sorted.p5dedupl.table | sort -k2,2 | uniq -c -f1 | awk '{sub(/^ +/, "", $0); sub(/ +/, "\t", $0); print $0}' | awk '$1==1' | cut -f2 | sort > Filtered."$barcode".p7sorted.p5dedupl.list

comm -12 Filtered."$barcode".p7sorted.p5dedupl.list Filtered."$barcode".p5dedupl.list > Filtered."$barcode".p5p7dedupl.list

seqtk subseq "$directory"Filtered.1."$barcode".fq.gz Filtered."$barcode".p5p7dedupl.list > "$directory"Filtered."$barcode".p5p7dedupl.1.fq
seqtk subseq "$directory"Filtered.2."$barcode".fq.gz Filtered."$barcode".p5p7dedupl.list > "$directory"Filtered."$barcode".p5p7dedupl.2.fq

cutadapt -j "$thr" --length "$len" -o "$directory"Short.Filtered."$barcode".p5p7dedupl.1.fq "$directory"Filtered."$barcode".p5p7dedupl.1.fq > FilterShort.1.txt
cutadapt -j "$thr" --length "$len" -e "$err" -g "$p7_5p_seq" -o "$directory"Short.Filtered."$barcode".p5p7dedupl.2.fq "$directory"Filtered."$barcode".p5p7dedupl.2.fq > FilterShort.2.txt #removing DBR with additional "GC" bases from P7 reads

cat "$directory"Short.Filtered."$barcode".p5p7dedupl.1.fq | wc -l | awk '{print $1/4}' > Short.Filtered."$barcode".p5p7dedupl.1.count
cat "$directory"Short.Filtered."$barcode".p5p7dedupl.2.fq | wc -l | awk '{print $1/4}' > Short.Filtered."$barcode".p5p7dedupl.2.count
echo "$barcode"."dedupl_filtered" > first.column
paste first.column Short.Filtered."$barcode".p5p7dedupl.1.count Short.Filtered."$barcode".p5p7dedupl.2.count > FilDedup."$barcode".p5p7dedupl.counts

echo "# deduplicating single Rescued P7 reads"

seqtk seq -A "$directory"Rescued.2."$barcode".fq.gz > "$directory"Rescued.2."$barcode".fasta

cat "$directory"Rescued.2."$barcode".fasta | tr '\n' '\t' | sed 's/\t>/\n/g' | sed 's/>//' | sed 's/\t$/\n/' | sed 's/ /\t/' | cut -f1,3 > "$directory"Rescued.2."$barcode".table

cat "$directory"Rescued.2."$barcode".table | awk -v p7len="$p7len" -v dbr="$dbr_pattern" -v suffix="$motif_suffix" -v adapter="$motif_cut_adapter" -v rerest="$motif_cut_rerest" 'BEGIN{OFS="\t"; regex="(" dbr ")" suffix} length($2) >= p7len { $2 = substr($2, 1, p7len); if (match($2, regex)) { dbr_part = substr($2, RSTART, RLENGTH - length(suffix)); $2 = rerest substr($2, RSTART + RLENGTH); print $1, dbr_part adapter, $2 } }' > "$directory"Rescued.2."$barcode".all.table

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

rm first.column
