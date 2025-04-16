#!/bin/bash

conda activate cutadapt

iftest="-0" # 1 if test "-0" if not test
directory="/mnt/qnap/projects/RafalWoycicki/" # "" - if localy
p5len=130 #length of p5 read with SBF1 cut site beginning with ^TGCA
p7len=140 #length of p7 read with MSE1 cut site and DBR.
len=130 # lengthof the final reads

cat barcodes.list | head -n "$iftest" | while read barcode; do # normal lociation of start of the loop

# deduplicating in Filtered in both reads

echo "$barcode"

seqtk seq -A "$directory"Filtered.1."$barcode".fq.gz > "$directory"Filtered.1."$barcode".fasta
seqtk seq -A "$directory"Filtered.2."$barcode".fq.gz > "$directory"Filtered.2."$barcode".fasta

cat "$directory"Filtered.1."$barcode".fasta | tr '\n' '\t' | sed 's/\t>/\n/g' | sed 's/>//' | sed 's/\t$/\n/' | sed 's/ /\t/' | cut -f1,3 > "$directory"Filtered.1."$barcode".table
cat "$directory"Filtered.2."$barcode".fasta | tr '\n' '\t' | sed 's/\t>/\n/g' | sed 's/>//' | sed 's/\t$/\n/' | sed 's/ /\t/' | cut -f1,3 > "$directory"Filtered.2."$barcode".table

paste "$directory"Filtered.1."$barcode".table "$directory"Filtered.2."$barcode".table | awk '$1==$3' | cut -f1,2,4 | sed "s/\t/_fc_/" | sed "s/\t/_sc_/" | sed "s/_fc_/\t_fc_\t/" | sed "s/_sc_/\t_sc_\t/" | awk -v p5len="$p5len" -v p7len="$p7len" 'length($3)>=p5len && length($5)>=p7len' | awk -v p5len="$p5len" -v p7len="$p7len" -F"\t" -vOFS="\t" '{ ($3=substr($3,1,p5len)) ($5=substr($5,1,p7len)) };1' | sed "s/\t_sc_\t\([AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT]\)GCTAA/\t\1GC\tTAA/" | sed "s/\t_fc_\t/\t/" > "$directory"Filtered."$barcode".all.table

cat "$directory"Filtered."$barcode".all.table | awk '{ print $1"\t"$3"_"$2 }' | sort -k2,2 > "$directory"Filtered."$barcode".p5sorted.table
cat "$directory"Filtered."$barcode".all.table | awk '{ print $1"\t"$3"_"$4 }' | sort -k2,2 > "$directory"Filtered."$barcode".p7sorted.table

cat "$directory"Filtered."$barcode".p5sorted.table | uniq -c -f1 | sed 's/^    //g' | sed 's/^  //g' | sed 's/^ //g' | sed 's/ /\t/' | awk '$1>0' | cut -f2 | sort > Filtered."$barcode".p5dedupl.list

grep -w -f Filtered."$barcode".p5dedupl.list "$directory"Filtered."$barcode".p7sorted.table > "$directory"Filtered."$barcode".p7sorted.p5dedupl.table

cat "$directory"Filtered."$barcode".p7sorted.p5dedupl.table | sort -k2,2 | uniq -c -f1 | sed 's/^    //g' | sed 's/^  //g' | sed 's/^ //g' | sed 's/ /\t/' | awk '$1==1' | cut -f2 | sort > Filtered."$barcode".p7sorted.p5dedupl.list

comm -12 Filtered."$barcode".p7sorted.p5dedupl.list Filtered."$barcode".p5dedupl.list > Filtered."$barcode".p5p7dedupl.list

seqtk subseq "$directory"Filtered.1."$barcode".fq.gz Filtered."$barcode".p5p7dedupl.list > "$directory"Filtered."$barcode".p5p7dedupl.1.fq
seqtk subseq "$directory"Filtered.2."$barcode".fq.gz Filtered."$barcode".p5p7dedupl.list > "$directory"Filtered."$barcode".p5p7dedupl.2.fq

#shortening the reads so they are equal in length

cutadapt -j 10 --length "$len" -o "$directory"Short.Filtered."$barcode".p5p7dedupl.1.fq "$directory"Filtered."$barcode".p5p7dedupl.1.fq > FilterShort.1.txt
cutadapt -j 10 --length "$len" -e 0 -g P7read5primDBR_MSE1=^NNNNNNNNGC -o "$directory"Short.Filtered."$barcode".p5p7dedupl.2.fq "$directory"Filtered."$barcode".p5p7dedupl.2.fq > FilterShort.2.txt #removing DBR with additional "GC" bases from P7 reads

cat "$directory"Short.Filtered."$barcode".p5p7dedupl.1.fq | wc -l | awk '{print $1/4}' > Short.Filtered."$barcode".p5p7dedupl.1.count
cat "$directory"Short.Filtered."$barcode".p5p7dedupl.2.fq | wc -l | awk '{print $1/4}' > Short.Filtered."$barcode".p5p7dedupl.2.count
echo "$barcode"."dedupl_filtered" > first.column
paste first.column Short.Filtered."$barcode".p5p7dedupl.1.count Short.Filtered."$barcode".p5p7dedupl.2.count > FilDedup."$barcode".p5p7dedupl.counts

# deduplicating in Rescued in P7 read only

seqtk seq -A "$directory"Rescued.2."$barcode".fq.gz > "$directory"Rescued.2."$barcode".fasta

cat "$directory"Rescued.2."$barcode".fasta | tr '\n' '\t' | sed 's/\t>/\n/g' | sed 's/>//' | sed 's/\t$/\n/' | sed 's/ /\t/' | cut -f1,3 > "$directory"Rescued.2."$barcode".table

cat "$directory"Rescued.2."$barcode".table | sed "s/\t/_fc_/" | sed "s/_fc_/\t_fc_\t/" | awk -v p7len="$p7len" 'length($3)>=p7len' | awk -v p7len="$p7len" -F"\t" -vOFS="\t" '{ $3=substr($3,1,p7len) };1' | sed "s/\t_fc_\t\([AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT]\)GCTAA/\t\1GC\tTAA/" > "$directory"Rescued.2."$barcode".all.table

cat "$directory"Rescued.2."$barcode".all.table | awk '{ print $1"\t"$2"_"$3 }' | sort -k2,2 > "$directory"Rescued.2."$barcode".p7sorted.table

cat "$directory"Rescued.2."$barcode".p7sorted.table | uniq -c -f1 | sed 's/^    //g' | sed 's/^  //g' | sed 's/^ //g' | sed 's/ /\t/' | awk '$1>0' | cut -f2 | sort > Rescued.2."$barcode".p7dedupl.list

seqtk subseq "$directory"Rescued.2."$barcode".fq.gz Rescued.2."$barcode".p7dedupl.list > "$directory"Rescued.2."$barcode".p7dedupl.fq

cutadapt -j 10 --length "$len" -e 0 -g P7read5primDBR_MSE1=^NNNNNNNNGC -o "$directory"Short.Rescued.2."$barcode".p7dedupl.fq "$directory"Rescued.2."$barcode".p7dedupl.fq > ResShort.2.txt

cat "$directory"Short.Rescued.2."$barcode".p7dedupl.fq | wc -l | awk '{print $1/4}' > Short.Rescued.2."$barcode".p7dedupl.count
echo "$barcode"."dedupl_rescued2" > first.column
paste first.column Short.Rescued.2."$barcode".p7dedupl.count > ResDedup.2."$barcode".p7dedupl.counts

done

cat Start.counts Adapters.counts Barcodes.*.counts Filtered.*.counts Rescued.*.counts FilDedup.*.counts ResDedup.*.counts | sed 's/ /\t/g'  > Counts_dedupl.stat

rm *.count

rm first.column


