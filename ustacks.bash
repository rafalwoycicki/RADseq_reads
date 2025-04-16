#!/bin/bash

directory="/mnt/qnap/projects/RafalWoycicki/MP/P1_MKDL240011297-1A_22LMFJLT4_L2/" # "" if locally

cat barcodes.list | while read barcode; do

mkdir "$directory""$barcode".ustacks.DIR

ustacks -f "$directory"Short.Filtered."$barcode".p5p7dedupl.1.fq -o "$directory""$barcode".ustacks.DIR -M 2 -m 3 -t 10 2>&1 | tee "$barcode".ustacks.out

done


