#!/bin/bash

directory="/mnt/qnap/projects/RafalWoycicki/" # "" if locally

cat barcodes.list | while read barcode; do

echo "$barcode"

cat "$directory"Short.Filtered."$barcode".p5p7dedupl.1.fq "$directory"Short.Rescued.2."$barcode".p7dedupl.fq > "$directory"Short.Fil1Res2."$barcode".p5p7dedupl.fq

mkdir "$directory""$barcode".ustacks.res.DIR

ustacks -f "$directory"Short.Fil1Res2."$barcode".p5p7dedupl.fq -o "$directory""$barcode".ustacks.res.DIR -M 2 -m 3 -t 10 2>&1 | tee "$barcode".ustacks.res.out

done

