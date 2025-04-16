#!/bin/bash

rm -rf ustacks.res.stat

cat barcodes.list | while read barcode; do

echo "$barcode"

cat "$barcode".ustacks.res.out | grep "Final" | tr "\n" "\t" | sed 's/.*built: //' | sed 's/\t.*mean=/\t/' | sed 's/;.*//' | sed 's/$/\n/' | tr '.' ',' > stat.tmp

echo "$barcode" > barcode.tmp

paste barcode.tmp stat.tmp >> ustacks.res.stat

done

