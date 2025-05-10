#!/bin/bash

# Usage:
# ./demux_header.bash R1.fastq.gz R2.fastq.gz barcode_mapping.txt first|second [swap]

R1=$1
R2=$2
MAPPING_FILE=$3
MODE=$4      # "first" or "second"
SWAP=$5      # "swap" or empty

if [[ -z "$R1" || -z "$R2" || -z "$MAPPING_FILE" || -z "$MODE" ]]; then
    echo "Usage: $0 <R1.fastq[.gz]> <R2.fastq[.gz]> <barcode_mapping.txt> <first|second> [swap]"
    exit 1
fi

[[ "$R1" == *.gz ]] && R1_CMD="gzip -cd $R1" || R1_CMD="cat $R1"
[[ "$R2" == *.gz ]] && R2_CMD="gzip -cd $R2" || R2_CMD="cat $R2"

# Load barcode to name mapping into memory (awk associative array)
awk -v mode="$MODE" -v swap="$SWAP" -v mapfile="$MAPPING_FILE" '
BEGIN {
    FS = "\t"; OFS = "\n";

    # Load barcode mapping: barcode \t name
    while ((getline < mapfile) > 0) {
        BARCODE_MAP[$1] = $2;
    }
}
{
    line_num = (NR - 1) % 4 + 1;

    if (line_num == 1) {
        header1 = $1;
        header2 = $5;

        # Extract barcodes from R1 header
        match(header1, / [^ ]+$/, bc_match);
        split(bc_match[0], bc_parts, "[+]");
        bc1 = bc_parts[1];
        bc2 = bc_parts[2];

        if (swap == "swap") {
            tmp = bc1;
            bc1 = bc2;
            bc2 = tmp;
        }

        selected_bc = (mode == "second") ? bc2 : bc1;

        # Check if barcode is in mapping
        if (selected_bc in BARCODE_MAP) {
            sample = BARCODE_MAP[selected_bc];
            out1 = sample "_R1.fq";
            out2 = sample "_R2.fq";

            # Update header barcodes
            sub(/[+ ].*/, " " bc1 "+" bc2, header1);
            sub(/[+ ].*/, " " bc1 "+" bc2, header2);

            write = 1;
        } else {
            write = 0;
        }
    }

    if (write) {
        if (line_num == 1) { print header1 >> out1; print header2 >> out2; }
        if (line_num == 2) { print $2 >> out1; print $6 >> out2; }
        if (line_num == 3) { print $3 >> out1; print $7 >> out2; }
        if (line_num == 4) { print $4 >> out1; print $8 >> out2; }
    }
}' < <(paste <($R1_CMD) <($R2_CMD))

echo "Demultiplexing completed."