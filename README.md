# RADseq_reads
pipeline for reads preprocessing from RADseq experiments done similar to one described in Schweyen et al. 2014 - but can be freely modified for diferenbt versions.
pipeline inspired by process_radtags pipeline from STACKS, outputing read files to be used by next STACKS steps.
pipeline uses cutadapt 5.0
when using specific conda envinorment like here for cutadapt 5.0 pipeline needs to be run: source <<file name>> to enable loading this specific envinorment. When not using CONDA you can comment the line with
conda activate cutadapt

You need to input proper PATH TO READS TEMPLATE at "$reads" variable

You need to input proper BARCODES file name at "$barcodes" variable

Default variables to be modifies by user prior to running the script:
# ***
reads="/mnt/qnap/projects/RADseq/MPabijan_RAD/X204SC23021387-Z01-F001/01.RawData/MP/MP_EKDL230001910-1A_HNW37DSX5_L2_" # template for reads
barcodes="old_barcodes_cutadapt.txt" # fasta file with barcodes
# ***
qualada="20,20" # quality trimming adapters
qualfil="20,20" # quality trimming filtering
errbar=1 # number of allowed errors for demultiplexing
errfil=0 # number of allowed errors for CUT site filtering
errfilanc=1 # number of allowed error for anchored CUT site filtering
thr=6 # number of processor threads to use
lenada=130 # minimum length for adapter removal
lenfil=100 # minimum length for filtering cut sites





