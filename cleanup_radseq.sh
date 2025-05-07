#!/bin/bash

# cleanup_radseq.sh
# Removes all intermediate files from ddRADseq processing.
# Keeps only:
# - Filtered.1.*.fq.gz
# - Filtered.2.*.fq.gz
# - Rescued.1.*.fq.gz
# - Rescued.2.*.fq.gz
# - Short.*.fq and Short.*.fq.gz
# - Counts*.stat

set -x  # Enable debug mode to show all executed commands

show_help() {
  echo "Usage: ./cleanup_radseq.sh [directory] [--dry-run|--help]"
  echo ""
  echo "Description:"
  echo "  This script deletes all intermediate files from a ddRADseq processing folder,"
  echo "  keeping only the essential output files:"
  echo "    - Filtered.1.*.fq.gz"
  echo "    - Filtered.2.*.fq.gz"
  echo "    - Rescued.1.*.fq.gz"
  echo "    - Rescued.2.*.fq.gz"
  echo "    - Short.*.fq and Short.*.fq.gz"
  echo "    - Counts*.stat"
  echo ""
  echo "Arguments:"
  echo "  [directory]     Optional. Path to the output directory (default: current directory)"
  echo "  --dry-run       Show what would be deleted without removing anything"
  echo "  --help          Show this help message and exit"
  echo ""
  echo "Examples:"
  echo "  ./cleanup_radseq.sh ./results/"
  echo "  ./cleanup_radseq.sh ./results/ --dry-run"
  echo ""
}

# Defaults
directory="."
dry_run=false

# Parse args
for arg in "$@"; do
  case "$arg" in
    --help)
      show_help
      exit 0
      ;;
    --dry-run)
      dry_run=true
      ;;
    *)
      directory="$arg"
      ;;
  esac
done

echo "Target directory: $directory"
echo "Dry run mode: $dry_run"

if [ "$dry_run" = true ]; then
  echo "Files that would be deleted:"
  find "$directory" -maxdepth 1 -type f \
    ! -name 'Filtered.1.*.fq.gz' \
    ! -name 'Filtered.2.*.fq.gz' \
    ! -name 'Rescued.1.*.fq.gz' \
    ! -name 'Rescued.2.*.fq.gz' \
    ! -name 'Short.*.fq' \
    ! -name 'Short.*.fq.gz' \
    ! -name 'Counts*.stat'
else
  echo "Deleting intermediate files from $directory..."
  find "$directory" -maxdepth 1 -type f \
    ! -name 'Filtered.1.*.fq.gz' \
    ! -name 'Filtered.2.*.fq.gz' \
    ! -name 'Rescued.1.*.fq.gz' \
    ! -name 'Rescued.2.*.fq.gz' \
    ! -name 'Short.*.fq' \
    ! -name 'Short.*.fq.gz' \
    ! -name 'Counts*.stat' \
    -exec rm -v {} +

  # Also remove common temporary files in the target directory
  rm -f "$directory"/*.count "$directory"/*.txt "$directory"/*.list "$directory"/*.table "$directory"/first.column 2>/dev/null

  echo "Cleanup complete. Essential files retained."
fi