#!/bin/bash

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <input BED file> <genome FASTA file> <genome size file>"
  exit 1
fi

# Input files
INPUT_BED=$1
GENOME_FASTA=$2
GENOME_SIZE_FILE=$3

# Output prefixes
FILTERED_BED_ONLY_READ1="${INPUT_BED%.bed}_only_read1.bed"
FLANKING_BED="${INPUT_BED%.bed}_flanking_regions.bed"
FLANKING_FASTA="${INPUT_BED%.bed}_flanking_sequences.fa"
CLEANED_FASTA="${INPUT_BED%.bed}_cleaned_flanking_sequences.fa"
FINAL_FASTA="${INPUT_BED%.bed}_final_flanking_sequences.fa"
FINAL_RC_FASTA="${INPUT_BED%.bed}_final_flanking_sequences_rc.fa"
DIPY_READS="${INPUT_BED%.bed}_reads_with_dipyrimidines.txt"
PREFIXED_DIPY="${INPUT_BED%.bed}_prefixed_reads_with_dipyrimidines.txt"
FILTERED_BED="${INPUT_BED%.bed}_dipy_filtered_reads.bed"
FILTERED_BAM="${INPUT_BED%.bed}_dipy_filtered_reads.bam"

###############################################################################
# Step 0: Filter the BED file to retain only entries with /1
###############################################################################
awk '$4 ~ /\/1$/' "$INPUT_BED" > "$FILTERED_BED_ONLY_READ1"
echo "Step 0: Filtered BED file with only /1 reads: $FILTERED_BED_ONLY_READ1"

###############################################################################
# Step 1: Generate 2 base pairs flanking regions (both upstream and downstream)
###############################################################################
bedtools flank -i "$FILTERED_BED_ONLY_READ1" -g "$GENOME_SIZE_FILE" -l 2 -r 0 -s > "$FLANKING_BED"
echo "Step 1: Flanking regions generated: $FLANKING_BED"

###############################################################################
# Step 2: Extract sequences for the flanking regions
###############################################################################
bedtools getfasta -fi "$GENOME_FASTA" -bed "$FLANKING_BED" -name -s -fo "$FLANKING_FASTA"
echo "Step 2: FASTA sequences extracted: $FLANKING_FASTA"

###############################################################################
# Step 3: Clean up the sequence headers
###############################################################################
sed 's/::.*//' "$FLANKING_FASTA" > "$CLEANED_FASTA"
echo "Step 3: Sequence headers cleaned: $CLEANED_FASTA"

###############################################################################
# Step 4: Convert all sequences to uppercase
###############################################################################
awk '{ 
  if (NR % 2 == 0) 
    print toupper($0); 
  else 
    print $0 
}' "$CLEANED_FASTA" > "$FINAL_FASTA"
echo "Step 4: Sequences converted to uppercase: $FINAL_FASTA"

###############################################################################
# Step 5: Reverse-complement the uppercase sequences
###############################################################################
awk '
  function revcomp(seq) {
    rev = ""
    for (i = length(seq); i > 0; i--) {
      c = substr(seq, i, 1)
      # Complement
      if      (c == "A") c = "T"
      else if (c == "T") c = "A"
      else if (c == "C") c = "G"
      else if (c == "G") c = "C"
      rev = rev c
    }
    return rev
  }
  NR % 2 == 1 {
    # Header line (e.g., >readName)
    print $0
  }
  NR % 2 == 0 {
    # Sequence line
    print revcomp($0)
  }
' "$FINAL_FASTA" > "$FINAL_RC_FASTA"
echo "Step 5: Sequences reverse-complemented: $FINAL_RC_FASTA"

###############################################################################
# Step 6: Identify reads with dipyrimidines (TT, CC, CT, TC)
#         (Now operating on the reverse-complemented FASTA)
###############################################################################
awk '
  NR % 2 == 1 {
    header=$0
  }
  NR % 2 == 0 && /TT|CC|CT|TC/ {
    print header
  }
' "$FINAL_RC_FASTA" \
  | sed 's/>//' \
  > "$DIPY_READS"

echo "Step 6: Reads with dipyrimidines identified: $DIPY_READS"

###############################################################################
# Step 7: Remove /1 or /2 suffix from read names
###############################################################################
sed 's:/[12]$::' "$DIPY_READS" > "$PREFIXED_DIPY"
echo "Step 7: Suffixes removed from read names: $PREFIXED_DIPY"

###############################################################################
# Step 8: Filter BED file to retain reads with dipyrimidines and their mates
###############################################################################
awk '
  NR==FNR {
    a[$1] = 1;
    next
  }
  {
    prefix = substr($4, 1, length($4)-2);
    if (prefix in a) {
      print
    }
  }
' "$PREFIXED_DIPY" "$INPUT_BED" > "$FILTERED_BED"
echo "Step 8: Filtered BED file created: $FILTERED_BED"

###############################################################################
# Step 9: Intersect BAM file with filtered BED file
###############################################################################
bedtools intersect -a "${INPUT_BED%.bed}.bam" -b "$FILTERED_BED" -u > "$FILTERED_BAM"
echo "Step 9: Filtered BAM file created: $FILTERED_BAM"

echo "Pipeline completed successfully."

####################################
# To run this file ./Dipyfilter.sh IMR_100J_UVC_1_aligned_outputE2E_dedup_mapq.bed hg38.fa hg38.chrom.sizes
####################################

