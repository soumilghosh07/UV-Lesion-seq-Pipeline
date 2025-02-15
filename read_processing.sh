#!/bin/bash



# Usage: ./process_reads.sh <reference_index> <R1.fastq> <R2.fastq> <output_prefix> <threads>

# Example: ./process_reads.sh GRCh38_noalt_as sample_R1.fastq sample_R2.fastq sample_output 16



# Input arguments

REFERENCE_INDEX=$1

R1_FASTQ=$2

R2_FASTQ=$3

OUTPUT_PREFIX=$4

THREADS=$5



# Load required modules

module load biology

module load samtools

module load bowtie2

module load umi_tools



# Alignment with Bowtie2

bowtie2 --end-to-end --very-sensitive -p "$THREADS" -x "$REFERENCE_INDEX" -1 "$R1_FASTQ" -2 "$R2_FASTQ" -S "${OUTPUT_PREFIX}_aligned_outputE2E.sam" 2> "${OUTPUT_PREFIX}_bowtie2_statsE2E.txt"



# Convert SAM to BAM

samtools view -bS "${OUTPUT_PREFIX}_aligned_outputE2E.sam" > "${OUTPUT_PREFIX}_aligned_outputE2E.bam"



# Fixmate

samtools fixmate -m "${OUTPUT_PREFIX}_aligned_outputE2E.bam" "${OUTPUT_PREFIX}_aligned_outputE2E_fixed.bam"



# Sort BAM file

samtools sort -@ "$THREADS" -o "${OUTPUT_PREFIX}_aligned_outputE2E_sorted.bam" "${OUTPUT_PREFIX}_aligned_outputE2E_fixed.bam"



# Index BAM file

samtools index -@ "$THREADS" "${OUTPUT_PREFIX}_aligned_outputE2E_sorted.bam"



# Deduplicate using UMI-tools

umi_tools dedup --extract-umi-method=read_id --umi-separator=: \

    -I "${OUTPUT_PREFIX}_aligned_outputE2E_sorted.bam" \

    -S "${OUTPUT_PREFIX}_aligned_outputE2E_dedup.bam" \

    --paired --chimeric-pairs=discard --unpaired-reads=discard



# Remove intermediate sorted BAM and index files

rm "${OUTPUT_PREFIX}_aligned_outputE2E_sorted.bam" "${OUTPUT_PREFIX}_aligned_outputE2E_sorted.bai"



# Filter reads by MAPQ

samtools view -q 20 -b -@ "$THREADS" "${OUTPUT_PREFIX}_aligned_outputE2E_dedup.bam" > "${OUTPUT_PREFIX}_aligned_outputE2E_dedup_mapq.bam"



# Index filtered BAM

samtools index -@ "$THREADS" "${OUTPUT_PREFIX}_aligned_outputE2E_dedup_mapq.bam"



# Count and calculate reads filtered out

INPUT_READS=$(samtools view -c "${OUTPUT_PREFIX}_aligned_outputE2E_dedup.bam")

OUTPUT_READS=$(samtools view -c "${OUTPUT_PREFIX}_aligned_outputE2E_dedup_mapq.bam")

FILTERED_OUT_READS=$((INPUT_READS - OUTPUT_READS))



echo "Number of reads filtered out due to MAPQ < 20: $FILTERED_OUT_READS"



# Generate stats

samtools idxstats "${OUTPUT_PREFIX}_aligned_outputE2E_dedup_mapq.bam" > "${OUTPUT_PREFIX}_aligned_outputE2E_dedup_mapq_idxstats.txt"

samtools flagstat "${OUTPUT_PREFIX}_aligned_outputE2E_dedup_mapq.bam" > "${OUTPUT_PREFIX}_aligned_outputE2E_dedup_mapq_flagstat.txt"



# Clean up

rm "${OUTPUT_PREFIX}_aligned_outputE2E.sam"



echo "Finished processing of ${OUTPUT_PREFIX} at $(date)"


