#!/usr/bin/env nextflow

// PARAMETERS
params.R1    = "/scratch/users/soumil/volcano2/wi38merged/align_outputs/raw/Undetermined_HT7VYDSXX_L2_1.fq"
params.R2    = "/scratch/users/soumil/volcano2/wi38merged/align_outputs/raw/Undetermined_HT7VYDSXX_L2_2.fq"
params.adapters = "/scratch/users/soumil/volcano2/wi38merged/align_outputs/raw/adapters.fa"
params.results = "/scratch/users/soumil/volcano2/wi38merged/align_outputs/raw/results"

log.info """\
================================================================================
DEMULTIPLEX & TRIM PIPELINE
================================================================================

|--------|
| INPUTS |
|--------|

 R1 fastq    : $params.R1
 R2 fastq    : $params.R2
 adapters fasta : $params.adapters
 results dir  : $params.results

================================================================================
"""

// WORKFLOW
workflow {

  // Input channels
  Channel
    .fromPath([params.R1, params.R2], checkIfExists: true)
    .collect()
    .set { reads }

  Channel
    .fromPath(params.adapters, checkIfExists: true)
    .set { adapters }

  // Run demultiplex / adapter trimming
  DEMULTIPLEX(reads, adapters)
}


// PROCESS
process DEMULTIPLEX {
  publishDir "${params.results}", mode: "copy"
  cpus 8

  input:
  tuple path(fastq_R1), path(fastq_R2)
  path adapters

  output:
  path "fastq/*R1.fastq", emit: fq1Files
  path "fastq/*R2.fastq", emit: fq2Files
  path "log/demux/*",  emit: logDemux

  """
  mkdir -p fastq
  mkdir -p log/demux

  # First trimming step
  cutadapt --json=log/demux/fragment_construction.cutadapt.json \
   --cores=${task.cpus} -m 50:55 \
   -a CTGTCTCTTATACACATCTAGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
   -A "CACTGCNNNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;max_error_rate=0.15;min_overlap=11" \
   -o R1_trim.fastq -p R2_trim.fastq \
   ${fastq_R1} ${fastq_R2}

  # Second trimming step using adapter files
  cutadapt --cores=${task.cpus} --no-indels --discard-untrimmed \
   -g file:${adapters} \
   -G AGATGTGTATAAGAGACAG \
   -o fastq/{name}.R1.fastq -p fastq/{name}.R2.fastq \
   R1_trim.fastq R2_trim.fastq
  """
}
