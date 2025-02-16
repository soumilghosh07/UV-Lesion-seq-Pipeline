# 🧬 **Multi-Step Sequencing Data Processing Pipeline**
This repository contains a complete **sequencing data processing pipeline**, combining **Nextflow-based workflows and Bash scripts** to handle:
- **Demultiplexing & Adapter Trimming** (Nextflow)
- **Read Processing** (Bowtie2 alignment, deduplication, filtering)
- **Dipyrimidine Filtering** (Identifying reads having undergone UV lesions in sequencing reads)

---

## 📂 **Repository Contents**
| File | Description |
|------|------------|
| `demultiplexing&adapterTrimming.nf` | Nextflow script for demultiplexing & adapter trimming |
| `read_processing.sh` | Bash script for read processing (alignment, deduplication, filtering) |
| `DipyrimidineFilter.sh` | Bash script for dipyrimidine filtering in sequencing reads |

---

## 🚀 **Installation & Dependencies**
This pipeline requires the following tools:

- **[Nextflow](https://www.nextflow.io/)** (workflow automation)
- **[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/)** (read alignment)
- **[SAMtools](http://www.htslib.org/)** (BAM file processing)
- **[UMI-tools](https://umi-tools.readthedocs.io/en/latest/)** (deduplication)
- **[BEDTools](https://bedtools.readthedocs.io/en/latest/)** (BED file processing)
- **[Cutadapt](https://cutadapt.readthedocs.io/en/stable/)** (adapter trimming)
