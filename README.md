# Pan-genome based Gene Family Analysis Pipeline

This pipeline provides a comprehensive approach for Pan-genome based gene family analysis by integrating BLASTP and HMM methods for screening candidate gene sequences. An optional feature allows users to construct a phylogenetic tree automatically. This README details the pipeline workflow, requirements, and usage instructions.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Usage](#usage)
- Pipeline Workflow
  - Parameter and Path Configuration
  - BLASTP Screening Module
  - HMM Screening Module
  - Sequence Extraction & HMM Rebuilding
  - Intersection of Candidate Genes
  - File Organization & Data Archiving
  - Optional Phylogenetic Tree Construction
  - Pan-genome Member Statistics (Optional)
- [Troubleshooting](#troubleshooting)
- [Contact Information](#contact-information)

---

## Overview

The gene family analysis pipeline is designed to identify and analyze candidate gene sequences by combining two independent screening methods:
1. **BLASTP Screening**: Compares input gene family FASTA sequences against a pre-built BLASTP database.
2. **HMM Screening**: Uses one or more HMM models to identify high-confidence candidate sequences in a protein database.

After the initial screening, candidate sequences are extracted and aligned using MAFFT. A new HMM model is then rebuilt to further refine candidate selection. Finally, the pipeline intersects the results from both methods, organizes the output files, and optionally constructs a phylogenetic tree. Pan-genome member statistics are generated based on a user-supplied sample ID file.

---

## Features

- **BLASTP Screening**: Filters candidate genes using an optional identity threshold.
- **HMM Screening**: Utilizes HMM models to screen for high-confidence candidate sequences (e-value < 1E-5).
- **Sequence Extraction & Alignment**: Extracts candidate sequences and performs multiple sequence alignment with MAFFT.
- **HMM Rebuilding**: Reconstructs a refined HMM model based on the aligned sequences.
- **Result Intersection**: Derives the final candidate gene set by intersecting BLASTP and HMM results.
- **File Organization**: Organizes intermediate and final output files into a user-specified directory.
- **Optional Phylogenetic Tree Construction**: Offers one-click tree construction using MAFFT and FastTree.
- **Pan-genome Statistics**: Computes candidate gene counts per sample based on a sample ID file (e.g., `110SampleID.txt`).

---

## Requirements

Ensure the following tools are installed and accessible in your system's PATH:
- **blast+** (from the NCBI BLAST+ suite)
- **hmmer** (includes `hmmscan`, `hmmsearch`, `hmmbuild`)
- **mafft** (multiple sequence alignment tool)
- **seqkit** (FASTA/Q file processing tool)
- **fasttree** (phylogenetic tree construction tool)

Additionally, prepare the following files:
- **Protein Database File**: e.g., `Pan110_xiaomi.fa`
- **BLASTP Database**: Pre-built BLASTP database files.
- **Gene Family FASTA File**: Input file for the gene family sequences (e.g., `AtSRS_10.fa`).
- **HMM Files**: One or more `.hmm` model files located in the script's directory.
- **Sample ID File**: e.g., `110SampleID.txt` for pan-genome member statistics.

---

## Usage

Run the script from the command line using:

```bash
bash genefamily_analysis.sh -i <gene_family_fasta_file> [-o <genefamily_folder_name>] [-t <threads>] [-m <i|u>] [-e <BLASTp_identity_threshold>] [yes/no]
```

-i <FASTA file>: Input FASTA file for gene family sequences.

-o <output folder>: The folder where output files will be saved (optional).

-t <threads>: Number of threads to use for parallel processing (optional, default is 1).

-m <mode>: Gene family mode (i for intersection, u for union, default is i).

-e <identity threshold>: BLASTP identity threshold (optional).

[yes/no]: Whether to automatically build a phylogenetic tree after gene family analysis. The default is no.

Example 1: Build the phylogenetic tree automatically:
```
bash genefamily_analysis.sh -i AtSRS_10.fa -o SRS -t 4 -m i -e 75 yes
```
Example 2: Skip tree construction:
```
bash genefamily_analysis.sh -i AtSRS_10.fa -o SRS -t 4 -m i -e 75 no
```
Example 3: Run without tree construction:
```
bash genefamily_analysis.sh -i AtSRS_10.fa
```

Troubleshooting
Environment Dependencies: Ensure all required tools (BLASTP, HMMER, MAFFT, Seqkit, FastTree) are installed and accessible in your system’s PATH.

Path Settings: Verify that the paths for PROTEIN_DB and BLASTP_DB are correctly set and match your local environment.

HMM Files: Ensure that all HMM model files are located in the same directory as the script.

Please Note: 
Sample ID File: if you wanna use genefamily Analysis Code for Pan-Genome, Please Confirm that the sample ID file (e.g., 110SampleID.txt) exists and is properly formatted.

Contact Information
Author: Ruimiao Li

Email: jiekexiaobai99@gmail.com or 1205630141@qq.com

Version: V1.3 (2025-04-26)

For any issues, questions, or further assistance, please contact the author via the provided email address.
