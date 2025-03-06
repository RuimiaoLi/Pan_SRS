Gene Family Analysis Pipeline
Overview
This pipeline provides a comprehensive approach to gene family analysis by integrating both BLASTP and HMM methods for screening candidate gene sequences. In addition, it offers an optional one-click phylogenetic tree construction feature. The pipeline is designed to facilitate the identification, alignment, and analysis of gene family members from a given protein database.

Features
BLASTP Screening:
Compares user-supplied gene family FASTA sequences against a pre-built BLASTP database to identify candidate genes. The results are filtered based on an identity threshold (optional parameter).

HMM Screening:
Uses one or more HMM models to search the protein database. Candidate sequences with high confidence (e-value < 1E-5) are extracted. If multiple HMM models are available, their results are intersected to increase specificity.

Sequence Extraction & HMM Rebuilding:
Candidate sequences are extracted from the protein database. These sequences are then aligned using MAFFT, and a new HMM model is rebuilt to refine candidate gene selection.

Intersection of Results:
The pipeline takes the intersection of the BLASTP and HMM results to obtain a final candidate gene set, ensuring higher accuracy in the screening process.

File Organization & Data Archiving:
Intermediate and final output files are organized into a user-specified folder. Raw data files (FASTA and HMM models) are moved to a designated directory for easy management and reproducibility.

Optional Phylogenetic Tree Construction:
Users have the option to build a phylogenetic tree from the final candidate gene sequences using MAFFT for alignment and FastTree for tree construction.

Pan-genome Member Statistics:
Based on a pre-supplied sample ID file (e.g., 110SampleID.txt), the pipeline calculates and categorizes the number of candidate genes present in each sample.

Requirements
Ensure that the following tools are installed and accessible from your system’s PATH:

BLASTP (from the NCBI BLAST+ suite)
HMMER (includes hmmscan, hmmsearch, and hmmbuild)
MAFFT (for multiple sequence alignment)
Seqkit (for FASTA/Q file manipulation)
FastTree (for phylogenetic tree construction)
Additionally, make sure you have the necessary input files:

Protein Database File: A FASTA file containing all protein sequences (e.g., Pan110_xiaomi.fa).
BLASTP Database: Pre-built BLASTP database files.
Gene Family FASTA File: The user-supplied file containing sequences for the gene family of interest.
HMM Files: One or more HMM model files (with .hmm extension) placed in the same directory as the script.
Sample ID File: A file (e.g., 110SampleID.txt) containing sample IDs for pan-genome statistics.
Usage
Run the script from the command line with the following format:

bash
复制
编辑
bash genefamily_analysis.sh <gene_family_folder_name> <gene_family_fasta_file> [BLASTP_identity_threshold]
Example
bash
复制
编辑
bash genefamily_analysis.sh SRS AtSRS_10.fa 30
SRS: The name of the folder where the output files will be stored.
AtSRS_10.fa: The FASTA file containing the gene family sequences to be analyzed.
30: (Optional) The identity threshold used for filtering BLASTP results. If not provided, the default value is 0.
During execution, the script will:

Check Environment: Verify that all required tools are installed.
Run BLASTP: Perform sequence comparisons and filter results based on the identity threshold.
Run HMM Analysis: Process each HMM file to extract high-confidence candidate sequences.
Extract Sequences: Retrieve candidate sequences from the protein database using the filtered results.
Perform Multiple Sequence Alignment: Align the candidate sequences with MAFFT.
Rebuild HMM: Reconstruct the HMM model based on the aligned sequences.
Intersect Results: Obtain the common candidate genes identified by both BLASTP and HMM methods.
Organize Output Files: Move intermediate and final output files to the designated folder.
Optional Phylogenetic Tree Construction: Ask the user if a phylogenetic tree should be built.
Pan-genome Statistics: Compute the number of candidate genes for each sample based on the provided sample ID file.
Pipeline Modules Explained
1. Parameter and Path Configuration
PROTEIN_DB & BLASTP_DB:
Define the paths to your protein FASTA file and BLASTP database, respectively.
Color Variables:
Used for enhancing terminal output with colors to indicate success or errors.
2. BLASTP Screening Module
Function: run_blastp()
Description:
Executes BLASTP using the gene family FASTA file as a query against the specified database. The output is filtered to create a candidate gene list based on an optional identity threshold.
3. HMM Screening Module
Function: run_hmm()
Description:
Iterates through available HMM files and runs hmmsearch on the protein database. High-confidence candidate sequences (e-value < 1E-5) are recorded for each HMM. If multiple HMM files are present, the script intersects their results to generate a more accurate candidate list.
4. Sequence Extraction and HMM Rebuilding
Functions:
extract_sequences_4_rebuild_HMM(), multi_align_seq(), and rebuild_hmm()
Description:
The candidate sequences identified by the HMM screening are extracted from the protein database. These sequences are then aligned using MAFFT, and a new HMM model is rebuilt using the alignment. A second round of hmmsearch further refines the candidate gene set.
5. Intersection of Candidate Gene Sets
Function: common_method()
Description:
Takes the intersection of the candidate lists produced by BLASTP and HMM methods, ensuring that only sequences identified by both approaches are considered in the final candidate set. The common candidate sequences are then extracted from the protein database using Seqkit.
6. File Organization & Data Archiving
Functions:
move_files() and move_rawData()
Description:
Organizes and moves all intermediate and final output files into a user-specified folder. Raw input files (gene family FASTA and HMM files) are archived in a designated subfolder for future reference.
7. Optional Phylogenetic Tree Construction
Functions:
ask_to_build_tree() and build_tree()
Description:
After candidate gene identification, the user is prompted to choose whether to construct a phylogenetic tree. If confirmed, the pipeline uses MAFFT for multiple sequence alignment and FastTree for tree construction, outputting the tree in Newick format.
8. Pan-genome Member Statistics
Function: count_Pan_members()
Description:
Reads a sample ID file to count and classify the candidate gene occurrences across samples. The script outputs a statistics file that categorizes samples based on whether they have more, less, or equal numbers of candidate genes compared to the mode (most common count).
Troubleshooting
Environment Dependencies:
If any required tool (BLASTP, HMMER, MAFFT, Seqkit, or FastTree) is not found, the script will terminate and display an error message. Verify that each tool is installed and accessible in your system PATH.

File Paths:
Ensure that the paths for PROTEIN_DB and BLASTP_DB are correctly set according to your directory structure.

HMM Files:
All HMM model files must reside in the same directory as the script to be processed automatically.

Sample ID File:
Make sure that the sample ID file (e.g., 110SampleID.txt) exists in the current directory and contains valid sample IDs for the pan-genome statistics.

Contact Information
Author: Ruimiao Li
Email: 1205630141@qq.com
Version: V1.3 (2024-11-19)
If you encounter any issues or have further questions regarding the pipeline, please feel free to reach out via the provided contact email.
