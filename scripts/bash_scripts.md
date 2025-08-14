## Table of Contents
  - [Table of Contents](#table-of-contents)
- [Artemisia Database Development: Bash Scripts](#artemisia-database-development-bash-scripts)
  - [Conversion of SRA Files to FASTQ Using fastq-dump (3.1.0)](#conversion-of-sra-files-to-fastq-using-fastq-dump-310)
  - [Converting Paired-End SRA Files to FASTQ](#converting-paired-end-sra-files-to-fastq)
  - [Converting Single-End SRA Files to FASTQ](#converting-single-end-sra-files-to-fastq)
  - [Compressing FASTQ Files to FASTQ.GZ](#compressing-fastq-files-to-fastqgz)
  - [Deleting Non-GZ Files](#deleting-non-gz-files)
  - [Quality Control and Adapter Removal Using fastp (0.23.2)](#quality-control-and-adapter-removal-using-fastp-0232)
  - [Quality Control for Short-Read Data (Paired-End and Single-End)](#quality-control-for-short-read-data-paired-end-and-single-end)
  - [Extracting Quality Metrics from fastp JSON Files](#extracting-quality-metrics-from-fastp-json-files)
  - [Quality Control for Long-Read RNA Sequencing Data](#quality-control-for-long-read-rna-sequencing-data)
  - [Alignment of RNA Sequencing Data](#alignment-of-rna-sequencing-data)
  - [Installation of Tools](#installation-of-tools)
  - [Alignment of Short-Read RNA Sequencing Data Using HISAT2](#alignment-of-short-read-rna-sequencing-data-using-hisat2)
  - [Alignment of Long-Read RNA Sequencing Data Using minimap2](#alignment-of-long-read-rna-sequencing-data-using-minimap2)
  - [Filtering and Transcriptome Assembly of RNA Sequencing Data](#filtering-and-transcriptome-assembly-of-rna-sequencing-data)
  - [Installation of Tools](#installation-of-tools)
  - [Filtering Uniquely Mapped Reads Using samtools](#filtering-uniquely-mapped-reads-using-samtools)
  - [Merging Long-Read Alignments](#merging-long-read-alignments)
  - [Transcriptome Assembly of Short-Read Sequences Using StringTie](#transcriptome-assembly-of-short-read-sequences-using-stringtie)
  - [Merging Short-Read StringTie Assemblies](#merging-short-read-stringtie-assemblies)
  - [Transcriptome Assembly of Long-Read Sequences Using StringTie](#transcriptome-assembly-of-long-read-sequences-using-stringtie)
  - [Transcriptome Assembly of Long-Read Sequences Using IsoQuant](#transcriptome-assembly-of-long-read-sequences-using-isoquant)
  - [Indexing the Genome for Faster Processing](#indexing-the-genome-for-faster-processing)
  - [Extracting Reference Transcriptome](#extracting-reference-transcriptome)
  - [Transcriptome Quality Assessment Using BUSCO (v5.7.1)](#transcriptome-quality-assessment-using-busco-v571)
  - [Installation of BUSCO](#installation-of-busco)
  - [BUSCO Analysis for Reference Genome](#busco-analysis-for-reference-genome)
  - [BUSCO Analysis for Transcriptome Assemblies](#busco-analysis-for-transcriptome-assemblies)

# Artemisia Database Development: Bash Scripts

## Conversion of SRA Files to FASTQ Using fastq-dump (3.1.0)
### Installation of sra-tools

Due to issues with an older version of `sra-tools` (2.8.0), version 3.1.0 is recommended. Download the latest version from Bioconda and install it using the following command:

```bash
conda install sra-tools-3.1.0-h9f5acd7_0.tar.bz2
```

## Converting Paired-End SRA Files to FASTQ

The following script converts paired-end SRA files to FASTQ format using `fastq-dump` with the `--split-files` option.

### `fastq-dump_PE.sh`

```bash
#!/bin/bash

#SBATCH --job-name=fastq-dump_PE
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mypy

# Define directories containing SRA files
directories=(
  "PRJNA1089479"
  "PRJNA371708"
  "PRJNA636791"
  "PRJNA752933"
  "PRJNA796810"
)

for dir in "${directories[@]}"; do
  files="<path_to_fastq_files>/$dir/SRR*"
  cd "<path_to_fastq_files>/$dir"
  fastq-dump --split-files $files
  cd -
done
```

**Notes**:
- Replace `<path_to_fastq_files>` with the path to your FASTQ files directory.
- The `--split-files` option ensures paired-end reads are split into separate FASTQ files.

## Converting Single-End SRA Files to FASTQ

The following script converts single-end SRA files to FASTQ format using `fastq-dump`.

### `fastq-dump_SE.sh`

```bash
#!/bin/bash

#SBATCH --job-name=fastq-dump_SE
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mypy

# Define directories containing SRA files
directories=(
  "PRJNA524494"
)

for dir in "${directories[@]}"; do
  files="<path_to_fastq_files>/$dir/SRR*"
  cd "<path_to_fastq_files>/$dir"
  fastq-dump $files
  cd -
done
```



## Compressing FASTQ Files to FASTQ.GZ

The following script compresses FASTQ files to `.gz` format using `gzip`.

### `gzip_fastq.sh`

```bash
#!/bin/bash

#SBATCH --job-name=gzip
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3

# Define directories containing FASTQ files
directories=(
  "PRJNA996091"
)

for dir in "${directories[@]}"; do
  files="<path_to_fastq_files>/$dir/*.fastq"
  cd "<path_to_fastq_files>/$dir"
  gzip -k --fast $files
  cd -
done
```

**Notes**:
- The `-k` option keeps the original `.fastq` files, and `--fast` prioritizes compression speed.

## Deleting Non-GZ Files

To remove all files in a directory except those with the `.gz` extension, use the following command:

```bash
find . -type f ! -name '*.gz' -delete
```

**Notes**:
- Run this command in the directory containing the files you wish to clean up.
- Use with caution, as it permanently deletes files without moving them to a trash/recycle bin.

## Quality Control and Adapter Removal Using fastp (0.23.2)

### Installation of fastp

Ensure `fastp` version 0.23.2 is installed in your environment. You can install it using Conda:

```bash
conda install fastp=0.23.2
```

Activate the environment before running the scripts:

```bash
source activate fastp
```

## Quality Control for Short-Read Data (Paired-End and Single-End)

The following script processes paired-end (PE) and single-end (SE) FASTQ files, removing reads with a mean read length < 40 bp or Q20 rate < 80%.

### `fastp.sh`

```bash
#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate fastp

# Define the main directory containing FASTQ files
main_dir="<path_to_fastq_files>"

# Define the output directory for quality-controlled files
output_dir="<path_to_qc_dir>"

# Define paired-end directories
pe_dir=(
  "PRJNA996091"
)

# Define single-end directories
se_dir=(
  "PRJNA524494"
  "PRJNA755452"
)

# Loop through the paired-end directories
for dir_name in "${pe_dir[@]}"; do
    # Construct the subdirectory path
    sub_dir="${main_dir}/${dir_name}"

    # Create the output subdirectory if it doesn't exist
    mkdir -p "${output_dir}/${dir_name}"

    # Loop through the files in the subdirectory
    for file_1 in "${sub_dir}"/*_1.fastq.gz; do
        # Get the file name without extension
        file_base=$(basename "$file_1" _1.fastq.gz)

        # Construct the corresponding file path for _2.fastq.gz
        file_2="${sub_dir}/${file_base}_2.fastq.gz"

        # Run fastp for the current files
        fastp -i "${file_1}" \
              -I "${file_2}" \
              -o "${output_dir}/${dir_name}/${file_base}_1.fastq.gz" \
              -O "${output_dir}/${dir_name}/${file_base}_2.fastq.gz" \
              -q 20 -l 40 --thread 8 --trim_front1 15 --trim_front2 15 --detect_adapter_for_pe -- fasseval_duplication \
              --html "${output_dir}/${dir_name}/${file_base}.html" \
              --json "${output_dir}/${dir_name}/${file_base}.json"
    done
done

# Loop through the single-end directories
for dir in "${se_dir[@]}"; do
    # Create the output directory for the current single-end directory
    mkdir -p "${output_dir}/${dir}"

    # Loop through the files in the current single-end directory
    for file in "${main_dir}/${dir}"/*.fastq.gz; do
        # Get the file name without extension
        file_base=$(basename "$file" .fastq.gz)

        # Run fastp for the current file
        fastp -i "${file}" \
              -o "${output_dir}/${dir}/${file_base}.fastq.gz" \
              -q 20 --thread 8 --dont_eval_duplication \
              --html "${output_dir}/${dir}/${file_base}.html" \
              --json "${output_dir}/${dir}/${file_base}.json"
    done
done
```

**Notes**:
- The `-q 20` option filters reads with a Q20 rate < 80%, and `-l 40` removes reads shorter than 40 bp.
- The `--dont_eval_duplication` option skips duplication evaluation to save time.
- For paired-end data, `--trim_front1 15` and `--trim_front2 15` trim 15 bases from the 5' end, and `--detect_adapter_for_pe` enables adapter detection.

### File Management

After running `fastp`, move the filtered `.fastq.gz` files to a designated directory (e.g., `02_filtered_FASTQ`):

```bash
cp -r <path_to_qc_dir>/* <path_to_filtered_fastq>/
```

Remove `.fastq.gz` files from the quality control directory, keeping only the HTML and JSON reports:

```bash
find <path_to_qc_dir> -type f -name "*.fastq.gz" -exec rm {} \;
```

Remove HTML and JSON reports from the filtered FASTQ directory:

```bash
find <path_to_filtered_fastq> -type f \( -name "*.html" -o -name "*.json" \) -exec rm {} \;
```

**Notes**:
- Use the `find` commands with caution, as they permanently delete files.

## Extracting Quality Metrics from fastp JSON Files

The following script extracts quality metrics from `fastp` JSON reports and compiles them into a CSV file.

### `fastp_status.sh`

```bash
#!/bin/bash

# Set the main directory containing JSON reports
main_dir="<path_to_qc_dir>"

# Create an empty CSV file
output_file="output.csv"
echo "sample,total_reads_before,total_bases_before,q20_rate_before,total_reads_after,total_bases_after,q20_rate_after" > "$output_file"

# Loop through each directory in the main directory
for dir in "$main_dir"/PRJ*; do
    if [ -d "$dir" ]; then
        # Find all the JSON files in the directory
        json_files=($(find "$dir" -type f -name "*.json"))

        for json_file in "${json_files[@]}"; do
            # Extract the sample name from the base name of the JSON file
            sample=$(basename "$json_file" .json)

            # Extract the desired information from the JSON file
            total_reads_before=$(jq -r '.summary.before_filtering.total_reads' "$json_file")
            total_bases_before=$(jq -r '.summary.before_filtering.total_bases' "$json_file")
            q20_rate_before=$(jq -r '.summary.before_filtering.q20_rate' "$json_file")
            total_reads_after=$(jq -r '.summary.after_filtering.total_reads' "$json_file")
            total_bases_after=$(jq -r '.summary.after_filtering.total_bases' "$json_file")
            q20_rate_after=$(jq -r '.summary.after_filtering.q20_rate' "$json_file")

            # Append the information to the CSV file
            echo "$sample,$total_reads_before,$total_bases_before,$q20_rate_before,$total_reads_after,$total_bases_after,$q20_rate_after" >> "$output_file"
        done
    fi
done

echo "CSV file created: $output_file"
```

**Notes**:
- Requires `jq` to parse JSON files. Install it if necessary (e.g., `conda install jq`).
- The output CSV file (`output.csv`) contains metrics such as total reads, bases, and Q20 rates before and after filtering.


## Quality Control for Long-Read RNA Sequencing Data

The following script processes long-read single-end FASTQ files, removing reads with a mean read length < 40 bp or Q20 rate < 80%.

### `fastp-long.sh`

```bash
#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate fastp

# Define the main directory containing long-read FASTQ files
main_dir="<path_to_fastq_files>/longreads"

# Define the output directory for quality-controlled files
output_dir="<path_to_qc_dir>"

# Define single-end directories
se_dir=(
  "PRJNA752933"
)

# Loop through the single-end directories
for dir in "${se_dir[@]}"; do
    # Create the output directory for the current single-end directory
    mkdir -p "${output_dir}/${dir}"

    # Loop through the files in the current single-end directory
    for file in "${main_dir}/${dir}"/*.fastq.gz; do
        # Get the file name without extension
        file_base=$(basename "$file" .fastq.gz)

        # Run fastp for the current file
        fastp -i "${file}" \
              -o "${output_dir}/${dir}/${file_base}.fastq.gz" \
              -q 20 -f 30 --thread 8 --dont_eval_duplication \
              --html "${output_dir}/${dir}/${file_base}.html" \
              --json "${output_dir}/${dir}/${file_base}.json"
    done
done
```

**Notes**:
- The `-f 30` option trims 30 bases from the 5' end of long reads.
- The `-q 20` option filters reads with a Q20 rate < 80%.


## Alignment of RNA Sequencing Data

## Installation of Tools

```bash
conda install hisat2=2.2.1 samtools minimap2
source activate mypy
```

## Alignment of Short-Read RNA Sequencing Data Using HISAT2

### Building the Genome Index

The following script builds the HISAT2 genome index for the *Artemisia* genome.

#### `hisat2-build.sh`

```bash
#!/bin/bash
#SBATCH --job-name=hisat2-build
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mypy

hisat2-build -p 32 <path_to_genome>/LQ-9_phase0_genome.fasta <path_to_genome>/LQ-9_phase0_genome_hisat2.idx
```

### Aligning Paired-End (PE) Short-Read Data

The following script aligns paired-end FASTQ files to the *Artemisia* genome using HISAT2.

#### `hisat2_PE.sh`

```bash
#!/bin/bash
#SBATCH --job-name=hisat2_PE
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mypy

# Set the paths
genome_index="<path_to_genome>/LQ-9_phase0_genome_hisat2.idx"
input_dir="<path_to_filtered_fastq>"
output_dir="<path_to_alignment>"

# Folders to process
target_folders=(
    "PRJNA1089479"
    "PRJNA371708"
    "PRJNA636791"
    "PRJNA752933"
    "PRJNA796810"
    "PRJNA996091"
)

# Iterate through the target folders
for target_folder in "${target_folders[@]}"; do
    # Create the output directory for the current target folder
    target_output_dir="$output_dir/$target_folder"
    mkdir -p "$target_output_dir"

    # Get the input subdirectory for the current target folder
    target_input_dir="$input_dir/$target_folder"

    # Iterate through the fastq.gz files in the current target folder
    for file_1 in "$target_input_dir"/*_1.fastq.gz; do
        # Get the file name without the extension
        file_base=$(basename "$file_1" _1.fastq.gz)

        # Get the corresponding fastq file for -2 option
        file_2="$target_input_dir/${file_base}_2.fastq.gz"

        # Run hisat2 command
        hisat2 -p 32 --dta -x "$genome_index" -1 "$file_1" -2 "$file_2" -S "$target_output_dir/${file_base}.sam" 2> "$target_output_dir/${file_base}.log"

        # Sort the output SAM file using samtools
        samtools sort -@ 32 -o "$target_output_dir/${file_base}_sorted.bam" "$target_output_dir/${file_base}.sam"

        # Remove the intermediate SAM file
        rm "$target_output_dir/${file_base}.sam"
    done
done
```

**Notes**:
- The `--dta` option optimizes for downstream transcript assembly.
- Output includes sorted BAM files and alignment logs.

### Aligning Single-End (SE) Short-Read Data

The following script aligns single-end FASTQ files to the *Artemisia* genome using HISAT2.

#### `hisat2_SE.sh`

```bash
#!/bin/bash
#SBATCH --job-name=hisat2_SE
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mypy

# Set the paths
genome_index="<path_to_genome>/LQ-9_phase0_genome_hisat2.idx"
input_dir="<path_to_filtered_fastq>"
output_dir="<path_to_alignment>"

# Folders to process
target_folders=(
    "PRJNA524494"
    "PRJNA755452"
)

# Iterate through the target folders
for target_folder in "${target_folders[@]}"; do
    # Create the output directory for the current target folder
    target_output_dir="$output_dir/$target_folder"
    mkdir -p "$target_output_dir"

    # Get the input subdirectory for the current target folder
    target_input_dir="$input_dir/$target_folder"

    # Iterate through the single-end fastq.gz files in the input directory
    for file in "${target_input_dir}"/*.fastq.gz; do
        # Get the file name without the extension
        file_base=$(basename "$file" .fastq.gz)

        # Run hisat2 command for single-end data
        hisat2 -p 32 --dta -x "${genome_index}" -U "${file}" -S "${target_output_dir}/${file_base}.sam" 2> "${target_output_dir}/${file_base}.log"

        # Sort the output SAM file using samtools
        samtools sort -@ 32 -o "${target_output_dir}/${file_base}_sorted.bam" "${target_output_dir}/${file_base}.sam"

        # Remove the intermediate SAM file
        rm "${target_output_dir}/${file_base}.sam"
    done
done
```

**Notes**:
- The `-U` option specifies single-end input.

### Filtering Alignments by Alignment Rate

The following script checks the overall alignment rate from HISAT2 logs and flags samples with rates < 50%.

#### `alignment_rate.sh`

```bash
#!/bin/bash

# Create a CSV file and write the header
echo "Project,Log file,Overall alignment rate,Results" > alignment_results.csv

# Set the paths
input_dir="<path_to_alignment>"

# Iterate through the subdirectories
for sub_dir in "$input_dir"/P*; do
    # Get the subdirectory name
    dir_name=$(basename "$sub_dir")
    subdirectory_path="${input_dir}/${dir_name}"

    # Iterate through the log files in the current subdirectory
    for log_file in "${subdirectory_path}"/*.log; do
        overall_alignment_rate=$(cat "$log_file" | grep -o -P '\d+\.\d+%(?=.*overall alignment rate)' | awk '{gsub(/%/, ""); print $1}')
        dir_name=$(basename "$subdirectory_path")
        log_file_name=$(basename "$log_file")
        # Determine the result based on the condition
        if (( $(bc <<< "$overall_alignment_rate <= 50") )); then
            result="Fail"
        else
            result="Pass"
        fi
        echo "$dir_name,$log_file_name,$overall_alignment_rate,$result" >> alignment_results.csv
    done
done
```

**Notes**:
- Requires `bc` for numerical comparisons (install via `conda install bc` if needed).
- The output CSV file (`alignment_results.csv`) lists projects, log files, alignment rates, and pass/fail status.

### Post-Alignment File Management

Samples with alignment rates < 50% were excluded from further analysis.

### Final Alignment Count

After filtering, the total number of alignments (sorted BAM files) is 266:

```bash
find <path_to_alignment> -name "*_sorted.bam" | wc -l
```

## Alignment of Long-Read RNA Sequencing Data Using minimap2

The following script aligns long-read FASTQ files to the *Artemisia* genome using minimap2.

#### `minimap2.sh`

```bash
#!/bin/bash
#SBATCH --job-name=minimap2
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mypy

# Set the paths
genome_index="<path_to_genome>/LQ-9_phase0_genome.fasta"
input_dir="<path_to_filtered_fastq>/longread"
output_dir="<path_to_alignment>/longread"

# Folders to process
target_folders=(
    "PRJNA752933"
)

# Iterate through the target folders
for target_folder in "${target_folders[@]}"; do
    # Create the output directory for the current target folder
    target_output_dir="$output_dir/$target_folder"
    mkdir -p "$target_output_dir"

    # Get the input subdirectory for the current target folder
    target_input_dir="$input_dir/$target_folder"

    # Iterate through the single-end fastq.gz files in the input directory
    for file in "${target_input_dir}"/*.fastq.gz; do
        # Get the file name without the extension
        file_base=$(basename "$file" .fastq.gz)

        # Run minimap2
        minimap2 -ax splice "${genome_index}" "${file}" > "${target_output_dir}/${file_base}.sam"

        # Sort the output SAM file using samtools
        samtools sort -@ 32 -o "${target_output_dir}/${file_base}_sorted.bam" "${target_output_dir}/${file_base}.sam"

        # Remove the intermediate SAM file
        rm "${target_output_dir}/${file_base}.sam"
    done
done
```

**Notes**:
- The `-ax splice` option optimizes for long-read RNA-seq alignment with splice awareness.

## Filtering and Transcriptome Assembly of RNA Sequencing Data
## Installation of Tools

Ensure the required tools are installed in your environment:

```bash
conda install samtools stringtie gffread
conda create -n isoquant
conda activate isoquant
conda install bioconda::isoquant
```

Activate the appropriate environment before running the scripts:

## Filtering Uniquely Mapped Reads Using samtools

Only uniquely mapped alignments are considered for further analysis. For paired-end (PE) data, alignments that are unmapped, have an unmapped mate, or are not primary alignments are filtered out (flag `-F 268`). For single-end (SE) data, alignments that are unmapped or not primary are filtered out (flag `-F 260`).

### Short-Read Data

#### `samtools_filter.sh`

```bash
#!/bin/bash
#SBATCH --job-name=samtools_filter
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mypy

# Set the paths
input_dir="<path_to_alignment>"
output_dir="<path_to_filtered_alignment>"

# Iterate through subdirectories starting with "PR"
for sub_dir in "$input_dir"/PR*; do
    # Get the subdirectory name
    dir_name=$(basename "$sub_dir")

    # Create the output directory for the current subdirectory
    subdirectory_output_dir="$output_dir/$dir_name"
    mkdir -p "$subdirectory_output_dir"

    # Set the appropriate samtools flag based on the subdirectory
    if [ "$dir_name" = "PRJNA755452" ]; then
        samtools_flag=260  # Single-end
        echo "Processing $dir_name as single-end data"
    else
        samtools_flag=268  # Paired-end
        echo "Processing $dir_name as paired-end data"
    fi

    # Iterate through the sorted.bam files in the current subdirectory
    for file in "$sub_dir"/*_sorted.bam; do
        # Get the file name without the extension
        file_base=$(basename "$file" _sorted.bam)

        # Define the output file path
        filtered_file="$subdirectory_output_dir/${file_base}_filtered.bam"

        # Run samtools with the appropriate flag
        samtools view -F "$samtools_flag" "$file" -O BAM -@ 64 > "$filtered_file"
    done
done
```

### Long-Read Data

#### `samtools_filter_long.sh`

```bash
#!/bin/bash
#SBATCH --job-name=samtools_filter
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mypy

# Set the paths
input_dir="<path_to_alignment>/longread"
output_dir="<path_to_filtered_alignment>/longread"

# Iterate through subdirectories starting with "P"
for sub_dir in "$input_dir"/P*; do
    # Get the subdirectory name
    dir_name=$(basename "$sub_dir")

    # Create the output directory for the current subdirectory
    subdirectory_output_dir="$output_dir/$dir_name"
    mkdir -p "$subdirectory_output_dir"

    # Set the appropriate samtools flag
    samtools_flag=260  # Single-end

    # Iterate through the sorted.bam files in the current subdirectory
    for file in "$sub_dir"/*_sorted.bam; do
        # Get the file name without the extension
        file_base=$(basename "$file" _sorted.bam)

        # Define the output file path
        filtered_file="$subdirectory_output_dir/${file_base}_filtered.bam"

        # Run samtools with the appropriate flag
        samtools view -F "$samtools_flag" "$file" -O BAM -@ 24 > "$filtered_file"
    done
done
```

**Notes**:
- The `-F 260` flag filters out unmapped reads and non-primary alignments for long-read single-end data.

## Merging Long-Read Alignments

### Generating BAM List

The following script generates a list of BAM files for merging long-read alignments.

#### `get_bam_list.sh`

```bash
#!/bin/bash

# Set the paths
input_dir="<path_to_filtered_alignment>/longread/PRJNA752933"
output_file="bam.list"

# Remove existing bam.list if it exists
rm -f "$output_file"

# Find all BAM files and print their absolute paths
find "$input_dir" -type f -name "*.bam" -print > "$output_file"
```

**Notes**:
- The output file `bam.list` contains the absolute paths to all BAM files in the specified directory.

### Merging Long-Read Alignments

The following script merges all long-read BAM files into a single sorted and indexed BAM file.

#### `merge_long_alignments.sh`

```bash
#!/bin/bash
#SBATCH --job-name=merge_long
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mypy

# Set the paths
dir="<path_to_filtered_alignment>/longread"

samtools merge -@ 32 -f -O BAM -b "$dir/bam.list" "$dir/all_long_alignments_merged.bam"
samtools sort -@ 32 "$dir/all_long_alignments_merged.bam" > "$dir/all_long_alignments_sortedByCoord.out.bam"
samtools index "$dir/all_long_alignments_sortedByCoord.out.bam"
```

**Notes**:
- The `-f` option forces overwriting if the output file exists.
- The final output is a sorted and indexed BAM file (`all_long_alignments_sortedByCoord.out.bam`).

## Transcriptome Assembly of Short-Read Sequences Using StringTie

The following script assembles transcriptomes from short-read alignments using StringTie.

#### `stringtie.sh`

```bash
#!/bin/bash
#SBATCH --job-name=stringtie
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mypy

# Set the paths
genome_annotation="<path_to_genome>/LQ-9_phase0_genome.gtf"
input_dir="<path_to_filtered_alignment>"
output_dir="<path_to_assembly>"

for sub_dir in "$input_dir"/P*; do
    # Get the subdirectory name
    dir_name=$(basename "$sub_dir")

    # Create the output directory for the current subdirectory
    subdirectory_output_dir="$output_dir/$dir_name"
    mkdir -p "$subdirectory_output_dir"

    # Iterate through the filtered.bam files in the current subdirectory
    for file in "${sub_dir}"/*_filtered.bam; do
        # Get the file name without the extension
        file_base=$(basename "$file" _filtered.bam)

        echo "Processing $file_base in $dir_name"

        # Run stringtie
        stringtie -p 32 -G "${genome_annotation}" -o "${subdirectory_output_dir}/${file_base}.gtf" "${input_dir}/${dir_name}/${file_base}_filtered.bam"
    done
done
```

## Merging Short-Read StringTie Assemblies

### Generating GTF List

The following script generates a list of GTF files for merging short-read assemblies.

#### `get_gtf_list.sh`

```bash
#!/bin/bash

# Variables
root_directory="<path_to_assembly>"
directory_pattern="PR*"
gft_pattern="*.gtf"
output_directory="<path_to_assembly>"
output_file="gtf_list.txt"

# Ensure the output directory exists
mkdir -p "$output_directory"

# Initialize the output file
> "$output_directory/$output_file"

# Validate the root directory
if [ ! -d "$root_directory" ]; then
    echo "Error: Root directory $root_directory does not exist."
    exit 1
fi

# Find and process GTF files in matching directories
find "$root_directory" -type f -path "*/$directory_pattern/$gft_pattern" -exec realpath {} \; >> "$output_directory/$output_file"

# Log completion
echo "GTF file paths written to: $output_directory/$output_file"
```

**Notes**:
- The output file `gtf_list.txt` contains the absolute paths to all GTF files.

### Running StringTie Merge

The following script merges short-read StringTie assemblies into a non-redundant set of transcripts.

#### `stringtie_merge.sh`

```bash
#!/bin/bash
#SBATCH --job-name=stringtie_merge
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mypy

# Set the paths
genome_annotation="<path_to_genome>/LQ-9_phase0_genome.gtf"
output_dir="<path_to_assembly>/stringtie_merge"

stringtie --merge -p 40 -c 5 -f 0.7 -m 200 -T 5 -G "${genome_annotation}" -o "${output_dir}/stringtie_merge.gtf" "${output_dir}/gtf_list.txt"
```

**Notes**:
- Replace `<path_to_genome>` and `<path_to_assembly>` with the appropriate paths.
- Parameters: `-c 5` (minimum coverage), `-f 0.7` (minimum isoform fraction), `-m 200` (minimum transcript length), `-T 5` (minimum TPM).

### Extracting Transcript Sequences (Short-Read)

Extract transcript sequences from the merged StringTie assembly:

```bash
gffread -w <path_to_assembly>/stringtie_merge/assembled_transcripts_stringtie.fa -g <path_to_genome>/LQ-9_phase0_genome.fasta <path_to_assembly>/stringtie_merge/stringtie_merge.gtf
```

**Notes**:
- Replace `<path_to_assembly>` and `<path_to_genome>` with the appropriate paths.
- The output is `assembled_transcripts_stringtie.fa`.

## Transcriptome Assembly of Long-Read Sequences Using StringTie

The following script assembles transcriptomes from long-read alignments using StringTie.

#### `stringtie_long.sh`

```bash
#!/bin/bash
#SBATCH --job-name=stringtie
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mypy

# Set the paths
genome_annotation="<path_to_genome>/LQ-9_phase0_genome.gtf"
input_dir="<path_to_filtered_alignment>/longread"
output_dir="<path_to_assembly>/longread"

stringtie -L -m 100 -p 32 -G "${genome_annotation}" -o "${output_dir}/stringtie_long_read.gtf" "${input_dir}/all_long_alignments_sortedByCoord.out.bam"
```

**Notes**:
- The `-L` option optimizes for long-read data, and `-m 100` sets the minimum transcript length.

### Extracting Transcript Sequences (Long-Read StringTie)

Extract transcript sequences from the long-read StringTie assembly:

```bash
gffread -w <path_to_assembly>/longread/assembled_transcripts_smrt_stringtie.fa -g <path_to_genome>/LQ-9_phase0_genome.fasta <path_to_assembly>/longread/stringtie_long_read.gtf
```

**Notes**:
- The output is `assembled_transcripts_smrt_stringtie.fa`.

## Transcriptome Assembly of Long-Read Sequences Using IsoQuant

The following script assembles transcriptomes from long-read alignments using IsoQuant.

#### `isoquant.sh`

```bash
#!/bin/bash
#SBATCH --job-name=isoquant
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=%j.out
#SBATCH --error=%j.err

# Increase the open file limit to prevent "Too many open files" errors
ulimit -n 8192

module load miniconda3
source activate isoquant

# Set the paths
genome_annotation="<path_to_genome>/LQ-9_phase0_genome.gtf"
genome="<path_to_genome>/LQ-9_phase0_genome.fasta"
input_dir="<path_to_filtered_alignment>/longread"
output_dir="<path_to_assembly>/longread"

isoquant.py --stranded none --threads 32 --reference "${genome}" --bam "${input_dir}/all_long_alignments_sortedByCoord.out.bam" --genedb "${genome_annotation}" --output "${output_dir}/IsoQuant" --check_canonical --data_type pacbio --count_exons
```

### Extracting Transcript Sequences (IsoQuant)

Extract transcript sequences from the IsoQuant assembly:

```bash
gffread -w <path_to_assembly>/longread/assembled_transcripts_isoseq.fa -g <path_to_genome>/LQ-9_phase0_genome.fasta <path_to_assembly>/longread/isoquant.gtf
```

## Indexing the Genome for Faster Processing

Create a genome index to improve the speed of `gffread`:

```bash
samtools faidx <path_to_genome>/LQ-9_phase0_genome.fasta
```

## Extracting Reference Transcriptome

Extract transcript sequences from the reference annotation:

```bash
gffread -w <path_to_assembly>/transcriptome.fa -g <path_to_genome>/LQ-9_phase0_genome.fasta <path_to_genome>/LQ-9_phase0_genome.gtf
```

## Transcriptome Quality Assessment Using BUSCO (v5.7.1)

## Installation of BUSCO

Install BUSCO and download the required lineage dataset:

```bash
conda create -n busco
conda activate busco
conda install -c conda-forge -c bioconda busco=5.7.1
busco --download embryophyta_odb10
```

## BUSCO Analysis for Reference Genome

Running BUSCO on the reference genome and its annotation provides insight into the effectiveness of the annotation method. The analysis is performed in genome mode (`-m genome`).

### `busco_genome.sh`

```bash
#!/bin/bash
#SBATCH --job-name=busco_genome
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate busco

# Set the paths
home_dir="<path_to_genome>"
buscodb="<path_to_busco_downloads>/embryophyta_odb10"

busco -f -i "${home_dir}/LQ-9_phase0_genome.fasta" -l "${buscodb}" -o busco_genome -m genome -c 24
```

## BUSCO Analysis for Transcriptome Assemblies

The following script evaluates the quality of multiple transcriptome assemblies (Reference transcriptome, StringTie, IsoSeq, and SMRT StringTie) in transcriptome mode (`-m tran`).

### `busco_slurm.sh`

```bash
#!/bin/bash

# Directories and resources
home_dir="<path_to_assembly>"
buscodb="<path_to_busco_downloads>/embryophyta_odb10"
output_dir="<path_to_busco_output>"

# List of assembly files
assemblies=(
    "assembled_transcripts_stringtie.fa"
    "transcriptome.fa"
    "assembled_transcripts_isoseq.fa"
    "assembled_transcripts_smrt_stringtie.fa"
)

# Loop through each assembly and create a separate SLURM job
for assembly in "${assemblies[@]}"; do
    # Extract base name for unique output
    base_name=$(basename "$assembly" .fa)
    output_name="busco_${base_name}"

    # Write a temporary SLURM script for the current assembly
    job_script="run_busco_${base_name}.sh"
    cat <<EOL > "$job_script"
#!/bin/bash
#SBATCH --job-name=busco_${base_name}
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --output=${output_dir}/${base_name}.out
#SBATCH --error=${output_dir}/${base_name}.err
module load miniconda3
source activate busco

busco -i "${home_dir}/${assembly}" \
      -l "$buscodb" \
      -o "${output_name}" \
      -m tran \
      -c 24 \
      -f

echo "BUSCO analysis completed for $assembly"
EOL

    # Submit the job to SLURM
    sbatch "$job_script"

    # Optional: Keep or remove the temporary script
    # rm "$job_script"
done

echo "All BUSCO jobs submitted."
```

### MMseqs2 Clustering

High numbers of duplicated genes in a BUSCO output often indicate redundancy in the transcriptome assembly. Clustering with a tool like MMseqs2 can help consolidate these sequences. We started with an assembly of **107,596 transcripts**.

#### Setting up the Environment

First, create a new `conda` environment and install **MMseqs2**.

```bash
conda create -n mmseqs2 -y
conda activate mmseqs2
conda install -c conda-forge -c bioconda mmseqs2 -y
```

#### Step-by-Step MMseqs2 Workflow

##### 1\. Create the MMseqs2 Database

The initial step is to create a database from the FASTA file of assembled transcripts. This is a prerequisite for all subsequent MMseqs2 commands.

```bash
#!/bin/bash
#SBATCH --job-name=mmseqs2_db
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mmseqs2

# Define input directory
input_dir="/path/to/your/assembly"

# Create a database from the FASTA file
mmseqs createdb "${input_dir}"/assembled_transcripts_stringtie.fa DB
```

##### 2\. Cluster the Transcripts

Next, we cluster the transcripts using `mmseqs easy-cluster`. A `min-seq-id` of `0.8` was chosen to cluster sequences with at least 80% identity. This step generates a representative FASTA file containing one sequence per cluster.

```bash
#!/bin/bash
#SBATCH --job-name=mmseqs2_cluster
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mmseqs2

# Define input directory
input_dir="/path/to/your/assembly"

# Perform clustering with 80% sequence identity threshold
mmseqs easy-cluster --min-seq-id 0.8 --threads 64 "${input_dir}"/assembled_transcripts_stringtie.fa clusterRes8 tmp
```

After clustering, the number of transcripts was reduced from 107,596 to **88,849**. This new, smaller set of representative sequences is stored in `clusterRes8_rep_seq.fasta`.

-----

### BUSCO Analysis of Clustered Transcripts

To evaluate the quality of the clustered assembly, we ran **BUSCO**. This provides a quantitative measure of gene content completeness.

```bash
#!/bin/bash
#SBATCH --job-name=busco_mmseqs2
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --output=%j.out
#SBATCH --error=%j.err

# Define directories and resources
input_dir="/path/to/your/assembly"
buscodb="/path/to/your/busco/db"
output_dir="/path/to/busco/output"
assembly_file="clusterRes8_rep_seq.fasta"

module load miniconda3
source activate busco

# Run BUSCO on the clustered assembly
busco -i "${input_dir}/${assembly_file}" \
      -l "$buscodb" \
      -o "busco_mmseqs2_clustered" \
      -m tran \
      -c 24 \
      -f

echo "BUSCO analysis completed for ${assembly_file}"
```

-----

### Filtering the GTF Annotation File

The final step for the clustering part is to create a new GTF file that only contains the annotations for the clustered transcripts.

```bash
# 1. Get the IDs of the clustered transcripts
grep ">" assembled_transcripts_stringtie_clustered.fa | sed 's/>//' > assembled_transcripts_stringtie_clustered.ids

# 2. Use gffread to filter the original GTF file
# Note: gffread can filter GFF3 files more directly, so we convert GTF to GFF3, filter, then convert back.
gffread stringtie_merge.gtf -O --ids assembled_transcripts_stringtie_clustered.ids > assembled_transcripts_stringtie_clustered.gff
gffread -T assembled_transcripts_stringtie_clustered.gff > assembled_transcripts_stringtie_clustered.gtf

# 3. Verify the transcript counts
echo "Original transcripts: 107596"
echo "Clustered transcripts:"
wc -l assembled_transcripts_stringtie_clustered.ids

echo "Transcripts in the filtered GTF file:"
awk -F "\t" '$3 == "transcript" {print}' assembled_transcripts_stringtie_clustered.gtf | wc -l
```


-----

### Portcullis Splice-Site Filtering

After clustering, we'll address potential false-positive splice sites using **Portcullis**. This is an important quality control step before merging annotations.

#### Setting up the Environment

Install `Portcullis` in a new `conda` environment.

```bash
conda create -n portcullis -y
conda activate portcullis
conda install bioconda::portcullis -y
```

#### Running Portcullis

The `portcullis full` command will identify and filter low-confidence splice junctions based on long-read alignments.

```bash
#!/bin/bash
#SBATCH --job-name=portcullis
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate portcullis

# Define directories
genome_dir="/path/to/your/genome"
bam_dir="/path/to/your/longread/alignment"
out_dir="/path/to/your/output"

# Run Portcullis
portcullis full -t 24 \
  "${genome_dir}"/LQ-9_phase0_genome.fasta \
  "${bam_dir}"/all_long_alignments_sortedByCoord.out.bam \
  -o "${out_dir}"
```

This command will generate a series of output files in the specified directory. The key output for downstream analysis is `portcullis_filtered.pass.junctions.bed`, which contains the high-confidence splice junctions.

Mikado is a pipeline designed to integrate multiple gene annotations into a single, high-quality, non-redundant gene set. This process involves several steps, from configuring the pipeline to running BLAST searches and predicting open reading frames (ORFs).

### 🛠️ Mikado Installation

First, install Mikado and its dependencies within a new Conda environment.

```bash
conda create -n mikado -y
conda activate mikado
git clone https://github.com/EI-CoreBioinformatics/mikado.git
cd mikado
pip3 install wheel
conda install -y -c conda-forge docutils==0.19
pip3 install -r requirements.txt
python3 setup.py bdist_wheel
pip3 install dist/*.whl
```

-----

### 📝 Step-by-Step Mikado Workflow

#### 1\. Create the Mikado Configuration File

A configuration file, `Mikado.conf`, is created to specify the input annotation files (GTF format) and their types (e.g., StringTie, IsoQuant).

```bash
# Define assembly directory
assembly_dir="/path/to/your/assembly"

# Create Mikado.conf with paths to your GTF files
echo -e "${assembly_dir}/assembled_transcripts_stringtie_clustered.gtf\tst\tFalse\t+1\tFalse\tTrue\tFalse" > Mikado.conf
echo -e "${assembly_dir}/stringtie_long_read.gtf\tstLR\tTrue\t\tFalse\tFalse\tFalse" >> Mikado.conf 
echo -e "${assembly_dir}/isoquant.gtf\tisoq\tTrue\t\tFalse\tFalse\tFalse" >> Mikado.conf 
echo -e "/path/to/your/reference_gtf/reference_genome.gtf\tref\tFalse\t\tFalse\tTrue\tFalse" >> Mikado.conf
```

-----

#### 2\. Run `mikado configure`

The `mikado configure` command generates a detailed `configuration.yaml` file that includes paths to the genome, splice junction files, and other parameters.

```bash
#!/bin/bash
#SBATCH --job-name=mikado
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --output=%j.out
#SBATCH --error=%j.err
module load miniconda3
source activate mikado

genome="/path/to/your/genome/LQ-9_phase0_genome.fasta"
input_dir="/path/to/your/annotations"
junction="/path/to/your/junctions/portcullis_filtered.pass.junctions.bed"
uniprot="/path/to/your/uniprot/uniprot_sprot.fasta"

mikado configure --list "${input_dir}"/Mikado.conf --scoring plant.yaml --reference "$genome" --codon-table 0 --junctions "$junction" -bt "$uniprot" configuration.yaml
```

-----

#### 3\. Run `mikado prepare`

This step parses the input GTF files based on the `configuration.yaml` and consolidates them into a single, non-redundant FASTA file named `mikado_prepared.fasta`.

```bash
#!/bin/bash
#SBATCH --job-name=mikado_prepare
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --output=%j.out
#SBATCH --error=%j.err
module load miniconda3
source activate mikado

input_dir="/path/to/your/annotations"

mikado prepare --json-conf "${input_dir}"/configuration.yaml --procs 8
```

-----

#### 4\. Blastx against Uniprot

The `mikado_prepared.fasta` file is then used as a query to search against the Uniprot database using **DIAMOND blastx**. This provides evidence for protein-coding genes.

```bash
#!/bin/bash
#SBATCH --job-name=blastx
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mikado

database="/path/to/your/uniprot/uniprot_sprot_new.dmnd"
input_dir="/path/to/your/annotations"

diamond blastx -d "${database}" -q "${input_dir}"/mikado_prepared.fasta --ultra-sensitive --threads "$SLURM_NTASKS_PER_NODE" --evalue 1e-5 --max-hsps 5 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop -o "${input_dir}"/mikado_prepared.blast.tsv
```

-----

#### 5\. ORF Calculation using TransDecoder

**TransDecoder** is used to identify potential protein-coding sequences within the prepared transcripts.

##### 5.1. Install TransDecoder

Install TransDecoder in the same `mikado` Conda environment.

```bash
conda activate mikado
conda install -y bioconda::transdecoder
```

##### 5.2. Prepare long ORFs

```bash
#!/bin/bash
#SBATCH --job-name=transdecoder_orfs
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=%j.out
#SBATCH --error=%j.err
module load miniconda3
source activate mikado

input_dir="/path/to/your/annotations"

TransDecoder.LongOrfs -t "${input_dir}"/mikado_prepared.fasta --output_dir "${input_dir}"/transdecoder
```

This step will create a directory named `transdecoder` containing the predicted ORFs.

##### 5.3. Blastp on Predicted ORFs

The predicted ORFs are then searched against the Uniprot database using **DIAMOND blastp** to find homologous proteins.

```bash
#!/bin/bash
#SBATCH --job-name=blastp
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mikado

database="/path/to/your/uniprot/uniprot_sprot_new.dmnd"
input_dir="/path/to/your/annotations"

diamond blastp -d "${database}" -q "${input_dir}"/transdecoder/mikado_prepared.fasta.transdecoder_dir/longest_orfs.pep --ultra-sensitive --threads "$SLURM_NTASKS_PER_NODE" -o "${input_dir}"/transdecoder/mikado_prepared.fasta.transdecoder_dir/longest_orfs.diamond.outfmt6 --outfmt 6 --evalue 1e-5 --max-target-seqs 1 
```

##### 5.4. Run `TransDecoder.Predict`

This final step uses the evidence from the BLASTp search to make a final prediction of the most likely coding regions within the transcripts.

```bash
#!/bin/bash
#SBATCH --job-name=transdecoder_predict
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=%j.out
#SBATCH --error=%j.err
module load miniconda3
source activate mikado

input_dir="/path/to/your/annotations"

TransDecoder.Predict -T 1000 -t "${input_dir}"/mikado_prepared.fasta --retain_blastp_hits "${input_dir}"/transdecoder/mikado_prepared.fasta.transdecoder_dir/longest_orfs.diamond.outfmt6
```

##### 5.5. Mikado `serialise`

The `mikado serialise` command creates a **SQLite database** (`mikado.db`) that stores all the information needed for Mikado's subsequent analysis. It integrates data from the configuration file, BLAST results, and ORF predictions.

```bash
#!/bin/bash
#SBATCH --job-name=serialise
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --output=%j.out
#SBATCH --error=%j.err
module load miniconda3
source activate mikado

# Define input directory
input_dir="/path/to/your/annotations"

mikado serialise -p "$SLURM_NTASKS_PER_NODE" \
  --json-conf "${input_dir}"/configuration.yaml \
  --xml "${input_dir}"/mikado_prepared.blast.tsv \
  --orfs "${input_dir}"/transdecoder/mikado_prepared.fasta.transdecoder_dir/longest_orfs.gff3
```

This command will create two files: `mikado.db` and `serialise.log`. You can inspect the database tables using a command like `sqlite3 mikado.db ".tables"`.

-----

##### 5.6. Mikado `pick`

The `mikado pick` command uses the SQLite database to select a final, high-quality set of transcripts from the input annotations. It generates the final annotation files in GFF3 format.

```bash
#!/bin/bash
#SBATCH --job-name=pick
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=%j.out
#SBATCH --error=%j.err
module load miniconda3
source activate mikado

# Define input directory
input_dir="/path/to/your/annotations"

mikado pick -p "$SLURM_NTASKS_PER_NODE" \
  --configuration "${input_dir}"/configuration.yaml \
  --subloci-out mikado.subloci.gff3 \
  --output-dir "${input_dir}"/mikado
```

#### Comparing different `pick` modes

Mikado offers different modes (e.g., `split`, `nosplit`, `stringent`, `lenient`, `permissive`) for gene model selection. To determine the best mode for your data, you can run a loop to generate a set of annotations for each mode and then compare them using a tool like BUSCO.

```bash
#!/bin/bash
#SBATCH --job-name=pick_loop
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mikado

# Define input directory
input_dir="/path/to/your/annotations"

# Loop through each mode and run mikado pick
for MODE in split nosplit stringent lenient permissive; do
  mikado pick --mode "$MODE" -p "$SLURM_NTASKS_PER_NODE" \
    --configuration "${input_dir}/configuration.yaml" \
    --subloci-out mikado.subloci.gff3 \
    --output-dir "${input_dir}/mikado.${MODE}"
done
```

The final chosen annotation set will be in the `mikado.loci.gff3` file within the selected output directory.

-----

##### 5.7. Finalizing the Annotation

After selecting the best annotation set, you must convert the GFF3 file to GTF format and extract the corresponding transcript sequences.

#### Convert GFF3 to GTF

The `mikado util convert` command is used to transform the GFF3 output into a GTF file.

```bash
# Convert the GFF3 file to GTF format
mikado util convert -of gtf -if gff3 mikado.loci.gff3 > mikado.loci.gtf
```

#### Extract Transcript Sequences

Use **gffread** to extract the transcript sequences from the final GTF file, using the reference genome as a guide.

```bash
# Extract sequences using gffread
gffread -w mikado.fasta -g /path/to/your/genome/LQ-9_phase0_genome.fasta mikado.loci.gtf
```

#### Count Transcript Groups and Isoforms

To understand the diversity of your final annotation, you can count the number of transcript groups (genes) and isoforms.

```bash
grep "^>" mikado.fasta | cut -d'.' -f1,2 | sort | uniq -c | sort -nr | head -1
```

### **Checking the Quality of Mikado Annotations using BUSCO**

To verify the quality and completeness of the Mikado-integrated annotation, we use **BUSCO** (Benchmarking Universal Single-Copy Orthologs). This analysis compares the newly generated transcripts against a known set of universal single-copy orthologs to assess the annotation's gene content.

-----

### **1. BUSCO Analysis of Mikado Output**

The following SLURM script runs BUSCO on the `mikado.fasta` file, which contains the final set of transcripts selected by Mikado. The analysis is performed in transcriptome mode (`-m tran`) using the `embryophyta_odb10` lineage database.

```bash
#!/bin/bash
#SBATCH --job-name=busco_mikado
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=%j.out
#SBATCH --error=%j.err
module load miniconda3
source activate busco

# Define directories
mikado_dir="/path/to/mikado/output"
buscodb="/path/to/your/busco_db/embryophyta_odb10"

# Run BUSCO
busco -i "${mikado_dir}"/mikado.fasta \
      -l "${buscodb}" \
      -o busco_mikado \
      -m tran \
      -c 32 \
      -f
```

This script will produce a summary report in the `busco_mikado` output directory, detailing the number of complete, fragmented, and missing BUSCO genes found in the Mikado annotation.

-----

### **2. Alternative Batch Submission Script**

If you need to run BUSCO on multiple files or prefer a more automated approach for submitting jobs to a SLURM cluster, this script can be used. It iterates through a list of assembly files and creates a dedicated SLURM script for each one, which is then submitted for execution.

```bash
#!/bin/bash

# Directories and resources
home_dir="/path/to/your/assembly/"
buscodb="/path/to/your/busco_db/embryophyta_odb10"
output_dir="/path/to/your/busco/output"

# List of assembly files
assemblies=(
    "mikado.fasta"
)

# Loop through each assembly and create a separate SLURM job
for assembly in "${assemblies[@]}"; do
    # Extract base name for unique output
    base_name=$(basename "$assembly" .fasta)
    output_name="busco_${base_name}"

    # Write a temporary SLURM script for the current assembly
    job_script="run_busco_${base_name}.sh"
    cat <<EOL > "$job_script"
#!/bin/bash
#SBATCH --job-name=busco_${base_name}
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --output=${output_dir}/${base_name}.out
#SBATCH --error=${output_dir}/${base_name}.err
module load miniconda3
source activate busco

busco -i "${home_dir}/${assembly}" \
      -l "${buscodb}" \
      -o "${output_name}" \
      -m tran \
      -c 24 \
      -f

echo "BUSCO analysis completed for $assembly"
EOL

    # Submit the job to SLURM
    sbatch "$job_script"

    # Optional: Keep or remove the temporary script
    # rm "$job_script"
done

echo "All BUSCO jobs submitted."
```

### **Estimating Transcript Abundances with Salmon**

This section outlines the process of using **Salmon** to estimate transcript abundances from RNA-seq data. It uses the mapping-based mode, which requires a transcriptome FASTA file (obtained from Mikado) and raw FASTQ files.

-----

### **1. Building the Salmon Index**

Before quantification, a Salmon index must be built. A crucial step for a high-quality index is to include **decoys** to handle reads that map to the genome but not to annotated transcripts. This helps prevent false-positive mappings to similar sequences.

#### **Generating the Decoy Transcriptome**

This script uses a custom script (`generateDecoyTranscriptome.sh`) to create a combined FASTA file (a "gentrome") from the genome and transcriptome, and a list of decoy sequences. It requires **mashmap** and **bedtools**.

```bash
#!/bin/bash
#SBATCH --job-name=decoy_transcriptome
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --output=%j.out
#SBATCH --error=%j.err
module load miniconda3
source activate salmon

# Define paths
main_directory="/path/to/your/scripts"
genome="/path/to/your/genome/LQ-9_phase0_genome.fasta"
transcriptome="/path/to/your/mikado/mikado.fasta"
annotation="/path/to/your/mikado/mikado.loci.gtf"
output="/path/to/your/salmon/output"

# Run the script to generate the gentrome and decoys
"${main_directory}"/generateDecoyTranscriptome.sh -j 64 -a "${annotation}" -g "${genome}" -t "${transcriptome}" -o "${output}"
```

This process generates two key files:

  - `gentrome.fa`: A FASTA file combining the transcriptome and the decoy sequences.
  - `decoys.txt`: A list of the decoy sequences.

-----

### **2. Transcript Quantification**

This pipeline script performs two main steps: building the Salmon index using the gentrome and decoys, and then quantifying transcript abundances for both single-end and paired-end reads.

```bash
#!/bin/bash
#SBATCH --job-name=salmon_pipeline
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --output=%j_pipeline.out
#SBATCH --error=%j_pipeline.err

# Load necessary modules and activate the environment
module load miniconda3
source activate salmon

# Define paths
salmon_directory="/path/to/your/salmon/bin"
gentrome="/path/to/your/salmon/output/gentrome.fa"
output="/path/to/your/salmon/output"
salmon_index="${output}/artemisia_index"
input_dir="/path/to/your/filtered_FASTQ"
output_dir="${output}/quants"
single_end_directory="PRJNA755452"

# Step 1: Create the Salmon index with decoys
echo "Starting Salmon Indexing..."
"${salmon_directory}"/salmon index -p 64 -t "${gentrome}" -i "${salmon_index}" --decoys "${output}/decoys.txt" -k 31
echo "Salmon Indexing Completed."

# Step 2: Quantify single-end reads
echo "Starting Salmon Quantification for Single-End Reads..."
single_end_output_dir="${output_dir}/${single_end_directory}"
mkdir -p "$single_end_output_dir"

for se_file in "${input_dir}/${single_end_directory}"/*.fastq.gz; do
    se_file_base=$(basename "$se_file" .fastq.gz)
    "${salmon_directory}"/salmon quant -i "$salmon_index" -l A -r "$se_file" -p 64 --validateMappings --dumpEq -o "${single_end_output_dir}/${se_file_base}_quant"
done
echo "Single-End Quantification Completed."

# Step 3: Quantify paired-end reads
echo "Starting Salmon Quantification for Paired-End Reads..."
for sub_dir in "$input_dir"/*; do
    dir_name=$(basename "$sub_dir")

    # Skip the single-end directory
    if [ "$dir_name" == "$single_end_directory" ]; then
        continue
    fi

    # Create output directory for paired-end reads
    subdirectory_output_dir="${output_dir}/${dir_name}"
    mkdir -p "$subdirectory_output_dir"

    # Process paired-end files
    for file_1 in "${sub_dir}"/*_1.fastq.gz; do
        file_base=$(basename "$file_1" _1.fastq.gz)
        file_2="${sub_dir}/${file_base}_2.fastq.gz"

        "${salmon_directory}"/salmon quant -i "${salmon_index}" -l A -1 "${file_1}" -2 "${file_2}" -p 64 --validateMappings --dumpEq -o "${subdirectory_output_dir}/${file_base}_quant"
    done
done
echo "Paired-End Quantification Completed."
```

-----

### **3. Checking Mapping Rates**

A crucial quality control step is to check the mapping rate for each sample. A mapping rate below 50% often suggests a problem with the sample or the reference transcriptome.

```bash
#!/bin/bash

# Set the main directory containing the quantification results
main_directory="/path/to/your/salmon/output/quants"
output_file="mapping_rates.csv"

# Remove the output file if it already exists
if [ -f "$output_file" ]; then
    rm "$output_file"
fi

# Add CSV header to the output file
echo "Run,mapping_rate,Status" >> "$output_file"

# Loop through project directories
for project_dir in "$main_directory"/*; do
    if [ -d "$project_dir" ]; then
        # Loop through read directories within each project directory
        for read_dir in "$project_dir"/*/; do
            if [ -d "$read_dir" ]; then
                # Locate the Salmon log file
                log_file="$read_dir/logs/salmon_quant.log"
                if [ -f "$log_file" ]; then
                    # Extract the mapping rate using grep and sed
                    rate=$(grep "Mapping rate" "$log_file" | sed -E 's/.*Mapping rate = ([0-9.]+)%.*/\1/')
                    if [ -z "$rate" ]; then
                        rate="NA"
                    fi

                    # Determine Pass/Fail status based on the mapping rate
                    if [[ "$rate" != "NA" && $(echo "$rate < 50" | bc -l) -eq 1 ]]; then
                        status="Fail"
                    else
                        status="Pass"
                    fi

                    # Extract the biosample name from the directory path
                    biosample=$(basename "$read_dir")

                    # Append the biosample, mapping rate, and status to the output file
                    echo "$biosample,$rate,$status" >> "$output_file"
                fi
            fi
        done
    fi
done

echo "Mapping rates extracted and saved to $output_file"
```

The script generates a CSV file named `mapping_rates.csv` which can be used to quickly identify samples with low mapping efficiency. In this case, all samples had a mapping rate of **\> 50%**.

### Functional Annotation of Mikado Transcripts

Functional analysis is a crucial step in understanding the biological roles of a newly assembled transcriptome. This section details the process of using **eggNOG-mapper** and **InterProScan** to annotate the transcripts from the Mikado assembly.

-----

### **1. eggNOG-mapper**

eggNOG-mapper is a tool for fast functional annotation of sequences based on the eggNOG database of orthologous groups.

#### **Installation and Setup**

First, create a dedicated Conda environment and clone the eggNOG-mapper repository.

```bash
conda create -n eggnogmapper -y
conda activate eggnogmapper
git clone https://github.com/eggnogdb/eggnog-mapper.git
cd eggnog-mapper
conda install --file requirements.txt -y
```

Next, download the necessary databases. We recommend downloading the HMMER-specific Viridiplantae database (`-H -d 33090`) as it was found to perform better than DIAMOND or MMseqs for this dataset.

```bash
python download_eggnog_data.py
python download_eggnog_data.py -P # Download Pfam database
python download_eggnog_data.py -H -d 33090 # Download Viridiplantae database for HMMER
```

To make the `emapper.py` command accessible from anywhere, add it to your system's `PATH`.

```bash
export PATH="/path/to/your/eggnog-mapper/eggnog-mapper:$PATH"
source ~/.bashrc # or source ~/.bash_profile
```

#### **Running eggNOG-mapper with HMMER**

The following SLURM script performs both the search and annotation steps in a single file. It uses HMMER as the search method and specifies the Viridiplantae database.

```bash
#!/bin/bash
#SBATCH --job-name=eggnogmapper_hmmer
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate eggnogmapper

# Define directories
mikado_fasta="/path/to/your/mikado.fasta"
output_dir="/path/to/your/eggnog-mapper/results"
temp_dir="/path/to/your/eggnog-mapper/temp_files"

# Run search step with HMMER
emapper.py --cpu 0 -m hmmer -d Viridiplantae --qtype seq --dbtype hmmdb --usemem --no_annot \
  -i "${mikado_fasta}" \
  --itype CDS -o hmmer --output_dir "${output_dir}" --evalue 1e-5 \
  --temp_dir "${temp_dir}"

# Run annotation step
emapper.py --cpu 0 -m no_search -i "${mikado_fasta}" --dbmem \
  --annotate_hits_table "${output_dir}/hmmer.emapper.seed_orthologs" \
  -o hmmer --output_dir "${output_dir}" --evalue 1e-5 \
  --target_taxa 33090 --pfam_realign none --temp_dir "${temp_dir}"
```

The final annotation results are stored in `hmmer.emapper.annotations` within the specified output directory.

-----

### **2. InterProScan**

InterProScan is a sequence analysis application that scans protein sequences against a library of protein family and domain databases.

#### **Installation and Setup**

InterProScan is typically installed by downloading the pre-compiled binaries and databases.

```bash
# Create environment and install dependencies
conda create -n interproscan -y
conda activate interproscan
conda install -y -c conda-forge openjdk

# Download InterProScan and database
# Use the provided FTP link to download the core software
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.67-99.0/interproscan-5.67-99.0-64-bit.tar.gz

# Check the checksum and extract the archive
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.67-99.0/interproscan-5.67-99.0-64-bit.tar.gz.md5
md5sum -c interproscan-5.67-99.0-64-bit.tar.gz.md5
tar xvzf interproscan-5.67-99.0-64-bit.tar.gz

# Prepare the HMM models
cd interproscan-5.67-99.0
python3 setup.py -f interproscan.properties
```

#### **Pre-processing Transcripts**

InterProScan works on protein sequences. Since the input from Mikado is a nucleotide FASTA file, it must be converted to protein sequences in six possible reading frames using `transeq` from the EMBOSS suite.

```bash
# Install EMBOSS
conda install -y bioconda::emboss

# Convert transcripts to protein sequences
transeq /path/to/your/mikado.fasta -frame=6 -clean -out mikado.proteins.fasta

# Split the large FASTA file into smaller chunks for parallel processing
seqkit split -s 80000 mikado.proteins.fasta
```

This will create a directory named `mikado.proteins.fasta.split` containing multiple FASTA files.

#### **Running InterProScan**

This SLURM script runs InterProScan on all the split FASTA files in parallel, specifying the Pfam and PANTHER databases.

```bash
#!/bin/bash
#SBATCH --job-name=interproscan
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --output=%j.out
#SBATCH --error=%j.err
module load miniconda3
source activate interproscan

# Define directories
interproscan_dir="/path/to/your/interproscan/interproscan-5.67-99.0"
input_dir="/path/to/your/mikado.proteins.fasta.split"
output_dir="/path/to/your/interproscan/output"

# Loop through each split file and run InterProScan
for file in "$input_dir"/*.fasta; do
    filename=$(basename "$file")
    "${interproscan_dir}"/interproscan.sh -cpu 64 -t p -dp -goterms \
      -appl Pfam,PANTHER \
      -i "$file" \
      --output-dir "$output_dir"
done
```

The individual `.tsv` output files from each run are then concatenated into a single file for downstream analysis.

```bash
cd "$output_dir"
cat *.tsv > mikado_interproscan.tsv
```

Functional annotation of transcripts involves using **DIAMOND blastx** to compare the assembled transcripts against established protein databases. This method translates the transcripts into six possible protein-coding frames and finds the most similar protein sequences in the database.

-----

### **BLASTx against Arabidopsis**

To find homologous genes in a closely related and well-studied plant, you can use the **Araport11** protein database. This SLURM script runs `diamond blastx` to find the best match for each of your transcripts in the Arabidopsis proteome.

```bash
#!/bin/bash
#SBATCH --job-name=blastx_arab
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate mikado

# Define paths for the database, query, and output
database="/path/to/your/araport11/Araport11_pep_20220914_new.dmnd"
query="/path/to/your/mikado.fasta"
output="/path/to/your/functional_annotation/arabidopsis_blastx"

# Run DIAMOND blastx with specified parameters
diamond blastx -d "${database}" -q "${query}" -k 1 --ultra-sensitive \
  --id 60 --max-target-seqs 1 --evalue 1e-5 --min-score 150 \
  -o "${output}"/transcriptome_vs_arabidopsis.outfmt6
```

The output file, `transcriptome_vs_arabidopsis.outfmt6`, contains the top hit for each transcript, providing a direct link to a known Arabidopsis gene and its annotation.

-----

### **BLASTx against Uniprot**

For a broader functional annotation, you can search your transcripts against the **UniprotKB/Swiss-Prot** database, which contains a high-quality, expertly curated set of protein sequences.

```bash
#!/bin/bash
#SBATCH --job-name=blastx_uni
#SBATCH --partition=64c512g
#SBATCH -N 1