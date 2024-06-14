#!/bin/bash

# purge any loaded modules
module purge

# change directories to ~/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/ so that all the analysis is stored there.
cd ~/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/

# make the output directories
# The -p option means mkdir will create the whole path if it
# does not exist and refrain from complaining if it does exist

mkdir -p results/01_fastp/
mkdir -p results/02_FastQC/
mkdir -p results/03_fastq_trimmed/
mkdir -p results/04_Kallisto/
mkdir -p results/05_Multiqc/

# initialize a variable with an intuitive name to store the name of the input fastq file
fq=$1

# Setup input file pathways
input_path=${fq%_1.fastq.gz}
R1_input_path=${input_path}_1.fastq.gz
R2_input_path=${input_path}_2.fastq.gz

# Setup output file names
output_dir=~/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/results/
sample_name=$(echo $input_path | rev | cut -d '/' -f1 | rev)
R1_sample_name=$(echo $R1_input_path | rev | cut -d '/' -f1 | rev)
R2_sample_name=$(echo $R2_input_path | rev | cut -d '/' -f1 | rev)
fastqc_trimmed_out=results/02_FastQC/
trimmed_checksums=trimmed_fastq_checksums.md5

# set up the software environment (use version numbers)
module load bio/fastp/0.23.2-GCC-10.3.0
module load bio/FastQC/0.11.7-Java-1.8.0_181
unset DISPLAY

# Write pathways and files to be trimmed
echo "Procession file path: $input_path"
echo "Starting fastp on $sample_name"
echo "Forward input file: $R1_sample_name"
echo "Reverse input file: $R2_sample_name"

# Run fastp on input file and move output to the appropriate folder
fastp \
--in1 ${R1_input_path} \
--in2 ${R2_input_path} \
--detect_adapter_for_pe \
--html "${output_dir}""01_fastp/""${sample_name}"_trim-report_fastp.html \
--json "${output_dir}""01_fastp/""${sample_name}"_trim-report_fastp.json \
--out1 "${output_dir}""03_fastq_trimmed/""${R1_sample_name%.fastq.gz}"_trimmed.fq.gz \
--out2 "${output_dir}""03_fastq_trimmed/""${R2_sample_name%.fastq.gz}"_trimmed.fq.gz \

# This script takes a fastq file of RNA-seq data in raw_data folder and runs fastp.
# USAGE: sh fastp_on_input_file.sh <name of fastq file>
