#!/bin/bash

# purge any loaded modules
module purge

# change directories to ~/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/ so that all the analysis is stored there.
cd ~/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/

# initialize a variable with an intuitive name to store the name of the input fastq file
fq=$1

# grab base of filename for naming outputs

samplename=`basename $fq .fq.gz`
echo "Sample name is $samplename"  

# specify the number of cores to use

cores=6

# make the output directory
# The -p option means mkdir will create the whole path if it 
# does not exist and refrain from complaining if it does exist

mkdir -p results/02_FastQC/

# set up output filenames and locations

fastqc_out=results/02_FastQC/

# set up the software environment (use version numbers)

module load bio/FastQC/0.11.7-Java-1.8.0_181
unset DISPLAY

echo "Processing file $fq"
echo "Starting QC for $samplename AFTER fastp quality filtering"

# Run FastQC and move output to the appropriate folder

fastqc -o $fastqc_out $fq

# This script takes an trimmed fastq file of paired fastq data in /results/02_trimmomatic folder, runs FastQC.
# USAGE: sh FastQC_on_trimmed_input_file.sh <name of fastq file>
