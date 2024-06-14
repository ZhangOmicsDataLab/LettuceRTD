#!/bin/bash
#SBATCH --job-name=multiqc               # Job name
#SBATCH --mail-type=END,FAIL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mfk513@york.ac.uk    # Where to send mail  
#SBATCH --ntasks=1                       # Run on a single CPU
#SBATCH --mem=1gb                        # Job memory request
#SBATCH --time=01:00:00                  # Time limit hrs:min:sec
#SBATCH --output=logs/%j.log             # Standard output and error log
#SBATCH --account=biol-gensac-2019       # Project account


# purge any loaded modules
module purge

# change directories to ~/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/ so that all the analysis is stored there.
cd ~/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/

# make the output directories
# The -p option means mkdir will create the whole path if it
# does not exist and refrain from complaining if it does exist

mkdir -p results/05_multiqc/
 
# setup output file directory
multiqc_out=results/05_multiqc/

# set up the software environment (use version numbers)
module load bio/MultiQC/1.13-foss-2021b
unset DISPLAY

# Write pathways and files to be quantified for transcript abundance via Kallisto
echo "Starting MultiQC for all stdout files from fastp and FastQC"

# Run MultiQC and move output to the appropriate folder
multiqc -d  results/04_salmon_quant/ \
--outdir ${multiqc_out} --filename MultiQC_salmon_report \
-f --config scripts/multiqc_salmon_config.yaml \

# This script takes a trimmed paired fastq files in results/02_trimmomatic folder and runs Kallisto quant.
# USAGE: sh MultiQC_on_all_stdout.sh

