#!/bin/bash/

echo My working directory is `pwd`
echo Running job on host:
echo -e '\t'`hostname` at `date`
echo

# purge any loaded modules
module purge

# change directories to /users/mfk513/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/ so that all the analysis is stored there.

cd ~/scratch/LsRTD_v1.0_valdn/Lsat_Bcin_timeseries/

# initialize a variable with an intuitive name to store the name of the input fastq file

fq=$1
idx=$2

# setup pathway to be analyse
pathway=${fq%_1.fastq.gz}

# grap base of input filename for naming outputs
samplename=$(echo $pathway | rev | cut -d '/' -f1 | rev)
idxname=$(echo $idx | rev | cut -d '/' -f1 | rev)
R1_untrimmed_pathway=${pathway}_1.fastq.gz
R2_untrimmed_pathway=${pathway}_2.fastq.gz
R1_untrimmed=$(echo $R1_untrimmed_pathway | rev | cut -d '/' -f1 | rev)
R2_untrimmed=$(echo $R2_untrimmed_pathway | rev | cut -d '/' -f1 | rev)

echo
echo Sample name is $samplename
echo Transcriptome index file is $idxname
echo File containing the 1th mates $R1_untrimmed
echo File containing the 2nd mates $R2_untrimmed

# specify the number of cores to use
thread=12

# make all of the output directories
# The -p option means mkdir will create the whole path if it 
# does not exist and refrain from complaining if it does exist

mkdir -p results/01_fastp/
mkdir -p results/02_fastqc/
mkdir -p results/03_fastq_trimmed/
mkdir -p results/04_salmon_quant/

# set up output filenames and locations
fastp_out=results/01_fastp/
fastqc_trimmed_out=results/02_fastqc/
fastq_trimmed_out=results/03_fastq_trimmed/
salmon_out=results/04_salmon_quant/

# set up the software environment
module load lang/Java/13.0.2
module load bio/fastp/0.23.2-GCC-11.3.0
module load bio/FastQC/0.11.7-Java-1.8.0_181
module load bio/Salmon/1.10.0-GCC-11.3.0
unset DISPLAY

echo
echo Processing file ${samplename}

echo
echo `date`
echo "Starting fastp for ${samplename}"

# Run fastp on input file and move output to the appropriate folder

fastp \
--in1 ${R1_untrimmed_pathway} \
--in2 ${R2_untrimmed_pathway} \
--detect_adapter_for_pe \
--html "${fastp_out}""${samplename}"_trim-report_fastp.html \
--json "${fastp_out}""${samplename}"_trim-report_fastp.json \
--out1 "${fastq_trimmed_out}""${R1_untrimmed%.fastq.gz}"_trimmed.fq.gz \
--out2 "${fastq_trimmed_out}""${R2_untrimmed%.fastq.gz}"_trimmed.fq.gz \

echo
echo `date`
echo fastp preprocessing on ${samplename} were completed

# Grap base of input filename for naming output
R1_trimmed_pathway="${fastq_trimmed_out}""${R1_untrimmed%.fastq.gz}"_trimmed.fq.gz
R2_trimmed_pathway="${fastq_trimmed_out}""${R2_untrimmed%.fastq.gz}"_trimmed.fq.gz
R1_trimmed=$(echo $R1_trimmed_pathway | rev | cut -d '/' -f1 | rev)
R2_trimmed=$(echo $R2_trimmed_pathway | rev | cut -d '/' -f1 | rev)

echo
echo `date`
echo Starting FastQC for trimmed ${samplename}

# Run FastQC and move output to the appropriate folder
fastqc -o ${fastqc_trimmed_out} ${R1_trimmed_pathway}
echo
fastqc -o ${fastqc_trimmed_out} ${R2_trimmed_pathway}

echo
echo `date`
echo FastQC analysis were completed 

echo
echo `date`
echo Starting `salmon -v` run for ${samplename} using ${idxname} transcriptome indices

# Run Salmon and move output to the appropriate folder

salmon quant -i ${idx} \
-p ${thread} \
-l A \
-1 ${R1_trimmed_pathway} \
-2 ${R2_trimmed_pathway} \
-o "${salmon_out}"/Quant_"${idxname}"/${samplename} \
--seqBias \
--posBias \
--gcBias \
--numBootstraps 100 \
--validateMappings

echo
echo `date`
echo Salmon mapping-based alignment were completed for ${samplename} using ${idxname} transcriptome indices

echo
echo Job completed at `date`


# This script takes a fastq file of RNA-seq data, runs fastp, FastQC, and Salmon.
# USAGE: sh rnaseq_on_input_file.sh <name of fastq file> <directory of transcriptome index file>
