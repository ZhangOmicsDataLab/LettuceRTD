Lettuce RTD Benchmarking Workflow

**Note:** This workflow was executed using `shell` scripts on the [Viking HPC](https://www.york.ac.uk/it-services/services/viking-computing-cluster/) facility at the University of York.

---

## Data Processing

To ensure the strandedness of RNA-Seq data, we utilized the [`how-are-we-stranded-here v1.0.1`](https://github.com/signalbash/how_are_we_stranded_here) Python package. At least 5 pairs of fastq files were checked with the following command:

```bash
check_strandedness --gtf GCF_002870075.3_Lsat_Salinas_v8_genomic.gtf --transcripts GCF_002870075.3_Lsat_Salinas_v8_cds_from_genomic.fna  --reads_1 *_1.fastq.gz --reads_2 *_2.fastq.gz
```
### Merge raw fastq files
Raw fastq files from sequencing replicates were merged using the following `bash` commands and metadata table:
```bash
cut -f1 metadata/metadata.tsv | \
	sort -k1,1 -u | \
	while read p; do
		grep -w "${p}" metadata/metadata.tsv | \
		awk '{print "rawdata/"$8"_1.fastq.gz"}' | \
		xargs cat > data/${p}_1.fastq.gz
		done < /dev/stdin
```

```bash
cut -f1 metadata/metadata.tsv | \
	sort -k1,1 -u | \
	while read p; do
		grep -w "${p}" metadata/metadata.tsv | \
		awk '{print "rawdata/"$8"_2.fastq.gz"}' | \
		xargs cat > data/${p}_2.fastq.gz
		done < /dev/stdin
```

### Reference Data
#### Lactuca sativa (garden lettuce)
For lettuce, along with LsRTDv1, [Lsat_Salinas_v11](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_002870075.4/) GenBank (GCA_002870075.4) and RefSeq(GCF_002870075.4) accessions were used as reference genome assemblies and annotations. The genome data sets were download as [describe](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#downloadservice) using the following rsync commands:

```bash
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/870/075/GCF_002870075.4_Lsat_Salinas_v11/ GCF_002870075.4_Lsat_Salinas_v11/
```

```bash
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/870/075/GCA_002870075.4_Lsat_Salinas_v11/ GCA_002870075.4_Lsat_Salinas_v11/
```

#### *Botrytis cinerea* B05.10
For *Botrytis cinerea*, [ASM14353v4](https://www.ncbi.nlm.nih.gov/labs/data-hub/genome/GCF_000143535.2/) GenBank (GCA_000143535.2) accession were used as reference assembly and annotation. The genome data sets were download as [describe](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#downloadservice) using the following rsync commands:

```bash
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/143/535/GCA_000143535.4_ASM14353v4/ GCA_000143535.4_ASM14353v4/
```

#### *Sclerotinia sclerotiorum* 1980 UF-70
For *Sclerotinia sclerotiorum*, [ASM185786v1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_001857865.1/) GenBank (GCA_001857865.1) accession were used as reference assembly and annotation. The genome data sets were download as [described](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#downloadservice) using the following rsync commands:

```bash
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/857/865/GCA_001857865.1_ASM185786v1/ GCA_001857865.1_ASM185786v1/
```

### Preprocessing
#### fastp quality control and trimming
We performed quality control and trimming on raw read files using [`fastp v0.23.2`](https://github.com/OpenGene/fastp) with the following command. The adapter sequence auto-detection is enabled via `--detect_adapter_for_pe`. The default settings were preserved for the other configuration options.

```bash
fastp \
--in1 ${R1_input_path} \
--in2 ${R2_input_path} \
--detect_adapter_for_pe \
--html "${output_dir}""01_fastp/""${sample_name}"_trim-report_fastp.html \
--json "${output_dir}""01_fastp/""${sample_name}"_trim-report_fastp.json \
--out1 "${output_dir}""03_fastq_trimmed/""${R1_sample_name%.fastq.gz}"_trimmed.fq.gz \
--out2 "${output_dir}""03_fastq_trimmed/""${R2_sample_name%.fastq.gz}"_trimmed.fq.gz \
```

#### FastQC quality control after trimming
FastQC quality control analysis was performed on trimmed read files after fastp using [`FastQC v0.11.7`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) with the following command:

```bash
fastqc *.fq.gz
```

### Read Quantification
All read quantification steps were carried out using [`Salmon v1.10.0`](https://github.com/COMBINE-lab/salmon) as described in its [documentation](https://salmon.readthedocs.io/en/latest/).

Salmon cDNA-only index files were generated for the combined transcriptome of lettuce (Lactuca sativa) and pathogens. A joint lettuce and pathogen (Botrytis cinerea and Sclerotinia sclerotiorum) transcriptome was created by concatenating the transcriptome files:

```bash
cat ./file1.fna.gz ./file2.fna.gz > file_merged.fna.gz
```

Salmon index files were generated using the `*_rna_from_genomic.fna.gz` file format.

```bash
salmon index -t *_rna_from_genomic.fna.gz -p 12 -i *_salmon_index
```

Salmon alignment and read quantification were performed in mapping-based mode using the combined transcript indices:

```bash
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
```

#### Summarize QC results using MultiQC
QC results from fastp, FastQC, and Salmon were aggregated into a single report using [`MultiQC v1.13`](https://github.com/MultiQC/MultiQC):

```bash
# Run MultiQC and move output to the appropriate folder
multiqc  results/ \
--outdir ${multiqc_out} \
-f --config scripts/multiqc_config.yaml \
```

### DE, DAS and DTU analysis
As a downstream analysis, we utilized the [3D-RNAseq App](https://3drnaseq.hutton.ac.uk/app_direct/3DRNAseq/) to perform differential expression (DE), differential alternative splicing (DAS), and differential transcript usage (DTU) on the Salmon quantification file (`quant.sf`).