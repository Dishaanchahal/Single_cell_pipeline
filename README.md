# Intro
This pipeline demonstrates the complete workflow for analyzing single-cell RNA sequencing data from ovarian cancer cells that developed resistance to carboplatin chemotherapy. The analysis progresses from raw FASTQ files through quality control, genome alignment with STAR Solo, and downstream analysis with Scanpy to identify distinct cell populations and expression patterns.

## Step 1
For the first step I dowloaded the data from the SRA database from NCBI's website. This data was the single cell RNA sequencing data from ovarian cancer cell lines that became resistant to the chemotherapy drug carboplatin.It was done by UC San Diego as part of their larger project, Evolution of platinum resistance. I chose this dataset because it was readily available and I wanted to work on a dataset related to cancer. Also this might be interesting because all cells will not develop cancer in the same way, some will be more resistant to chemotherapy with carboplatin and some less.

Sequencing machine used for this experiment was HiSeq4000 

The raw sra file was about 1.94 GB = SRR6334436. This is the run accession number of the file.

After I downloaded this into my computer I ran fastq on this which gave me two files - fastq_1.fastq and fastq_2.fastq

```bash
- fasterq-dump SRR6334436 --split-files -O ./fastq/
```

Output files
- fastq_1.fastq(26 bp reads) -  COntains cell barcodes + UMI sequences
- fastq_2.fastq(98 bp reads) - Contains actual cDNA.

This can be confusing first because we thought that we were sequencing RNA, but we start with RNA, then do reverse transcription on it to make complimentary cDNA and then send it for PCR amplification to make more compies. We use DNA becuse it is more stable for sequencing reactions. RNA is unstable and degrades quickly. 

After running this, This is the output I got
spots read      : 28,768,375
reads read      : 57,536,750
reads written   : 57,536,750

## Step 2
Now after we have our fastq files, we will perform quality control on those by running 

```bash
fastqc - fastqc fastq_1.fastq fastq_2.fastq -o results/    
'''

This will output an html and zip file ( attached in the repo).

html file is the interactive quality report for reads and zip is the raw data. 


## Step 3: Reference Genome Preparation

After confirming FASTQ file quality, we prepared the reference genome for alignment using STAR.

### Reference Files Used:
- **Genome sequence**: `GRCh38.primary_assembly.genome.fa` - Human reference genome (GRCh38/hg38)
- **Gene annotations**: `gencode.v47.annotation.gtf` - Gene Transfer Format file containing gene locations, exon boundaries, transcript isoforms, and UTR regions

### STAR Genome Index Generation:

```bash
STAR --runThreadN 20 \
     --runMode genomeGenerate \
     --genomeDir genome_index \
     --genomeFastaFiles ref/GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile ref/gencode.v47.annotation.gtf \
     --sjdbOverhang 99
```

**Parameters explained:**
- `--runThreadN 20`: Use 20 CPU threads for faster processing
- `--runMode genomeGenerate`: Create genome index for alignment
- `--genomeDir genome_index`: Output directory for index files
- `--genomeFastaFiles`: Input reference genome sequence
- `--sjdbGTFfile`: Gene annotation file for splice junction detection
- `--sjdbOverhang 99`: Optimal for 100bp reads (read length - 1)

This step creates an indexed version of the human genome that STAR can efficiently search during alignment, incorporating known splice junction information for accurate RNA-seq mapping.

## Step 4: Single-Cell RNA-seq Quantification with STAR Solo

With the genome index prepared, we used STAR Solo to simultaneously align reads and quantify gene expression at the single-cell level.

### Initial STAR Solo Command:
```bash
STAR --runThreadN 20 \
     --genomeDir genome_index \
     --readFilesIn fastq_1.fastq fastq_2.fastq \
     --outFileNamePrefix results/STARsolo/sample1_ \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 80000000000 \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 16 \
     --soloUMIstart 17 --soloUMIlen 10 \
     --soloFeatures Gene \
     --soloCellFilter None \
     --soloCBwhitelist None \
     --soloBarcodeReadLength 0
```

### Parameters Explained:
- `--runThreadN 20`: Use 20 CPU threads for parallel processing
- `--readFilesIn`: Input FASTQ files (R1: barcodes, R2: cDNA)
- `--soloType CB_UMI_Simple`: Process cell barcodes and UMIs for single-cell data
- `--soloCBstart/len`: Extract 16bp cell barcodes starting at position 1
- `--soloUMIstart/len`: Extract 10bp UMIs starting at position 17
- `--soloCellFilter None`: No initial cell filtering applied
- `--soloCBwhitelist None`: No barcode whitelist used (initial approach)

**Note**: This initial run was wrong, it produced suboptimal results (only 9 cells detected) due to missing barcode validation. The approach was later refined with proper 10x barcode whitelists and corrected file ordering.

To correct these mistakes, Here's the optimal STAR Solo command that worked:

## Step 4 (Revised): Optimized Single-Cell Quantification

### Why the Initial Run Failed:
The first attempt detected only 9 cells due to two critical issues:
1. **Incorrect file order**: STAR Solo expected cDNA reads first, then barcode reads
2. **Missing barcode whitelist**: Without validated 10x barcodes, most cellular barcodes were rejected as invalid

### The Need for Barcode Whitelists:
10x Genomics uses precisely designed barcode sequences to minimize sequencing errors. The whitelist contains 737,280 validated barcodes that:
- Are spaced apart to avoid confusion from sequencing errors
- Enable distinction between real cells and empty droplets
- Filter out random sequences that aren't genuine cell identifiers

### Optimal STAR Solo Command:
```bash
# Download 10x v2 barcode whitelist
curl -O https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/737K-april-2014_rc.txt

# Corrected STAR Solo run
STAR --runThreadN 20 \
     --genomeDir genome_index \
     --readFilesIn fastq_2.fastq fastq_1.fastq \  # Swapped file order
     --outFileNamePrefix results/STARsolo/sample1_swapped_ \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 80000000000 \
     --soloType CB_UMI_Simple \
     --soloCBstart 2 --soloCBlen 14 \              # Adjusted for actual barcode structure
     --soloUMIstart 16 --soloUMIlen 11 \
     --soloFeatures Gene \
     --soloCellFilter CellRanger2.2 \              # Enabled filtering
     --soloCBwhitelist 737K-april-2014_rc.txt \    # Added whitelist
     --soloBarcodeReadLength 0
```

### Results:
- **Raw cells detected**: 306,301 (vs. 9 in initial run)
- **Valid barcodes**: 100% (vs. 5.48% initially)
- **Reads mapped to genes**: 58.8% (vs. 0.05% initially)


## Step 5: Scanpy

After successful quantification with STAR Solo, we performed downstream analysis using Scanpy to identify cell populations and expression patterns.

### Analysis Pipeline:

```python
# Load STAR Solo output
adata = sc.read_10x_mtx("results2/STARsolo/sample1_swapped_Solo.out/Gene/raw/")

# Quality control and filtering
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"])
sc.pp.filter_cells(adata, min_genes=200)  # Remove low-quality cells
sc.pp.filter_genes(adata, min_cells=3)   # Remove rarely expressed genes

# Normalization and feature selection
sc.pp.normalize_total(adata, target_sum=1e4)  # Normalize to 10,000 counts/cell
sc.pp.log1p(adata)  # Log transform
sc.pp.highly_variable_genes(adata)  # Identify variable genes

# Dimensionality reduction and clustering
sc.tl.pca(adata)      # Principal component analysis
sc.pp.neighbors(adata) # Build neighborhood graph
sc.tl.umap(adata)     # UMAP visualization
sc.tl.leiden(adata)   # Leiden clustering
```

### Key Outputs:
- Quality control plots showing cell and gene metrics
- UMAP visualization revealing cell clusters
- Leiden clustering identifying distinct cell populations
- All figures automatically saved to `results/figures/`

This analysis transforms raw count data into interpretable cell clusters, enabling identification of different cell states within the platinum-resistant ovarian cancer sample.