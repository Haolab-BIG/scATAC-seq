# scATAC-seq Quality Control and Data Cleaning Pipeline

This repository provides a comprehensive workflow for single-cell ATAC-seq (scATAC-seq) data cleaning and quality control (QC), with support for **10x Genomics** data using **Cell Ranger ATAC** and **Signac** in **R**.

---

## 📌 Pipeline Overview

The pipeline consists of the following stages:

1. **Primary Processing**: FASTQ → Filtered peak–barcode matrix  
2. **Raw Data Quality Control**: Sequencing, cell-level, library, mapping, and targeting metrics  
3. **Signac-Based Data Cleaning**: Creating Seurat objects and applying filtering thresholds

---

## 🧬 Part I: Primary Processing (`FASTQ` → `filtered_peak_bc_matrix/`)

We use **Cell Ranger ATAC** (by 10x Genomics) to process raw sequencing data into a usable matrix for downstream analysis.

### 🔧 Required Inputs

| Option         | Description                        |
|----------------|------------------------------------|
| `--id`         | Output directory name              |
| `--reference`  | Path to reference genome           |
| `--fastqs`     | Path to FASTQ files                |
| `--sample`     | FASTQ filename prefix              |
| `--localcores` | Number of CPU threads              |
| `--localmem`   | RAM in GB                          |

### ✅ Example Command

```bash
cellranger-atac count \
  --id 10k_pbmc \
  --reference refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
  --fastqs 10k-human-pbmcs-atac \
  --sample atac_pbmc_10k_nextgem \
  --localcores 64 \
  --localmem 128
```

### 📁 Output Structure

```
filtered_peak_bc_matrix/
├── barcodes.tsv
├── matrix.mtx
└── peaks.bed
```

---

## 🧪 Part II: Raw Data Quality Control

### 2.1 Sequencing Metrics

| Metric                         | Threshold     |
|-------------------------------|---------------|
| Valid barcodes                | > 85%         |
| Q30 bases in barcode (i2)     | > 65%         |
| Q30 bases in Read 1 & Read 2  | > 65%         |
| Q30 in sample index (i1)      | > 90%         |

### 2.2 Cell-Level Metrics

| Metric                                | Threshold        |
|---------------------------------------|------------------|
| Estimated number of cells             | 500–10,000 ±20%  |
| Mean raw read pairs per cell          | > 5,000          |
| High-quality fragments in cells       | > 40%            |
| Transposition events in peaks         | > 25%            |
| Median HQ fragments per cell          | > 100            |

### 2.3 Library Complexity

| Metric           | Threshold     |
|------------------|---------------|
| Percent duplicates | Low expected |

### 2.4 Mapping Metrics

| Metric                           | Threshold   |
|----------------------------------|-------------|
| Confidently mapped read pairs    | > 80%       |
| Unmapped read pairs              | < 5%        |
| Non-nuclear read pairs           | < 20%       |
| Fragments in nucleosome-free regions | > 40%   |

### 2.5 Targeting Metrics

| Metric                                         | Threshold |
|------------------------------------------------|-----------|
| Fraction of genome in peaks                    | < 75%     |
| TSS enrichment score                           | > 5       |
| HQ fragments overlapping TSS                   | > 25%     |
| HQ fragments overlapping peaks                 | > 25%     |

---

## 🧰 Part III: Signac-Based QC and Filtering (R Pipeline)

We use **Signac** and **Seurat** to perform downstream QC and analysis.

### 📦 Required Packages

```r
library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
```

### 🔹 Step 1: Create Seurat Object

```r
counts <- Read10X_h5("filtered_peak_bc_matrix.h5")
metadata <- read.csv("singlecell.csv", row.names = 1)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  fragments = "fragments.tsv.gz",
  sep = c(":", "-"),
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
```

### 🔹 Step 2: Add Genome Annotation

```r
ah <- AnnotationHub()
ensdb <- ah[["AH75011"]]  # Ensembl v98 (hg38)
annotations <- GetGRangesFromEnsDb(ensdb)
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(pbmc) <- annotations

# Filter to standard chromosomes
keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc <- subset(pbmc, features = keep)
```

### 🔹 Step 3: Compute QC Metrics

```r
pbmc <- NucleosomeSignal(object = pbmc)
pbmc <- TSSEnrichment(object = pbmc)

pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

pbmc$blacklist_ratio <- FractionCountsInRegion(
  object = pbmc,
  assay = "peaks",
  regions = blacklist_hg38_unified
)
```

### 🔹 Step 4: Visualize QC

```r
DensityScatter(pbmc, x = "nCount_peaks", y = "TSS.enrichment", log_x = TRUE)

VlnPlot(pbmc,
  features = c("nCount_peaks", "TSS.enrichment", "pct_reads_in_peaks", 
               "nucleosome_signal", "blacklist_ratio"),
  pt.size = 0.1, ncol = 5
)

FragmentHistogram(pbmc, group.by = ifelse(pbmc$nucleosome_signal > 4, "High NS", "Low NS"))
```

### 🔹 Step 5: Filter Low-Quality Cells

```r
pbmc <- subset(pbmc,
  subset = nCount_peaks >= 1000 &
           TSS.enrichment > 4 &
           pct_reads_in_peaks > 15 &
           nucleosome_signal < 4 &
           blacklist_ratio < 0.01
)
```

---

## 📚 References

- 10x Genomics Cell Ranger ATAC Documentation  
- Stuart et al., Signac: single-cell chromatin analysis toolkit  
- Satija Lab Seurat Documentation  

---

## 🧵 To Do

- [ ] Add multi-sample support  
- [ ] Integrate ArchR-based alternative  
- [ ] Export cleaned fragment files for downstream trajectory analysis  
