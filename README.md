# GeneForge v2.1
**A Nextflow Pipeline for Eukaryotic Gene Prediction and Functional Annotation**

GeneForge is a high-throughput Nextflow pipeline for comprehensive structural and functional annotation of eukaryotic genomes. It orchestrates parallel execution of BRAKER3 and FunAnnotate, evaluates their performance via BUSCO, and provides a unified functional annotation suite.

---

## What's New in v2.1

### Containerized Architecture
All modules now run inside purpose-built Docker/Singularity containers, eliminating conda environment conflicts across tools. Each container ships with its own isolated software stack:

| Container | Module(s) |
| :--- | :--- |
| `abdoallahsharaf/geneforge-trnascan:2.1` | tRNAscan-SE |
| `abdoallahsharaf/geneforge-rnaseq:2.1` | STAR, StringTie, Trimmomatic |
| `abdoallahsharaf/geneforge-braker3:2.1` | BRAKER3, AUGUSTUS, GeneMark-ETP, AGAT, BUSCO |
| `abdoallahsharaf/geneforge-funannotate:2.1` | FunAnnotate, PASA, Trinity, EvidenceModeler |
| `abdoallahsharaf/geneforge-funannotate-func:2.1` | EggNOG-mapper, InterProScan, Phobius, SignalP6 |

### Long-Read RNA-seq Support
`--nanopore_mrna` and `--pacbio_isoseq` inputs are now fully integrated into the FunAnnotate training step. Long reads are automatically preprocessed: FASTQ→FASTA conversion, U→T substitution, and length filtering (≥200 bp) to prevent seqclean `IndexError` downstream.

### Protein Evidence Pre-filtering
A diamond blastp + seqtk subseq step is now inserted between `funannotate train` and `funannotate predict`. PASA TransDecoder peptides are used as query against the full protein database, and only matching proteins are passed to prediction. This avoids the multi-week exonerate runtimes caused by large databases (e.g., full Metazoa UniProt).

### Strand-Aware BAM Splitting
RNASEQ_PROCESSING now produces strand-specific BAM files (`_plus_strand.bam`, `_minus_strand.bam`) for stranded libraries. These are passed directly to BRAKER via `--bam=plus,minus` with `--stranded=+,-`, improving intron hint accuracy for stranded protocols.

### Robust Dummy-File Handling
Dummy placeholder files are now written to `workflow.workDir` (not `projectDir`) and cleaned up via `workflow.onComplete`. Nanopore and PacBio dummies (indices 11 and 12) were added to prevent Nextflow cache invalidation on re-runs without long reads.

### InterProScan Self-Installation
The functional annotation module now auto-downloads and unpacks the full InterProScan 5.67-99.0 64-bit distribution on first run if not already present, with flock-based protection against parallel downloads.

### Improved Error Resilience
- BRAKER and FUNANNOTATE processes use `errorStrategy = 'ignore'` so a single-tool failure does not abort the pipeline in `both` mode.
- GeneMark-ET failures in FunAnnotate automatically retry with GeneMark-ES.
- All processes initialize required output files before execution to prevent Nextflow tracking crashes on early failures.

---

## Features

**Dual Annotation Engine**: Parallel execution of BRAKER3 (evidence-based) and FunAnnotate (RNA-seq/Protein-guided).  
**Intelligent Selection**: Automatically compares predictions via BUSCO scores to select the highest-quality gene set.  
**Multi-Omics Integration**: Short-read RNA-seq (forward, reverse, or unstranded), Long-read RNA-seq (Nanopore/PacBio), and protein homology evidence.  
**Protein Pre-filtering**: Diamond-based filtering of protein databases before prediction to prevent exonerate bottlenecks.  
**Comprehensive Functional Suite**: InterProScan, EggNOG-mapper, Phobius, and SignalP6.  
**Containerized & Reproducible**: All modules run in versioned Singularity/Apptainer containers.

---

## Prerequisites

- **Nextflow**: Version ≥ 22.04
- **Singularity/Apptainer**: Required for containerized execution
- **Hardware**: Linux-based system, 32GB+ RAM recommended, 16+ CPU cores
- **Required Proprietary Files**:
  - **GeneMark**: License key (`gm_key_64.gz`) and tarball (`gmes_linux_64_4.tar.gz`)
  - **Functional Tools**: Phobius (`phobius101_linux.tgz`) and SignalP (`signalp-6.0h.fast.tar.gz`)

---

## Installation

```bash
git clone https://github.com/yourusername/GeneForge.git
cd GeneForge
nextflow -v
```

---

## Run Modes

GeneForge uses the `--mode` flag to define the structural annotation strategy. Every mode concludes with AGAT-based cleanup to resolve overlaps and ensure GFF3 compliance.

| Mode | Prediction Strategy | Final Annotation Logic |
| :--- | :--- | :--- |
| **`braker`** | BRAKER3 only | Uses BRAKER3 as the reference; merges tRNAs from tRNAscan-SE |
| **`funannotate`** | FunAnnotate only | Uses FunAnnotate as the reference; merges tRNAs from tRNAscan-SE |
| **`both`** | Dual Engine | Uses the BUSCO winner as the backbone; unique non-overlapping models from the runner-up complement it, followed by tRNA integration |

### Functional Annotation (Optional)
Triggered by `--func_annotation`. Runs Phobius, SignalP6, EggNOG-mapper, and InterProScan on the final merged consensus. Databases are downloaded and cached automatically on first run.

---

## Workflow Overview

```
tRNAscan-SE ──────────────────────────────────────────────────────────┐
                                                                       │
RNA-seq (short + long) ──► RNASEQ_PROCESSING                          │
                               │                                       │
                   ┌───────────┴───────────┐                          │
                   ▼                       ▼                          │
              BRAKER3                 FUNANNOTATE                      │
           (evidence-based)      (protein pre-filtered)               │
                   │                       │                          │
                   └───────────┬───────────┘                          │
                               ▼                                      │
                        COMPARE_BUSCO                                  │
                               │                                      │
                               ▼                                      │
                      MERGE_ANNOTATIONS ◄────────────────────────────┘
                               │
                               ▼
                    FUNCTIONAL_ANNOTATION (optional)
```

1. **tRNA Scanning**: `TRNASCAN_SE` identifies eukaryotic tRNAs and generates high-confidence `.tbl` and `.gff` outputs.
2. **RNA-seq Processing**: `RNASEQ_PROCESSING` trims reads, aligns with STAR, assembles with StringTie, and optionally splits BAMs by strand.
3. **Gene Prediction**:
   - `BRAKER`: Evidence-based prediction using RNA-seq BAMs, strand-specific BAMs (when applicable), and protein homology.
   - `FUNANNOTATE`: PASA-based training → diamond protein pre-filtering → EvidenceModeler prediction → PASA update.
4. **Comparison**: `MERGE_ANNOTATIONS` evaluates BUSCO scores to select the backbone annotation.
5. **Merge & Complement**: The backbone is complemented with unique models from the alternative tool; tRNAs are integrated; AGAT resolves overlaps.
6. **Functional Annotation**: `FUNCTIONAL_ANNOTATION` adds InterProScan, EggNOG, Phobius, and SignalP6 functional descriptors.

---

## Usage

```bash
nextflow run GeneForge/main.nf \
  --mandatory_csv mandatory.csv \
  --optional_csv optional.csv \
  --mode both \
  --func_annotation
```

---

## Input Configuration

### `mandatory.csv`

Core genomic data and species information.

**Format**: `name,species,organism,busco_db,busco_db_fun,genome_masked,genome_unmasked,protein_evidence,genemark_dir`

```csv
name,species,organism,busco_db,busco_db_fun,genome_masked,genome_unmasked,protein_evidence,genemark_dir
Aip,Exaiptasia diaphana,other,metazoa_odb10,metazoa,/path/to/genome.fasta.masked,/path/to/genome.fasta,/path/to/proteins.fa,/path/to/genemark
```

| Field | Description |
| :--- | :--- |
| `name` | Sample identifier used as prefix for all outputs |
| `species` | Full species name (spaces allowed; internally converted to underscores) |
| `organism` | Funannotate organism type: `other`, `fungus`, `vertebrate`, etc. |
| `busco_db` | BUSCO lineage for structural evaluation (e.g., `metazoa_odb10`) |
| `busco_db_fun` | BUSCO lineage for funannotate internal use (e.g., `metazoa`) |
| `genome_masked` | Path to soft-masked genome FASTA |
| `genome_unmasked` | Path to unmasked genome FASTA |
| `protein_evidence` | Protein database FASTA (pre-filtered internally via diamond) |
| `genemark_dir` | Directory containing `gm_key_64.gz` and `gmes_linux_64_4.tar.gz` |

### `optional.csv`

RNA-seq data, databases, and third-party tool paths.

**Format**: `rnaseq_dir,funanno_DB,eggnog_DB,stranded,nanopore_mrna,pacbio_isoseq,gc_probability,func_tool_dir`

```csv
rnaseq_dir,funanno_DB,eggnog_DB,stranded,nanopore_mrna,pacbio_isoseq,gc_probability,func_tool_dir
/path/to/rnaseq,/path/to/funannotate_DB,/path/to/eggnog_DB,reverse,/path/to/ONT.fastq.gz,,0.6377,/path/to/tools
```

| Field | Description |
| :--- | :--- |
| `rnaseq_dir` | Directory of paired-end FASTQ files (`*_R1*` / `*_R2*` naming) |
| `funanno_DB` | Path to existing funannotate database (auto-installed if absent) |
| `eggnog_DB` | Path to EggNOG database (auto-downloaded if absent) |
| `stranded` | Library strandedness: `forward`, `reverse`, or `no` |
| `nanopore_mrna` | Path to ONT direct-RNA FASTQ/FASTA (optional; leave blank if absent) |
| `pacbio_isoseq` | Path to PacBio IsoSeq FASTQ/FASTA (optional; leave blank if absent) |
| `gc_probability` | GC content prior for GeneMark (optional) |
| `func_tool_dir` | Directory containing Phobius and SignalP tarballs (required for `--func_annotation`) |

**Note**: For `--func_annotation`, `func_tool_dir` must contain `phobius101_linux.tgz` and `signalp-6.0h.fast.tar.gz`.

---

## Outputs

Results are organized under `results/`:

| Directory | Contents |
| :--- | :--- |
| `geneforge/` | Final merged annotation: GFF3, proteins FASTA, BUSCO summary |
| `braker/` | Raw BRAKER3 outputs: GFF3, proteins, BUSCO summary, error log |
| `funannotate/` | Raw FunAnnotate outputs: GFF3, proteins, BUSCO summary |
| `tRNA_scan/` | High-confidence tRNA annotations (GFF3, `.tbl`) |
| `RNASeq/` | STAR BAM, StringTie GTF, transcripts FASTA, strand-split BAMs |
| `busco_comparison/` | BRAKER3 vs FunAnnotate comparison report |
| `functional_annotation/` | InterProScan XML, EggNOG annotations, Phobius/SignalP results, funannotate synthesis |

---

## Resource Configuration

Default resource labels (adjustable in `nextflow.config`):

| Label | CPUs | Memory |
| :--- | :--- | :--- |
| `process_low` | 8 | 16 GB |
| `process_medium` | 30 | 50 GB |
| `process_high` | 50 | 100 GB |

Gene prediction processes (`BRAKER`, `FUNANNOTATE`) use `process_high` and run with `errorStrategy = 'ignore'` so a failure in one does not abort the other.

---

## License & Attribution

### Primary License
This pipeline is licensed under the **MIT License**.

### Third-Party Licenses
GeneForge automates the use of third-party tools. Users are responsible for complying with their respective licenses:

- **GeneMark**: Proprietary (academic use only; commercial license required)
- **BRAKER3 / tRNAscan-SE / BUSCO**: GPL-3.0
- **STAR / Samtools / StringTie**: MIT/BSD
- **InterProScan**: Apache 2.0
- **EggNOG-mapper**: LGPL
- **Phobius / SignalP6**: Proprietary (academic use; register and download separately)

Individual license files are in `third_party_licenses/`.

---

## Citation

If you use GeneForge in your research, please cite:

**Sharaf, A., & Voolstra, C. R. (2026)**. GeneForge v2.1: A Nextflow Pipeline for Gene Prediction and Functional Annotation. Zenodo. https://doi.org/10.5281/zenodo.18592773

---

## Acknowledgments

Supported by the Sequencing Analysis (SequAna) Core Facility at the University of Konstanz.

Contact: [abdoallah.sharaf@uni-konstanz.de](mailto:abdoallah.sharaf@uni-konstanz.de)
