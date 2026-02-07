# GeneForge v2.0
**A Nextflow Pipeline for Eukaryotic Gene Prediction and Functional Annotation**
GeneForge is a high-throughput Nextflow pipeline designed for the comprehensive structural and functional annotation of eukaryotic genomes. It orchestrates the parallel execution of BRAKER3 and FunAnnotate, evaluates their performance using BUSCO, and provides a unified functional annotation suite.

## Features
**Dual Annotation Engine**: Parallel execution of BRAKER3 (evidence-based) and FunAnnotate (RNA-seq/Protein-guided).
**Intelligent Selection**: Automatically compares predictions via BUSCO scores to select the highest-quality gene set for downstream annotation.
**Multi-Omics Integration**: Supports Short-read RNA-seq (forward, reverse, or unstranded), Long-read RNA-seq (Nanopore/PacBio), and protein homology evidence.
**Comprehensive Functional Suite**: Integrates InterProScan, EggNOG-mapper, Phobius, and SignalP.
**Containerized & Reproducible**: Built using Nextflow DSL2 with full Singularity/Apptainer support for HPC environments.

## Prerequisites
- **Nextflow**: Version ≥ 22.04
- **Singularity/Apptainer**: Required for containerized execution.
- **Hardware**: Linux-based system, 16GB+ RAM (32GB+ recommended), 8+ CPU cores.
- **Required Proprietary Files**:
    - **GeneMark**: License key (```gm_key_64.gz```) and Tarball (```gmes_linux_64_4.tar.gz```).
    - **Functional Tools**: Phobius (```phobius101_linux.tgz```) and SignalP (```signalp-6.0h.fast.tar.gz```).


## Installation

````bash
# Clone the repository
git clone https://github.com/yourusername/GeneForge.git
cd GeneForge

# Verify Nextflow installation
nextflow -v
````

## Run Modes

GeneForge uses the `--mode` flag to define the structural annotation strategy. Every mode concludes with an **AGAT-based cleanup** to resolve overlaps and ensure GFF3 compliance.

| Mode | Prediction Strategy | Final Annotation Logic |
| :--- | :--- | :--- |
| **`braker`** | BRAKER3 only | Uses BRAKER3 as the reference; merges tRNAs from tRNAscan-SE. |
| **`funannotate`** | FunAnnotate only | Uses FunAnnotate as the reference; merges tRNAs from tRNAscan-SE. |
| **`both`** | Dual Engine | Uses the BUSCO winner as the **Backbone**. Unique, non-overlapping models from the runner-up are added to "complement" the backbone, followed by tRNA integration. |

### Functional Annotation (Optional)
Triggered by the `--func_annotation` flag:
- **Scope**: Includes Phobius, SignalP, EggNOG-mapper, and InterProScan.
- **Target**: Runs exclusively on the **final merged consensus** (the output of the modes above).

---

## Workflow Overview

The pipeline follows a modular architecture using Nextflow DSL2:

1. **tRNA Scanning**: `TRNASCAN_SE` identifies eukaryotic tRNAs.
2. **RNA-seq Processing**: `RNASEQ_PROCESSING` aligns FASTQ files and generates BAM/GTF evidence.
3. **Gene Prediction**:
    * `BRAKER_RUN`: Evidence-based gene prediction using RNA-seq and protein homology.
    * `FUNANNOTATE`: Parallel prediction integrating RNA-seq, protein, and tRNA data.
4. **Post-processing**: `BRAKER_POST` standardizes formats and prepares outputs for comparison.
5. **Comparison**: `COMPARE_BUSCO` evaluates results to select the highest-quality "Backbone."
6. **Merge & Complement (`MERGE_ANNOTATIONS`)**: 
    * The Backbone is complemented with missing models from the alternative tool.
    * tRNA annotations are merged into the consensus.
    * **AGAT** is used to resolve overlaps and finalize GFF3 coordinates.
7. **Functional Annotation**: `FUNCTIONAL_ANNOTATION` (Optional) adds functional descriptors to the finalized consensus.

![GeneForge Workflow Overview](workflow_sct.jpg)

## Usage

GeneForge uses CSV files to manage complex metadata and file paths.
````bash
nextflow run GeneForge/main.nf \
  --mandatory_csv mandatory.csv \
  --optional_csv optional.csv \
  --mode both \
  --func_annotation
````

## Input Configuration
**1.** ```mandatory.csv```
Used for core genomic data and species information. 

````
name,species,organism,busco_db,genome_masked,genome_unmasked,genemark_dir
Cther,Cladocopium thermophilum,other,protists,Cther.masked.fasta,Cther.fasta,/path/to/genemark
````

**2.** ```optional.csv```
Used for RNA-seq data, databases, and third-party tool directories. 

````
rnaseq_dir,funanno_DB,eggnog_DB,stranded,nanopore_mrna,pacbio_isoseq,gc_probability,func_tool_dir
/path/to/RNA_Cther,/path/to/funannotate_DB,/path/to/eggnog_DB,reverse,,,0.6377,/path/to/tools
````
## Outputs
Results are organized in the ```results/``` directory:

- ```geneforge/```: The final "Best" annotation set (GFF3, Proteins, BUSCO summary).
- ```braker/``` **&** ```funannotate/```: Raw prediction outputs from each individual tool.
- ```tRNA_scan/```: High-confidence tRNA annotations.
- ```busco_comparison/```: A summary report comparing BRAKER3 vs. FunAnnotate scores.
- ```functional_annotation/```: Integrated results from Phobius, EggNOG, and InterProScan.

## License & Attribution
### Primary License
This pipeline is licensed under the **MIT License**.
### Third-Party Licenses
GeneForge automates the use of third-party tools. Users are responsible for complying with their respective licenses:

- **GeneMark**: Proprietary (Academic use only; commercial license required).
- **BRAKER3/tRNAscan-SE/BUSCO**: GPL-3.0.
- **STAR/Samtools**: MIT/BSD.
- Individual license files for these dependencies can be found in ```/third_party_licenses/```.

## Citations
If you use GeneForge in your research, please cite:

**Sharaf, A., & Voolstra, C. R. (2025)**. GeneForge: A Nextflow Pipeline for Gene Prediction and Functional Annotation. Zenodo. https://doi.org/10.5281/zenodo.16631467

## Acknowledgments
Supported by the Sequencing Analysis (SequAna) Core Facility at the University of Konstanz.

Contact: [abdoallah.sharaf@uni-konstanz.de]

Would you like me to help you generate a CHANGELOG.md or a CONTRIBUTING.md file to further professionalize the repository?









