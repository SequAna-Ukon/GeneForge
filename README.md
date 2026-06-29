# GeneForge v2.1
**A Nextflow Pipeline for Eukaryotic Gene Prediction and Functional Annotation**

GeneForge is a high-throughput Nextflow pipeline designed for the comprehensive structural and functional annotation of eukaryotic genomes. It orchestrates the parallel execution of BRAKER3 and FunAnnotate, evaluates their performance, and provides a unified functional annotation suite.

## Features
* **Dual Annotation Engine**: Parallel execution of BRAKER3 (evidence-based) and FunAnnotate (RNA-seq/Protein-guided).
* **Unified Selection & Integration**: The alternative evaluation layer is fully merged downstream. The highest-quality gene set is dynamically isolated and integrated directly during the consensus build.
* **Multi-Omics Integration**: Supports Short-read RNA-seq (forward, reverse, or unstranded), Long-read RNA-seq (Nanopore/PacBio with length thresholds to optimize `SeqClean` container performance), and protein homology evidence.
* **Comprehensive Functional Suite**: Integrates InterProScan, EggNOG-mapper, Phobius, and SignalP.
* **Containerized & Reproducible**: Built using Nextflow DSL2 with full Singularity/Apptainer support for HPC environments via public registries.

## Prerequisites
* **Nextflow**: Version ≥ 22.04
* **Singularity/Apptainer**: Required for containerized execution.
* **Hardware**: Linux-based system, 16GB+ RAM (32GB+ recommended), 8+ CPU cores.
* **Required Proprietary Files**:
    * **GeneMark**: License key (`gm_key_64.gz`) and Tarball (`gmes_linux_64_4.tar.gz`).
    * **Functional Tools**: Phobius (`phobius101_linux.tgz`) and SignalP (`signalp-6.0h.fast.tar.gz`).

## Installation

```bash
# Clone the repository
git clone [https://github.com/yourusername/GeneForge.git](https://github.com/yourusername/GeneForge.git)
cd GeneForge

# Verify Nextflow installation
nextflow -v
Run ModesGeneForge uses the --mode flag to define the structural annotation strategy. Every mode concludes with an AGAT-based cleanup to resolve overlaps and ensure GFF3 compliance.ModePrediction StrategyFinal Annotation LogicbrakerBRAKER3 onlyUses BRAKER3 as the reference; merges tRNAs from tRNAscan-SE.funannotateFunAnnotate onlyUses FunAnnotate as the reference; merges tRNAs from tRNAscan-SE.bothDual EngineUses the BUSCO winner as the Backbone. Unique, non-overlapping models from the runner-up are added to "complement" the backbone, followed by tRNA integration.Functional Annotation (Optional)Triggered by the --func_annotation flag:Scope: Includes Phobius, SignalP, EggNOG-mapper, and InterProScan.Target: Runs exclusively on the final merged consensus (the output of the modes above).Workflow OverviewThe pipeline follows a modular architecture using Nextflow DSL2 backed by dedicated, public container layers:tRNA Scanning: TRNASCAN_SE identifies eukaryotic tRNAs using abdoallahsharaf/geneforge-trnascan:2.0.12.RNA-seq Processing: RNASEQ_PROCESSING aligns FASTQ files and generates BAM/GTF evidence using abdoallahsharaf/geneforge-rnaseq:2.0.Gene Prediction:BRAKER: Evidence-based gene prediction using RNA-seq and protein homology via abdoallahsharaf/braker3:v2.FUNANNOTATE: Parallel prediction integrating RNA-seq, protein, and tRNA data via abdoallahsharaf/geneforge-funannotate:2.0.Unified Merge and Evaluation (MERGE_ANNOTATIONS):Evaluates metrics on the fly from both prediction pipelines to automatically isolate the highest-quality structural Backbone.The Backbone is complemented with missing models from the alternative tool via abdoallahsharaf/braker3:v2.tRNA annotations are merged into the consensus.AGAT is used to resolve overlaps and finalize GFF3 coordinates.Functional Annotation: FUNCTIONAL_ANNOTATION (Optional) adds functional descriptors to the finalized consensus using abdoallahsharaf/geneforge-funannotate-func:1.8.17.UsageGeneForge uses CSV files to manage complex metadata and file paths.Bashnextflow run GeneForge/main.nf \
  --mandatory_csv mandatory.csv \
  --optional_csv optional.csv \
  --mode both \
  --func_annotation
Input Configuration1. mandatory.csvUsed for core genomic data and species information.Format: name,species,organism,busco_db,busco_db_fun,genome_masked,genome_unmasked,protein_evidence,genemark_dirExample:Code snippetname,species,organism,busco_db,busco_db_fun,genome_masked,genome_unmasked,protein_evidence,genemark_dir
Cther,Cladocopium thermophilum,other,alveolata_odb10,protists,/path/to/Cther.fasta.masked,/path/to/Cther.fasta,/path/to/Alveolata.fa,/path/to/genemark
2. optional.csvUsed for RNA-seq data, databases, and third-party tool directories.Format:rnaseq_dir,funanno_DB,eggnog_DB,stranded,nanopore_mrna,pacbio_isoseq,gc_probability,func_tool_dirExample:Code snippetrnaseq_dir,funanno_DB,eggnog_DB,stranded,nanopore_mrna,pacbio_isoseq,gc_probability,func_tool_dir
/path/to/RNA_Cther,/path/to/funannotate_DB,/path/to/eggnog_DB,reverse,/path/to/ONT.fastq.gz,/path/to/pacbio.fastq.gz,0.6377,/path/to/tools
Note: For functional annotation, func_tool_dir must contain:phobius101_linux.tgz (Phobius tarball).signalp-6.0h.fast.tar.gz (SignalP tarball).Container RegistriesThe framework dynamically coordinates task execution inside the following public images:tRNAscan-SE: docker://abdoallahsharaf/geneforge-trnascan:2.0.12RNA-Seq Suite: docker://abdoallahsharaf/geneforge-rnaseq:2.0BRAKER Framework: docker://abdoallahsharaf/braker3:v2Funannotate Engine: docker://abdoallahsharaf/geneforge-funannotate:2.0Functional Suite: docker://abdoallahsharaf/geneforge-funannotate-func:1.8.17OutputsResults are organized in the results/ directory:geneforge/: The final "Best" annotation set (GFF3, Proteins, BUSCO summary).braker/ & funannotate/: Raw prediction outputs from each individual tool.tRNA_scan/: High-confidence tRNA annotations.functional_annotation/: Integrated results from Phobius, EggNOG, and InterProScan.License & AttributionPrimary LicenseThis pipeline is licensed under the MIT License.Third-Party LicensesGeneForge automates the use of third-party tools. Users are responsible for complying with their respective licenses:GeneMark: Proprietary (Academic use only; commercial license required).BRAKER3/tRNAscan-SE/BUSCO: GPL-3.0.STAR/Samtools: MIT/BSD.Individual license files for these dependencies can be found in /third_party_licenses/.CitationsIf you use GeneForge in your research, please cite:Sharaf, A., & Voolstra, C. R. (2026). GeneForge v2.0: A Nextflow Pipeline for Gene Prediction and Functional Annotation. Zenodo. https://doi.org/10.5281/zenodo.18592773AcknowledgmentsSupported by the Sequencing Analysis (SequAna) Core Facility at the University of Konstanz.Contact: [abdoallah.sharaf@uni-konstanz.de]
