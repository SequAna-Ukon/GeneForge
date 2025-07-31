# GeneForge
A Nextflow Pipeline for Gene Prediction and Functional Annotation

**GeneForge** is a flexible Nextflow pipeline for gene prediction and optional functional annotation of eukaryotic genomes. It integrates BRAKER3 and FunAnnotate for gene prediction, compares their quality using BUSCO, and supports functional annotation with tools such as Phobius, InterProScan, and eggNOG-mapper. The pipeline offers three run modes: BRAKER3 only, FunAnnotate only, or both, with customizable parameters for species, input paths, and annotation preferences.

## Features

- **Dual Annotation**: Combines BRAKER3 (evidence-based gene prediction) and Funannotate (RNA-seq and protein-guided annotation).
- **Flexible Modes**: Run BRAKER3 (```braker```), Funannotate (```funannotate```), or both (```both```).
- **Stranded RNA-seq Support**: Handles ```forward```, ```reverse```, or unstranded data with automatic BAM file processing.
- **Functional Annotation**: Integrates Funannotate's functional annotation with databases like EggNOG and BUSCO.
- **Error Handling**: Robust validation of inputs and detailed logging (```*_error.log```).
- **Modular Design**: Uses Nextflow DSL2 for scalability and reproducibility.

## Prerequisites

- **Nextflow**: Version ≥ 22.04 
- **Singularity**: Required for containerized execution. Install via your package manager or see [Singularity documentation](https://sylabs.io/docs/).
- **System**: Linux/Unix-based system with sufficient RAM (≥16GB recommended) and CPU cores (≥8 recommended).
- **Disk Space**: Depends on genome size and RNA-Seq data; 100GB+ recommended.
- **Required Input Files**:
    - Masked and unmasked genome FASTA.
    - Protein evidence FASTA.
    - **GeneMark License Key** (```gm_key_64.gz```) and **GeneMark Tarball** (```gmes_linux_64_4.tar.gz```).
    - For functional annotation:
        - **Phobius Tarball** (```phobius101_linux.tgz```).
        - **SignalP Tarball** (```signalp-6.0h.fast.tar.gz```).
    - Optional: RNA-seq FASTQ, Funannotate/EggNOG databases.

## Installation

1. **Clone the Repository**:
````bash
git clone https://github.com/yourusername/GeneForge.git
````
2. Ensure Singularity is installed and configured.
3. Provide paths to required input files (see Input Files).
## Usage

Run the pipeline with:

```bash
nextflow run main_full.nf \
  --mandatory_csv mandatory.csv \
  --optional_csv optional.csv \
  --mode both \
  --func_annotation
```
## Parameters

- ```--mandatory_csv```: CSV with required inputs (see Input Files).
- ```--optional_csv```: CSV with optional inputs.
- ```--mode```: ```both```, ```braker```, or ```funannotate```.
- ```--func_annotation```: Enable or disable functional annotation.


## Input Files

```mandatory.csv```

Format: 
```name,species,organism,busco_db,busco_db_fun,genome_masked,genome_unmasked,protein_evidence,genemark_dir```

Example:
````csv
name,species,organism,busco_db,busco_db_fun,genome_masked,genome_unmasked,protein_evidence,genemark_dir
Cther,Cladocopium thermophilum,other,alveolata_odb10,protists,/path/to/Cther.fasta.masked,/path/to/Cther.fasta,/path/to/Alveolata.fa,/path/to/genemark
````
**Note**: ```genemark_dir``` must contain:
- ```gm_key_64.gz``` (GeneMark license key).
- ```gmes_linux_64_4.tar.gz``` (GeneMark tarball).

```optional.csv```

Format: 
```rnaseq_dir,funanno_DB,eggnog_DB,stranded,nanopore_mrna,pacbio_isoseq,gc_probability,func_tool_dir```

Example:
````csv
rnaseq_dir,funanno_DB,eggnog_DB,stranded,nanopore_mrna,pacbio_isoseq,gc_probability,func_tool_dir
/path/to/RNA_Cther,/path/to/funannotate_DB,/path/to/eggnog_DB,reverse,/path/to/ONT.fastq.gz,/path/to/pacbio.fastq.gz,0.6377,/path/to/tools
````
**Note**: For functional annotation, ```func_tool_dir``` must contain:
- ```phobius101_linux.tgz``` (**Phobius tarball**).
- ```signalp-6.0h.fast.tar.gz``` (**SignalP tarball**).
  
## Outputs
- **RNASeq** (```results/RNASeq/```):[optional: if RNASeq data were used]
    - ```${name}_RNASeqAll.Stringtie.gtf```: Prepared Stringtie evidence.
    - ```${name}_RNASeqAll.STAR.bam```: Prepared Smapping file evidence.
    - ```${name}_RNASeqAll.transcripts.fasta```: Prepared transcripts evidence.
    - ```${name}_plus_strand.bam```: Prepared Smapping file (forward strand) evidence. [optional: if the ```stranded``` option is used]
    - ```${name}_minus_strand.bam```: Prepared Smapping file (reverse strand) evidence. [optional: if the ```stranded``` option is used]

- **tRNAScan_SE** (```results/tRNA_scan/```):
    - ```${name}_trna_annotation.gff```: Eukaryotic tRNAs annotation.
    - ```${name}_highconf.tbl```: Eukaryotic tRNAs.
    
- **Funannotate** (```results/funannotate/```):
    - ```${name}_funannotate.gff3```: Funannotate Gene predictions.
    - ```${name}.funannotate.prot.fasta```: Funannotate Protein sequences.
    - ```${name}_busco_funannotate.txt```: Funannotate BUSCO summary.
      
- **BRAKER3** (```results/braker/```):
    - ```${name}_braker.gff3```: Braker Gene predictions.
    - ```${name}.braker.prot.fasta```: Braker Protein sequences.
    - ```${name}_busco_braker.txt```: Braker BUSCO summary.
      
- **BUSCO Comparison** (```results/busco_comparison/```): [optional: if the ```both``` option is used]
    - ```busco_comparison.txt```: A report that summeries the BUSCO scores and which gene prediction was used for the functional annotation if  the ```--func_annotation``` option is used.
      
- **Functional Annotation** (```results/functional_annotation/```): [optional: if the ```--func_annotation``` option is used]
    - ```${name}_functional_annotation/annotate_results/*```: Directory contains all funannotate annotattion outputs.

## Workflow Overview

1. **tRNA Scanning**: ```TRNASCAN_SE``` identifies tRNAs.
2. **RNA-seq Processing**: ```RNASEQ_PROCESSING``` aligns FASTQ files and generates BAM/GTF.
3. **Gene Prediction**:
    - ```BRAKER_RUN```: Uses RNA-seq and protein evidence for gene prediction.
    - ```FUNANNOTATE```: Integrates RNA-seq, protein, and tRNA data for annotation.
4. **Post-processing**: ```BRAKER_POST``` refines BRAKER3 outputs.
5. **Comparison**: ```COMPARE_BUSCO``` evaluates Funannotate vs. BRAKER3 using BUSCO scores.
6. **Functional Annotation**: ```FUNCTIONAL_ANNOTATION``` adds functional annotations using funannotate.

![GeneForge Workflow Overview](workflow_sct.png)

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Cite Us

If you use this software, please cite it as below:

``` commandline
Sharaf, A., & Voolstra, C. R. (2025).
SequAna-Ukon/GeneForge: GeneForge v1.0 (v1.0). Zenodo.
https://doi.org/10.5281/zenodo.16631467
```

## Acknowledgments
This work is supported through the Sequencing analysis (SequAna) core facility at the University of Konstanz [https://www.biologie.uni-konstanz.de/sequana/]

**Tools**: BRAKER3, funannotate, BUSCO, Trimmomatic, STAR, StringTie, tRNAscan-SE, Phobius, InterProScan, eggNOG-mapper.
**Nextflow**: For workflow orchestration.


For more information or help, please contact [abdoallah.sharaf@uni-konstanz.de].
