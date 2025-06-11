nextflow.enable.dsl=2

// Pipeline Description:
// GeneForge performs gene prediction and optional functional annotation using BRAKER3 and funannotate.
// - Prediction: Runs BRAKER3 and funannotate, compares BUSCO scores.
// - Annotation: Optionally annotates the selected prediction using Phobius, InterProScan, eggNOG-mapper, and funannotate.
// - RNA-Seq: Optional; if not provided, predictions rely on genome and protein evidence.
// - Run Modes: prediction_only, annotation_only, both.

// Default parameters
params {
    rnaseq_dir = null  // Optional RNA-Seq directory
    genome_fasta = "${params.docker_mount_base}/genome.fasta"
    masked_genome = "${params.docker_mount_base}/genome_masked.fasta"
    protein_evidence = "${params.docker_mount_base}/proteins.fa"
    adapters = "${params.docker_mount_base}/adapters.fa"
    outdir = "${params.docker_mount_base}/results"
    threads = 8
    memory = "16G"
    busco_db = "protists"
    funannotate_db = "/opt/funannotate_db"
    species = "Unknown_species"
    min_protlen = 100
    star_index_nbases = 10
    bam_sort_memory = "10G"
    gc_probability = false
    busco_lineage = "protists_odb10"
    enable_functional_annotation = false
    run_mode = "both"  // Options: prediction_only, annotation_only, both
    annotation_prediction = "best_busco"  // Options: braker, funannotate, best_busco
    docker_mount_base = "/workspace"  // Docker mount point inside container
}

// Validate inputs
if (!file(params.genome_fasta).exists()) { error "Genome FASTA file not found: ${params.genome_fasta}" }
if (!file(params.masked_genome).exists()) { error "Masked genome FASTA file not found: ${params.masked_genome}" }
if (!file(params.protein_evidence).exists()) { error "Protein evidence file not found: ${params.protein_evidence}" }
if (!file(params.funannotate_db).isDirectory()) { error "Funannotate database directory not found: ${params.funannotate_db}" }
if (params.rnaseq_dir != null && !file(params.rnaseq_dir).isDirectory()) { error "RNA-Seq directory not found: ${params.rnaseq_dir}" }
if (params.rnaseq_dir != null && !file(params.adapters).exists()) { error "Adapter FASTA file not found: ${params.adapters}" }

// Validate species name (no special characters except underscore)
if (!params.species.matches("[a-zA-Z0-9_ ]+")) { error "Invalid species name: ${params.species}. Use alphanumeric, underscore, or space only." }

// Validate run_mode
if (!["prediction_only", "annotation_only", "both"].contains(params.run_mode)) { error "Invalid run_mode: ${params.run_mode}. Valid options: prediction_only, annotation_only, both" }

// Validate annotation_prediction
if (!["braker", "funannotate", "best_busco"].contains(params.annotation_prediction)) { error "Invalid annotation_prediction: ${params.annotation_prediction}. Valid options: braker, funannotate, best_busco" }

// Validate annotation_only requirements
if (params.run_mode == "annotation_only" && !params.enable_functional_annotation) { error "run_mode annotation_only requires enable_functional_annotation = true" }
if (params.run_mode == "annotation_only" && params.annotation_prediction == "best_busco") { error "annotation_only mode requires annotation_prediction to be 'braker' or 'funannotate', not 'best_busco'" }

// Check if RNA-Seq data is provided
def has_rnaseq = params.rnaseq_dir != null && file(params.rnaseq_dir).isDirectory()

// Collect paired-end RNA-Seq reads if provided
read_pairs = has_rnaseq ? Channel.fromFilePairs("${params.rnaseq_dir}/*_{R1,R2}_*.{fastq,fq}.gz", checkIfExists: true) : null

// Process 1: Clean RNA-Seq reads with Trimmomatic
process CleanReads {
    tag "${sample_id}"
    publishDir "${params.outdir}/cleaned_reads", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_1P.fq.gz"), path("${sample_id}_2P.fq.gz"), emit: paired_reads
    path "${sample_id}_1U.fq.gz", emit: unpaired_r1
    path "${sample_id}_2U.fq.gz", emit: unpaired_r2

    script:
    """
    java -jar /opt/trimmomatic/trimmomatic-0.39.jar PE -threads ${params.threads} -phred33 \
        ${reads[0]} ${reads[1]} \
        ${sample_id}_1P.fq.gz ${sample_id}_1U.fq.gz \
        ${sample_id}_2P.fq.gz ${sample_id}_2U.fq.gz \
        ILLUMINACLIP:${params.adapters}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Process 2: Index genome with STAR
process IndexGenome {
    publishDir "${params.outdir}/index", mode: 'copy'

    input:
    path genome

    output:
    path "star_index/*", emit: index

    script:
    """
    mkdir star_index
    STAR --runThreadN ${params.threads} \
         --runMode genomeGenerate \
         --genomeDir star_index \
         --genomeFastaFiles ${genome} \
         --genomeSAindexNbases ${params.star_index_nbases}
    """
}

// Process 3: Map RNA-Seq reads to genome with STAR
process MapReads {
    tag "${sample_id}"
    publishDir "${params.outdir}/mapped", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path index

    output:
    path "${sample_id}_Aligned.sortedByCoord.out.bam", emit: bam

    script:
    """
    STAR --runThreadN ${params.threads} \
         --genomeDir ${index} \
         --readFilesIn ${read1} ${read2} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMstrandField intronMotif \
         --outFilterIntronMotifs RemoveNoncanonical \
         --outFileNamePrefix ${sample_id}_ \
         --limitBAMsortRAM ${params.bam_sort_memory}
    """
}

// Process 4: Merge BAM files
process MergeBAMs {
    publishDir "${params.outdir}/merged", mode: 'copy'

    input:
    path bams

    output:
    path "merged.bam", emit: merged_bam

    script:
    """
    samtools merge -@ ${params.threads} merged.bam ${bams}
    """
}

// Process 5: Generate expression evidence with StringTie
process StringTie {
    publishDir "${params.outdir}/stringtie", mode: 'copy'

    input:
    path bam

    output:
    path "stringtie.gtf", emit: gtf
    path "stringtie_stats.txt", emit: stats

    script:
    """
    stringtie -p ${params.threads} -o stringtie.gtf ${bam}
    grep -v "#" stringtie.gtf | cut -f3 | sort | uniq -c > stringtie_stats.txt
    """
}

// Process 6: Convert GTF to transcript FASTA
process GTFtoFasta {
    publishDir "${params.outdir}/transcripts", mode: 'copy'

    input:
    path gtf
    path genome

    output:
    path "transcripts.fasta", emit: transcript_fasta

    script:
    """
    gtf_genome_to_cdna_fasta.pl ${gtf} ${genome} > transcripts.fasta
    """
}

// Process 7: Predict tRNAs
process tRNAScan {
    publishDir "${params.outdir}/trnascan", mode: 'copy'

    input:
    path masked_genome

    output:
    path "trnascan/eukconf/*", emit: trna_results
    path "trnascan/trnascan-se.out", emit: trna_out

    script:
    """
    tRNAscan-SE -E -I -H --detail --thread ${params.threads} \
        -o trnascan-se.out -f trnascan-se.tbl -m trnascan-se.log \
        ${masked_genome}
    EukHighConfidenceFilter -i trnascan-se.out -s trnascan-se.tbl -o eukconf -p filt
    """
}

// Process 8: Train funannotate
process FunannotateTrain {
    publishDir "${params.outdir}/funannotate_train", mode: 'copy'

    input:
    path genome
    path paired_reads
    val species

    output:
    path "funannotate_train/*", emit: training_data

    script:
    def reads = paired_reads.collect { it.toString() }.join(' ')
    """
    funannotate train -i ${genome} \
        -o funannotate_train \
        --species "${species}" \
        --left ${reads.findAll { it.contains('_1P.') }} \
        --right ${reads.findAll { it.contains('_2P.') }} \
        --cpus ${params.threads} \
        --memory ${params.memory}
    """
}

// Process 9: Gene prediction with funannotate
process FunannotatePredict {
    publishDir "${params.outdir}/funannotate_predict", mode: 'copy'

    input:
    path masked_genome
    path bam, stageAs: 'bam/*'
    path gtf
    path transcript_fasta
    path trna_results
    val species

    output:
    path "funannotate_predict/*", emit: prediction_results
    path "funannotate_predict/update_results/*.gff3", emit: gff3
    path "funannotate_predict/update_results/*.proteins.fa", emit: proteins

    script:
    def bam_opt = bam ? "--rna_bam ${bam}" : ""
    def gtf_opt = gtf ? "--stringtie ${gtf}" : ""
    def transcript_opt = transcript_fasta ? "--transcript_evidence ${transcript_fasta}" : ""
    """
    funannotate predict -i ${masked_genome} \
        -s "${species}" \
        -o funannotate_predict \
        --name ${species.replace(' ', '_')} \
        ${bam_opt} \
        ${gtf_opt} \
        --protein_evidence ${params.protein_evidence} \
        ${transcript_opt} \
        --trnascan ${trna_results}/filt \
        --organism other \
        --database ${params.funannotate_db} \
        --busco_db ${params.busco_db} \
        --min_protlen ${params.min_protlen} \
        --cpus ${params.threads}
    """
}

// Process 10: Add UTRs with funannotate
process FunannotateUpdate {
    publishDir "${params.outdir}/funannotate_update", mode: 'copy'

    input:
    path prediction_results

    output:
    path "funannotate_update/*", emit: updated_results
    path "funannotate_update/update_results/*.gff3", emit: gff3
    path "funannotate_update/update_results/*.proteins.fa", emit: proteins

    script:
    """
    funannotate update -i ${prediction_results} --cpus ${params.threads}
    """
}

// Process 11: Gene prediction with BRAKER3
process BRAKER3 {
    publishDir "${params.outdir}/braker3", mode: 'copy'
    container 'teambraker/braker3:latest'

    input:
    path masked_genome
    path bam, stageAs: 'bam/*'
    path protein_evidence
    val species

    output:
    path "braker3/*", emit: braker_results
    path "braker3/braker.gtf", emit: braker_gtf, optional: true

    script:
    def gc_option = params.gc_probability != false ? "--gc_probability=${params.gc_probability}" : ""
    def bam_opt = bam ? "--bam=${bam}" : ""
    """
    braker.pl --species=${species.replace(' ', '_')} \
        --genome=${masked_genome} \
        ${bam_opt} \
        ${gc_option} \
        --threads ${params.threads} \
        --prot_seq=${protein_evidence} \
        --busco_lineage=${params.busco_lineage} \
        --workingdir=braker3
    """
}

// Process 12: Add UTRs to BRAKER3 GTF using stringtie2utr.py
process AddUTRs {
    publishDir "${params.outdir}/braker3_utr", mode: 'copy'

    input:
    path braker_gtf
    path stringtie_gtf

    output:
    path "braker_with_utrs.gtf", emit: braker_utr_gtf

    script:
    """
    wget https://raw.githubusercontent.com/Gaius-Augustus/BRAKER/utr_from_stringtie/scripts/stringtie2utr.py
    python3 stringtie2utr.py -g ${braker_gtf} -s ${stringtie_gtf} -o braker_with_utrs.gtf
    """
}

// Process 13: Convert tRNAscan-SE output to GFF3
process Convert_tRNA_to_GFF {
    publishDir "${params.outdir}/trnascan", mode: 'copy'

    input:
    path trna_out

    output:
    path "trna_annotation.gff", emit: trna_gff

    script:
    """
    wget https://raw.githubusercontent.com/Gaius-Augustus/BRAKER/master/scripts/convert_tRNAScanSE_to_gff3.pl
    perl convert_tRNAScanSE_to_gff3.pl --input=${trna_out} > trna_annotation.gff
    """
}

// Process 14: Convert BRAKER3 GTF with UTRs to GFF3
process Convert_GTF_to_GFF {
    publishDir "${params.outdir}/braker3_utr", mode: 'copy'

    input:
    path braker_utr_gtf

    output:
    path "braker.gff3", emit: braker_gff3

    script:
    """
    wget https://raw.githubusercontent.com/Gaius-Augustus/BRAKER/master/scripts/gtf2gff.pl
    cat ${braker_utr_gtf} | perl gtf2gff.pl --gff3 -o braker.gff3
    """
}

// Process 15: Merge BRAKER3 GFF3 with tRNA GFF
process MergeGFFs {
    publishDir "${params.outdir}/merged_annotations", mode: 'copy'

    input:
    path braker_gff3
    path trna_gff

    output:
    path "merged.gff", emit: merged_gff

    script:
    """
    agat_sp_merge_annotations.pl --gff ${braker_gff3} --gff ${trna_gff} --out merged.gff
    """
}

// Process 16: Fix overlapping genes in merged GFF
process FixOverlappingGenes {
    publishDir "${params.outdir}/merged_annotations", mode: 'copy'

    input:
    path merged_gff

    output:
    path "fixed_merged.gff", emit: fixed_gff

    script:
    """
    agat_sp_fix_overlaping_genes.pl --gff ${merged_gff} --out fixed_merged.gff
    """
}

// Process 17: Validate GFF3 file
process ValidateGFF {
    publishDir "${params.outdir}/merged_annotations", mode: 'copy'

    input:
    path fixed_gff

    output:
    path "gff_validation.txt", emit: validation_report

    script:
    """
    gt gff3validator ${fixed_gff} > gff_validation.txt 2>&1
    """
}

// Process 18: Export protein sequences for BRAKER3
process ExportProteinsBRAKER {
    publishDir "${params.outdir}/proteins/braker", mode: 'copy'

    input:
    path fixed_gff
    path genome

    output:
    path "braker_proteins.fasta", emit: proteins

    script:
    """
    gffread ${fixed_gff} -g ${genome} -y braker_proteins.fasta
    """
}

// Process 19: Quality control with BUSCO for funannotate
process BUSCO_Funannotate {
    publishDir "${params.outdir}/busco/funannotate", mode: 'copy'

    input:
    path proteins

    output:
    path "busco_funannotate/*", emit: busco_results
    path "busco_funannotate/short_summary.txt", emit: summary

    script:
    """
    busco -i ${proteins} \
        -m proteins \
        -l ${params.busco_lineage} \
        -c ${params.threads} \
        -o busco_funannotate
    """
}

// Process 20: Quality control with BUSCO for BRAKER3
process BUSCO_BRAKER3 {
    publishDir "${params.outdir}/busco/braker3", mode: 'copy'

    input:
    path proteins

    output:
    path "busco_braker3/*", emit: busco_results
    path "busco_braker3/short_summary.txt", emit: summary

    script:
    """
    busco -i ${proteins} \
        -m proteins \
        -l ${params.busco_lineage} \
        -c ${params.threads} \
        -o busco_braker3
    """
}

// Process 21: Compare BUSCO results
process CompareBUSCO {
    publishDir "${params.outdir}/busco", mode: 'copy'

    input:
    path funannotate_summary
    path braker_summary

    output:
    path "busco_comparison.txt", emit: comparison
    tuple val(best_prediction), path(best_gff), path(best_proteins), emit: best_prediction

    script:
    """
    #!/usr/bin/env python3
    import re
    import shutil

    def parse_busco_summary(file_path):
        with open(file_path, 'r') as f:
            content = f.read()
            match = re.search(r'(\d+\.\d+)%\[\d+\]\s*Complete', content)
            if match:
                return float(match.group(1))
            return 0.0

    funannotate_score = parse_busco_summary("${funannotate_summary}")
    braker_score = parse_busco_summary("${braker_summary}")

    with open("busco_comparison.txt", 'w') as f:
        f.write(f"Funannotate BUSCO Completeness: {funannotate_score}%\n")
        f.write(f"BRAKER3 BUSCO Completeness: {braker_score}%\n")
        if funannotate_score > braker_score:
            f.write("Best prediction: Funannotate\n")
            best_pred = "funannotate"
        elif braker_score > funannotate_score:
            f.write("Best prediction: BRAKER3\n")
            best_pred = "braker"
        else:
            f.write("Equal BUSCO scores; defaulting to Funannotate\n")
            best_pred = "funannotate"

    with open("best_prediction.txt", 'w') as f:
        f.write(best_pred)

    if best_pred == "funannotate":
        shutil.copy("${params.outdir}/funannotate_update/update_results/${params.species.replace(' ', '_')}.gff3", "best.gff3")
        shutil.copy("${params.outdir}/funannotate_update/update_results/${params.species.replace(' ', '_')}.proteins.fa", "best.proteins.fa")
    else:
        shutil.copy("${params.outdir}/merged_annotations/fixed_merged.gff", "best.gff3")
        shutil.copy("${params.outdir}/proteins/braker/braker_proteins.fasta", "best.proteins.fa")
    """
}

// Process 22: Run Phobius for transmembrane topology and signal peptide prediction
process Phobius {
    publishDir "${params.outdir}/functional_annotation/phobius", mode: 'copy'

    input:
    path proteins

    output:
    path "phobius.results.txt", emit: phobius_results

    script:
    """
    phobius.pl -short ${proteins} > phobius.results.txt
    """
}

// Process 23: Run InterProScan
process InterProScan {
    publishDir "${params.outdir}/functional_annotation/interproscan", mode: 'copy'

    input:
    path proteins

    output:
    path "Cther_iprscan.xml", emit: iprscan_results

    script:
    """
    funannotate iprscan -i ${proteins} -m docker -c ${params.threads} -o Cther_iprscan.xml
    """
}

// Process 24: Run eggNOG-mapper
process EggNOGMapper {
    publishDir "${params.outdir}/functional_annotation/eggnog", mode: 'copy'

    input:
    path proteins

    output:
    path "Cther_eggnog.emapper.annotations", emit: eggnog_results

    script:
    """
    emapper.py --cpu ${params.threads} -m mmseqs --data_dir ${params.funannotate_db} \
        -i ${proteins} -o Cther_eggnog
    """
}

// Process 25: Functional annotation with funannotate
process FunannotateAnnotate {
    publishDir "${params.outdir}/functional_annotation/anno", mode: 'copy'

    input:
    path gff
    path genome
    path proteins
    path phobius_results
    path iprscan_results
    path eggnog_results
    val species

    output:
    path "anno/*", emit: annotation_results

    script:
    """
    funannotate annotate --gff ${gff} \
        --fasta ${genome} \
        -s "${species}" \
        --busco_db ${params.busco_db} \
        --eggnog ${eggnog_results} \
        --iprscan ${iprscan_results} \
        --phobius ${phobius_results} \
        --cpus ${params.threads} \
        -o anno
    """
}

// Workflow
workflow {
    if (params.run_mode in ["prediction_only", "both"]) {
        // Core processing
        trna_results = tRNAScan(params.masked_genome)

        def merged_bam = null
        def stringtie_results = null
        def transcript_fasta = null
        def training_data = null

        if (has_rnaseq) {
            cleaned_reads = CleanReads(read_pairs)
            genome_index = IndexGenome(params.genome_fasta)
            mapped_bams = MapReads(cleaned_reads.paired_reads, genome_index)
            merged_bam = MergeBAMs(mapped_bams.bam.collect())
            stringtie_results = StringTie(merged_bam)
            transcript_fasta = GTFtoFasta(stringtie_results.gtf, params.genome_fasta)
            training_data = FunannotateTrain(params.genome_fasta, cleaned_reads.paired_reads.flatten().collect(), params.species)
        }

        // Funannotate prediction
        funannotate_results = FunannotatePredict(
            params.masked_genome,
            merged_bam ?: [],
            stringtie_results ? stringtie_results.gtf : [],
            transcript_fasta ?: [],
            trna_results.trna_results,
            params.species
        )
        funannotate_updated = FunannotateUpdate(funannotate_results.prediction_results)

        // BRAKER3 prediction
        braker_results = BRAKER3(
            params.masked_genome,
            merged_bam ?: [],
            params.protein_evidence,
            params.species
        )

        // Add UTRs to BRAKER3 GTF if RNA-Seq is available
        def braker_utr = null
        if (has_rnaseq && braker_results.braker_gtf) {
            braker_utr = AddUTRs(braker_results.braker_gtf, stringtie_results.gtf)
        } else {
            braker_utr = braker_results.braker_gtf
        }

        // Convert tRNAscan-SE to GFF
        trna_gff = Convert_tRNA_to_GFF(trna_results.trna_out)

        // Convert BRAKER3 GTF to GFF3
        def braker_gff3 = null
        if (braker_utr) {
            braker_gff3 = Convert_GTF_to_GFF(braker_utr)
        } else {
            // Create a placeholder GFF3 if no GTF (minimal impact, as BUSCO will handle)
            file("${params.outdir}/braker3_utr/placeholder.gff3").text = ""
            braker_gff3 = file("${params.outdir}/braker3_utr/placeholder.gff3")
        }

        // Merge BRAKER3 GFF3 with tRNA GFF if valid
        def merged_gff = null
        if (braker_gff3.name != "placeholder.gff3") {
            merged_gff = MergeGFFs(braker_gff3, trna_gff)
        } else {
            merged_gff = trna_gff
        }

        // Fix overlapping genes
        fixed_gff = FixOverlappingGenes(merged_gff)

        // Validate GFF3
        validation_report = ValidateGFF(fixed_gff)

        // Export proteins from BRAKER3
        braker_proteins = ExportProteinsBRAKER(fixed_gff, params.genome_fasta)

        // Run BUSCO for quality control
        busco_funannotate = BUSCO_Funannotate(funannotate_updated.proteins)
        busco_braker3 = BUSCO_BRAKER3(braker_proteins.proteins)

        // Compare BUSCO scores
        compare_busco = CompareBUSCO(busco_funannotate.summary, busco_braker3.summary)
    }

    if (params.run_mode in ["annotation_only", "both"] && params.enable_functional_annotation) {
        // Select GFF and proteins based on annotation_prediction
        if (params.run_mode == "annotation_only") {
            if (params.annotation_prediction == "braker") {
                gff = file("${params.outdir}/merged_annotations/fixed_merged.gff")
                proteins = file("${params.outdir}/proteins/braker/braker_proteins.fasta")
                if (!gff.exists() || !proteins.exists()) {
                    error "BRAKER3 GFF or proteins not found for annotation_only mode"
                }
            } else {
                gff = file("${params.outdir}/funannotate_update/update_results/${params.species.replace(' ', '_')}.gff3")
                proteins = file("${params.outdir}/funannotate_update/update_results/${params.species.replace(' ', '_')}.proteins.fa")
                if (!gff.exists() || !proteins.exists()) {
                    error "Funannotate GFF or proteins not found for annotation_only mode"
                }
            }
        } else {
            // For 'both' mode, use user-specified or best BUSCO prediction
            if (params.annotation_prediction == "braker") {
                gff = fixed_gff
                proteins = braker_proteins.proteins
            } else if (params.annotation_prediction == "funannotate") {
                gff = funannotate_updated.gff3
                proteins = funannotate_updated.proteins
            } else {
                // Use best BUSCO prediction
                gff = compare_busco.best_prediction[1]
                proteins = compare_busco.best_prediction[2]
            }
        }

        // Functional annotation
        phobius_results = Phobius(proteins)
        iprscan_results = InterProScan(proteins)
        eggnog_results = EggNOGMapper(proteins)
        annotation_results = FunannotateAnnotate(
            gff,
            params.genome_fasta,
            proteins,
            phobius_results,
            iprscan_results,
            eggnog_results,
            params.species
        )
    }
}
