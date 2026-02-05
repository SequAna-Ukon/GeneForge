#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =============== PARAMETERS ===============
params.mandatory_csv = 'mandatory.csv'
params.optional_csv = 'optional.csv'
params.outdir = 'results'
params.script_dir = "$projectDir/scripts"
params.rnaseq_stranded = 'no'
params.mode = 'both'
params.func_annotation = false

// =============== VALIDATION ===============
def valid_modes = ['both', 'braker', 'funannotate']
if (!valid_modes.contains(params.mode)) {
    exit 1, "ERROR: Invalid mode '${params.mode}'. Valid modes are: ${valid_modes.join(', ')}"
}

// =============== CSV PARSING ===============
def parseCsvToMap(filePath, csvType = 'Mandatory') {
    if (!file(filePath).exists()) exit 1, "ERROR: ${csvType} CSV file ${filePath} not found!"
    def text = file(filePath).text.trim()
    def lines = text.split('\n')
    if (lines.size() < 2) exit 1, "ERROR: ${csvType} CSV file ${filePath} needs at least 2 lines"
    def headers = lines[0].split(',').collect { it.trim() }
    def values = lines[1].split(',', -1).collect { it.trim() }
    
    log.info "=" * 60
    log.info "${csvType} CSV Headers: ${headers}"
    log.info "${csvType} CSV Values:  ${values}"
    log.info "=" * 60
    
    def map = [:]
    headers.eachWithIndex { h, i -> map[h] = (i < values.size() && values[i]) ? values[i] : null }
    return map
}

// =============== DUMMY FILES FOR CACHING ===============
// Stable hidden dummies for strand-specific BAMS (ensures perfect caching)
def hidden_dummy_dir = file("$projectDir/.dummies")
hidden_dummy_dir.mkdirs()

def DUMMY_PLUS_BAM  = file("$hidden_dummy_dir/NO_PLUS_BAM.bam")
def DUMMY_MINUS_BAM = file("$hidden_dummy_dir/NO_MINUS_BAM.bam")

if (!DUMMY_PLUS_BAM.exists())  DUMMY_PLUS_BAM.text = ''
if (!DUMMY_MINUS_BAM.exists()) DUMMY_MINUS_BAM.text = ''

// Other workDir dummies (only used when no RNA-Seq)
def dummy_files = [
    "${workflow.workDir}/NO_GTF_FILE.gtf",
    "${workflow.workDir}/NO_BAM_FILE.bam",
    "${workflow.workDir}/NO_TRANSCRIPTS.fasta",
    "${workflow.workDir}/NO_R1.fastq.gz",
    "${workflow.workDir}/NO_R2.fastq.gz",
    "${workflow.workDir}/NO_FA_BUSCO.txt",
    "${workflow.workDir}/NO_BR_BUSCO.txt",
    "${workflow.workDir}/NO_PROTEINS.fa",
    "${workflow.workDir}/NO_GFF3.gff3",
    "${workflow.workDir}/NO_FUNANNO_DB.empty",
    "${workflow.workDir}/NO_EGGNOG_DB.empty",
    "${workflow.workDir}/NO_PHOBIUS_TARBALL.empty",
    "${workflow.workDir}/NO_SIGNALP_TARBALL.empty"
].collect { file(it, checkIfExists: false) }

dummy_files.each { f ->
    if (!f.exists()) {
        f.parent.mkdirs()
        f.text = ''
    }
}

// =============== MODULE INCLUDES ===============
include { TRNASCAN_SE } from './modules/trnascan_se.nf'
include { RNASEQ_PROCESSING } from './modules/rnaseq_processing.nf'
include { FUNANNOTATE } from './modules/funannotate.nf'
include { BRAKER_RUN } from './modules/braker_run.nf'
include { BRAKER_POST } from './modules/braker_post.nf'
include { COMPARE_BUSCO } from './modules/compare_busco.nf'
include { MERGE_ANNOTATIONS } from './modules/merge_annotations.nf'
include { FUNCTIONAL_ANNOTATION } from './modules/func_annotate.nf'

// =============== MAIN WORKFLOW ===============
workflow {
    log.info "Starting GeneForge Pipeline"
    log.info "Mode: ${params.mode}"
    log.info "Output directory: ${params.outdir}"
    
    // =============== PARSE CONFIGURATION FILES ===============
    def mandatoryParams = parseCsvToMap(params.mandatory_csv, 'Mandatory')
    def optionalParams = parseCsvToMap(params.optional_csv, 'Optional')
    
    // =============== VALIDATE MANDATORY INPUTS ===============
    if (!file(mandatoryParams.genome_masked).exists()) exit 1, "ERROR: genome_masked missing"
    if (!file(mandatoryParams.genome_unmasked).exists()) exit 1, "ERROR: genome_unmasked missing"
    if (!file(mandatoryParams.protein_evidence).exists()) exit 1, "ERROR: protein_evidence missing"
    if (!file(mandatoryParams.genemark_dir).exists()) exit 1, "ERROR: genemark_dir missing"
    
    // =============== PROCESS OPTIONAL PARAMETERS ===============
    def valid_stranded = ['no', 'forward', 'reverse']
    if (optionalParams.stranded && !(optionalParams.stranded in valid_stranded))
        exit 1, "Invalid stranded value: ${optionalParams.stranded}. Valid values: ${valid_stranded.join(', ')}"
    optionalParams.stranded = optionalParams.stranded ?: params.rnaseq_stranded
    
    // =============== CREATE METADATA ===============
    def fullMeta = [
        id: mandatoryParams.name,
        species: mandatoryParams.species.replaceAll(" ", "_"),
        organism: mandatoryParams.organism,
        busco_db: mandatoryParams.busco_db,
        busco_db_fun: mandatoryParams.busco_db_fun ?: '',
        stranded: optionalParams.stranded,
        gc_probability: optionalParams.gc_probability ?: ''
    ]
    
    // =============== tRNA SCAN ===============
    TRNASCAN_SE(Channel.of([fullMeta, file(mandatoryParams.genome_masked)]), Channel.of(file(params.script_dir)))
    def trna_gff_channel = TRNASCAN_SE.out.gff
    
    // =============== RNA-SEQ PROCESSING ===============
    def has_rnaseq = optionalParams.rnaseq_dir && file(optionalParams.rnaseq_dir).exists()
    def rnaseq_outputs
    
    if (has_rnaseq) {
        def rnaseq_process = RNASEQ_PROCESSING(
            Channel.of(fullMeta),
            Channel.of(file(mandatoryParams.genome_unmasked)),
            Channel.of(file(optionalParams.rnaseq_dir)),
            Channel.of(file(params.script_dir)),
            Channel.of(optionalParams.stranded)
        )
        rnaseq_outputs = [
            gtf: rnaseq_process.gtf,
            bam: rnaseq_process.bam,
            transcripts: rnaseq_process.transcripts,
            r1: rnaseq_process.trimmed_r1,
            r2: rnaseq_process.trimmed_r2,
            plus_bam: rnaseq_process.plus_strand,
            minus_bam: rnaseq_process.minus_strand
        ]
    } else {
        rnaseq_outputs = [
            gtf: Channel.of([fullMeta, file("${workflow.workDir}/NO_GTF_FILE.gtf")]),
            bam: Channel.of([fullMeta, file("${workflow.workDir}/NO_BAM_FILE.bam")]),
            transcripts: Channel.of([fullMeta, file("${workflow.workDir}/NO_TRANSCRIPTS.fasta")]),
            r1: Channel.of([fullMeta, file("${workflow.workDir}/NO_R1.fastq.gz")]),
            r2: Channel.of([fullMeta, file("${workflow.workDir}/NO_R2.fastq.gz")]),
            plus_bam: Channel.of([fullMeta, DUMMY_PLUS_BAM]),
            minus_bam: Channel.of([fullMeta, DUMMY_MINUS_BAM])
        ]
    }
    
    // =============== FIND GENEMARK FILES ===============
    def findGenemarkFiles = { genemarkDir ->
        // Find gm_key*.gz file
        def gmKeyFiles = new File(genemarkDir).listFiles({ f -> 
            f.name.startsWith('gm_key') && f.name.endsWith('.gz')
        } as FileFilter)
        
        if (!gmKeyFiles || gmKeyFiles.size() == 0) {
            exit 1, "ERROR: No gm_key*.gz file found in ${genemarkDir}"
        }
        
        def gmKeyFile = file(gmKeyFiles[0].path)
        
        // Find gmes_linux*.tar.gz file
        def gmesFiles = new File(genemarkDir).listFiles({ f -> 
            f.name.startsWith('gmes_linux') && f.name.endsWith('.tar.gz')
        } as FileFilter)
        
        if (!gmesFiles || gmesFiles.size() == 0) {
            exit 1, "ERROR: No gmes_linux*.tar.gz file found in ${genemarkDir}"
        }
        
        def gmesFile = file(gmesFiles[0].path)
        
        return [gmKeyFile, gmesFile]
    }

    // Get Genemark files once
    def (gmKeyFile, gmesFile) = findGenemarkFiles(mandatoryParams.genemark_dir)
    
    // =============== FUNANNOTATE ===============
    if (params.mode == 'both' || params.mode == 'funannotate') {
        def funannotate_input = TRNASCAN_SE.out.highconf
            .combine(rnaseq_outputs.gtf)
            .combine(rnaseq_outputs.bam)
            .combine(rnaseq_outputs.transcripts)
            .combine(rnaseq_outputs.r1)
            .combine(rnaseq_outputs.r2)
            .map { tuple ->
                def meta = tuple[0]
                def safeFile = { p -> p && file(p).exists() ? file(p) : file("${workflow.workDir}/NO_FILE") }
                [meta, file(mandatoryParams.genome_masked), file(mandatoryParams.genome_unmasked), tuple[1],
                 file(mandatoryParams.protein_evidence), safeFile(tuple[3]), safeFile(tuple[5]), safeFile(tuple[7]),
                 safeFile(tuple[9]), safeFile(tuple[11]), gmKeyFile, gmesFile, '', '']
            }

        FUNANNOTATE(funannotate_input)
    }
    
    // =============== BRAKER ===============
    if (params.mode == 'both' || params.mode == 'braker') {
        // BRAKER — stable dummies for perfect caching
        def braker_input = Channel.of(fullMeta)
            .combine(Channel.of(file(mandatoryParams.genome_masked)))
            .combine(Channel.of(file(mandatoryParams.protein_evidence)))
            .combine(rnaseq_outputs.bam.map { it[1] })
            .combine(fullMeta.stranded in ['forward', 'reverse'] ? rnaseq_outputs.plus_bam.map { it[1] } : Channel.of(DUMMY_PLUS_BAM))
            .combine(fullMeta.stranded in ['forward', 'reverse'] ? rnaseq_outputs.minus_bam.map { it[1] } : Channel.of(DUMMY_MINUS_BAM))
            .combine(Channel.of(fullMeta.gc_probability))
            .map { meta, genome, proteins, bam, plus_bam, minus_bam, gc_prob ->
                [meta, genome, proteins, bam, plus_bam, minus_bam, gc_prob]
            }

        BRAKER_RUN(braker_input)
        
        if (params.mode == 'both' || params.mode == 'braker') {
            def braker_post_input = BRAKER_RUN.out.braker_gtf.map { meta, gtf -> [meta, gtf] }
                .combine(BRAKER_RUN.out.braker_dir.map { _, dir -> dir })
                .combine(Channel.of(file(mandatoryParams.genome_unmasked)))
                .combine(trna_gff_channel.map { _, gff -> gff })
                .combine(Channel.of(file(params.script_dir)))
                .map { meta, gtf, dir, genome, trna, script -> [meta, gtf, dir, genome, trna, script] }
            BRAKER_POST(braker_post_input)
        }
    }
    
    // =============== DUMMY FALLBACKS FOR OUTPUTS ===============
    def no_fa_proteins = file("${workflow.workDir}/NO_PROTEINS.fa")
    def no_fa_gff = file("${workflow.workDir}/NO_GFF3.gff3")
    def no_br_proteins = file("${workflow.workDir}/NO_PROTEINS.fa")
    def no_br_gff = file("${workflow.workDir}/NO_GFF3.gff3")

    def fa_proteins = (params.mode in ['both', 'funannotate'] ? FUNANNOTATE.out.proteins : Channel.of([fullMeta, no_fa_proteins])).ifEmpty([fullMeta, no_fa_proteins])
    def fa_gff3 = (params.mode in ['both', 'funannotate'] ? FUNANNOTATE.out.gff3 : Channel.of([fullMeta, no_fa_gff])).ifEmpty([fullMeta, no_fa_gff])
    def br_proteins = (params.mode in ['both', 'braker'] ? BRAKER_POST.out.proteins : Channel.of([fullMeta, no_br_proteins])).ifEmpty([fullMeta, no_br_proteins])
    def br_gff3 = (params.mode in ['both', 'braker'] ? BRAKER_POST.out.gff3 : Channel.of([fullMeta, no_br_gff])).ifEmpty([fullMeta, no_br_gff])
    
    // =============== MERGE ANNOTATIONS (BOTH MODES) ===============
    if (params.mode == 'both') {
        def compare_input = FUNANNOTATE.out.busco_summary
            .combine(BRAKER_POST.out.busco_summary, by: 0)
            .combine(FUNANNOTATE.out.gff3, by: 0)
            .combine(FUNANNOTATE.out.proteins, by: 0)
            .combine(BRAKER_POST.out.gff3, by: 0)
            .combine(BRAKER_POST.out.proteins, by: 0)
            .map { meta, fa_busco, br_busco, fa_gff, fa_prot, br_gff, br_prot ->
                [meta, fa_busco, br_busco, fa_prot, fa_gff, br_prot, br_gff]
            }
        
        COMPARE_BUSCO(compare_input, Channel.of(file(params.script_dir)))
        
        def merge_input = COMPARE_BUSCO.out.best_choice
            .map { meta, best_file_path ->
                [meta, file(best_file_path)]
            }
            .combine(FUNANNOTATE.out.gff3, by: 0)
            .combine(FUNANNOTATE.out.proteins, by: 0)
            .combine(BRAKER_POST.out.gff3, by: 0)
            .combine(BRAKER_POST.out.proteins, by: 0)
            .combine(trna_gff_channel, by: 0)
            .combine(Channel.of([fullMeta, file(mandatoryParams.genome_unmasked)]), by: 0)
            .map { meta, best_file, fa_gff, fa_prot, br_gff, br_prot, trna_gff_file, genome ->
                [meta, best_file, fa_gff, fa_prot, br_gff, br_prot, trna_gff_file, genome]
            }
        
        MERGE_ANNOTATIONS(merge_input)
    }
    
    // =============== FUNCTIONAL ANNOTATION ===============
    if (params.func_annotation) {
        def functional_annotation_input = fa_gff3
            .combine(fa_proteins, by: 0)
            .combine(br_gff3, by: 0)
            .combine(br_proteins, by: 0)
            .map { meta, fa_gff, fa_prot, br_gff, br_prot ->
                def protein_fasta = no_br_proteins
                def gff_file = no_br_gff
                if (params.mode == 'funannotate') {
                    protein_fasta = file(fa_prot).exists() && file(fa_prot).size() > 0 ? file(fa_prot) : no_fa_proteins
                    gff_file = file(fa_gff).exists() && file(fa_gff).size() > 0 ? file(fa_gff) : no_fa_gff
                } else if (params.mode == 'braker' || params.mode == 'both') {
                    protein_fasta = file(br_prot).exists() && file(br_prot).size() > 0 ? file(br_prot) : no_br_proteins
                    gff_file = file(br_gff).exists() && file(br_gff).size() > 0 ? file(br_gff) : no_br_gff
                }
                [meta, protein_fasta, gff_file, file(mandatoryParams.genome_unmasked),
                 optionalParams.funanno_DB ? file(optionalParams.funanno_DB) : file("${workflow.workDir}/NO_FUNANNO_DB.empty"),
                 optionalParams.eggnog_DB ? file(optionalParams.eggnog_DB) : file("${workflow.workDir}/NO_EGGNOG_DB.empty"),
                 gmKeyFile,
                 gmesFile,
                 optionalParams.func_tool_dir ? file("${optionalParams.func_tool_dir}/phobius101_linux.tgz") : file("${workflow.workDir}/NO_PHOBIUS_TARBALL.empty"),
                 optionalParams.func_tool_dir ? file("${optionalParams.func_tool_dir}/signalp-6.0h.fast.tar.gz") : file("${workflow.workDir}/NO_SIGNALP_TARBALL.empty")]
            }
        
        FUNCTIONAL_ANNOTATION(functional_annotation_input)
    }
}

// =============== CLEANUP ===============
workflow.onComplete {
    log.info "Pipeline completed successfully!"
    dummy_files.each { f -> if (f.exists()) f.delete() }
    log.info "Results available in: ${params.outdir}"
}
