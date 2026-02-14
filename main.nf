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
    def map = [:]
    headers.eachWithIndex { h, i -> map[h] = (i < values.size() && values[i]) ? values[i] : null }
    return map
}

// =============== DUMMY FILES ===============
def hidden_dummy_dir = file("$projectDir/.dummies")
hidden_dummy_dir.mkdirs()
def DUMMY_PLUS_BAM  = file("$hidden_dummy_dir/NO_PLUS_BAM.bam")
def DUMMY_MINUS_BAM = file("$hidden_dummy_dir/NO_MINUS_BAM.bam")
if (!DUMMY_PLUS_BAM.exists())  DUMMY_PLUS_BAM.text = ''
if (!DUMMY_MINUS_BAM.exists()) DUMMY_MINUS_BAM.text = ''

def dummy_files = [
    "${workflow.workDir}/NO_GTF_FILE.gtf", "${workflow.workDir}/NO_BAM_FILE.bam",
    "${workflow.workDir}/NO_TRANSCRIPTS.fasta", "${workflow.workDir}/NO_R1.fastq.gz",
    "${workflow.workDir}/NO_R2.fastq.gz", "${workflow.workDir}/NO_PROTEINS.fa",
    "${workflow.workDir}/NO_GFF3.gff3", "${workflow.workDir}/NO_FUNANNO_DB.empty",
    "${workflow.workDir}/NO_EGGNOG_DB.empty", "${workflow.workDir}/NO_PHOBIUS_TARBALL.empty",
    "${workflow.workDir}/NO_SIGNALP_TARBALL.empty"
].collect { file(it, checkIfExists: false) }

dummy_files.each { f -> if (!f.exists()) { f.parent.mkdirs(); f.text = '' } }

// =============== MODULE INCLUDES ===============
include { TRNASCAN_SE } from './modules/trnascan_se.nf'
include { RNASEQ_PROCESSING } from './modules/rnaseq_processing.nf'
include { FUNANNOTATE } from './modules/funannotate.nf'
include { BRAKER_RUN } from './modules/braker_run.nf'
include { BRAKER_POST } from './modules/braker_post.nf'
include { COMPARE_BUSCO } from './modules/compare_busco.nf'
include { MERGE_ANNOTATIONS } from './modules/merge_annotations.nf'
include { FUNCTIONAL_ANNOTATION } from './modules/func_annotate.nf'

workflow {
    def mandatoryParams = parseCsvToMap(params.mandatory_csv, 'Mandatory')
    def optionalParams = parseCsvToMap(params.optional_csv, 'Optional')
    optionalParams.stranded = optionalParams.stranded ?: params.rnaseq_stranded
    
    def fullMeta = [
        id: mandatoryParams.name,
        species: mandatoryParams.species.replaceAll(" ", "_"),
        organism: mandatoryParams.organism,
        busco_db: mandatoryParams.busco_db,
        busco_db_fun: mandatoryParams.busco_db_fun ?: '',
        stranded: optionalParams.stranded,
        gc_probability: optionalParams.gc_probability ?: '',
        funanno_DB: optionalParams.funanno_DB ?: ''
    ]
    
    TRNASCAN_SE(Channel.of([fullMeta, file(mandatoryParams.genome_masked)]), Channel.of(file(params.script_dir)))
    def trna_gff_channel = TRNASCAN_SE.out.gff
    
    def has_rnaseq = optionalParams.rnaseq_dir && file(optionalParams.rnaseq_dir).exists()
    def rnaseq_outputs
    
    if (has_rnaseq) {
        def proc = RNASEQ_PROCESSING(Channel.of(fullMeta), Channel.of(file(mandatoryParams.genome_unmasked)), Channel.of(file(optionalParams.rnaseq_dir)), Channel.of(file(params.script_dir)), Channel.of(optionalParams.stranded))
        rnaseq_outputs = [
            gtf: proc.gtf, bam: proc.bam, transcripts: proc.transcripts,
            r1: proc.trimmed_r1, r2: proc.trimmed_r2,
            plus_bam: proc.plus_strand, minus_bam: proc.minus_strand
        ]
    } else {
        rnaseq_outputs = [
            gtf: Channel.of([fullMeta, dummy_files[0]]), bam: Channel.of([fullMeta, dummy_files[1]]),
            transcripts: Channel.of([fullMeta, dummy_files[2]]), r1: Channel.of([fullMeta, dummy_files[3]]),
            r2: Channel.of([fullMeta, dummy_files[4]]), plus_bam: Channel.of([fullMeta, DUMMY_PLUS_BAM]),
            minus_bam: Channel.of([fullMeta, DUMMY_MINUS_BAM])
        ]
    }

    def findGenemarkFiles = { dir ->
        def f = new File(dir).listFiles()
        def key = file(f.find{ it.name.startsWith('gm_key') && it.name.endsWith('.gz') }.path)
        def tar = file(f.find{ it.name.startsWith('gmes_linux') && it.name.endsWith('.tar.gz') }.path)
        return [key, tar]
    }
    def (gmKeyFile, gmesFile) = findGenemarkFiles(mandatoryParams.genemark_dir)

    if (params.mode == 'both' || params.mode == 'funannotate') {
        def funannotate_input = TRNASCAN_SE.out.highconf
            .combine(rnaseq_outputs.gtf).combine(rnaseq_outputs.bam).combine(rnaseq_outputs.transcripts).combine(rnaseq_outputs.r1).combine(rnaseq_outputs.r2)
            .map { tuple ->
                def meta = tuple[0]
                def safeF = { p -> (p && p.toString() != "" && !p.toString().contains('NO_') && file(p).exists()) ? file(p) : file("${workflow.workDir}/NO_FILE") }
                def nano = optionalParams.nanopore_mrna ? safeF(optionalParams.nanopore_mrna) : file("${workflow.workDir}/NO_FILE")
                def pb   = optionalParams.pacbio_isoseq ? safeF(optionalParams.pacbio_isoseq) : file("${workflow.workDir}/NO_FILE")
                [meta, file(mandatoryParams.genome_masked), file(mandatoryParams.genome_unmasked), tuple[1], file(mandatoryParams.protein_evidence), tuple[3], tuple[5], tuple[7], tuple[9], tuple[11], gmKeyFile, gmesFile, nano, pb]
            }
        FUNANNOTATE(funannotate_input)
    }

    if (params.mode == 'both' || params.mode == 'braker') {
        def braker_input = Channel.of(fullMeta)
            .combine(Channel.of(file(mandatoryParams.genome_masked)))
            .combine(Channel.of(file(mandatoryParams.protein_evidence)))
            .combine(rnaseq_outputs.bam.map { it[1] })
            .combine(fullMeta.stranded in ['forward', 'reverse'] ? rnaseq_outputs.plus_bam.map { it[1] } : Channel.of(DUMMY_PLUS_BAM))
            .combine(fullMeta.stranded in ['forward', 'reverse'] ? rnaseq_outputs.minus_bam.map { it[1] } : Channel.of(DUMMY_MINUS_BAM))
            .combine(Channel.of(fullMeta.gc_probability))
            .map { meta, genome, proteins, bam, plus, minus, gc -> [meta, genome, proteins, bam, plus, minus, gc] }
        BRAKER_RUN(braker_input)
        
        def braker_post_input = BRAKER_RUN.out.braker_gtf.map { m, g -> [m, g] }
            .combine(BRAKER_RUN.out.braker_dir.map { _, d -> d })
            .combine(Channel.of(file(mandatoryParams.genome_unmasked)))
            .combine(trna_gff_channel.map { _, g -> g })
            .combine(Channel.of(file(params.script_dir)))
            .map { meta, gtf, dir, genome, trna, script -> [meta, gtf, dir, genome, trna, script] }
        BRAKER_POST(braker_post_input)
    }

    // =============== DUMMY FALLBACKS ===============
    def no_prot = dummy_files[5]
    def no_gff = dummy_files[6]

    def fa_proteins = (params.mode in ['both', 'funannotate'] ? FUNANNOTATE.out.proteins : Channel.of([fullMeta, no_prot])).ifEmpty([fullMeta, no_prot])
    def fa_gff3 = (params.mode in ['both', 'funannotate'] ? FUNANNOTATE.out.gff3 : Channel.of([fullMeta, no_gff])).ifEmpty([fullMeta, no_gff])
    def br_proteins = (params.mode in ['both', 'braker'] ? BRAKER_POST.out.proteins : Channel.of([fullMeta, no_prot])).ifEmpty([fullMeta, no_prot])
    def br_gff3 = (params.mode in ['both', 'braker'] ? BRAKER_POST.out.gff3 : Channel.of([fullMeta, no_gff])).ifEmpty([fullMeta, no_gff])

    // =============== MERGE ANNOTATIONS ===============
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
            .map { meta, best_file_path -> [meta, file(best_file_path)] }
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
                def protein_fasta = no_prot
                def gff_file = no_gff
                if (params.mode == 'funannotate') {
                    protein_fasta = (fa_prot && file(fa_prot).exists() && file(fa_prot).size() > 0) ? file(fa_prot) : no_prot
                    gff_file = (fa_gff && file(fa_gff).exists() && file(fa_gff).size() > 0) ? file(fa_gff) : no_gff
                } else if (params.mode == 'braker' || params.mode == 'both') {
                    protein_fasta = (br_prot && file(br_prot).exists() && file(br_prot).size() > 0) ? file(br_prot) : no_prot
                    gff_file = (br_gff && file(br_gff).exists() && file(br_gff).size() > 0) ? file(br_gff) : no_gff
                }
                [meta, protein_fasta, gff_file, file(mandatoryParams.genome_unmasked),
                 optionalParams.funanno_DB ? file(optionalParams.funanno_DB) : dummy_files[7],
                 optionalParams.eggnog_DB ? file(optionalParams.eggnog_DB) : dummy_files[8],
                 gmKeyFile, gmesFile,
                 optionalParams.func_tool_dir ? file("${optionalParams.func_tool_dir}/phobius101_linux.tgz") : dummy_files[9],
                 optionalParams.func_tool_dir ? file("${optionalParams.func_tool_dir}/signalp-6.0h.fast.tar.gz") : dummy_files[10]]
            }
        FUNCTIONAL_ANNOTATION(functional_annotation_input)
    }
}

workflow.onComplete {
    log.info "Pipeline completed! Results in: ${params.outdir}"
}
