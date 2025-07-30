#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.mandatory_csv = 'mandatory.csv'
params.optional_csv = 'optional.csv'
params.outdir = 'results'
params.script_dir = "$projectDir/scripts"
params.rnaseq_stranded = 'no'
params.mode = 'both' // Options: both, braker, funannotate
params.func_annotation = false // Control FUNCTIONAL_ANNOTATION

// Validate mode parameter
def valid_modes = ['both', 'braker', 'funannotate']
if (!valid_modes.contains(params.mode)) {
    exit 1, "ERROR: Invalid mode '${params.mode}'. Valid modes are: ${valid_modes.join(', ')}"
}

// Utility function to parse CSV files
def parseCsvToMap(filePath) {
    if (!file(filePath).exists()) {
        exit 1, "ERROR: CSV file ${filePath} not found!"
    }

    def text = file(filePath).text.trim()
    def lines = text.split('\n')
    if (lines.size() < 2) {
        exit 1, "ERROR: CSV file ${filePath} needs at least 2 lines (header + data)"
    }

    def headers = lines[0].split(',').collect { it.trim() }
    def values = lines[1].split(',', -1).collect { it.trim() }

    println "CSV Headers: $headers"
    println "CSV Values: $values"

    def map = [:]
    headers.eachWithIndex { h, i -> map[h] = (i < values.size() && values[i]) ? values[i] : null }
    return map
}

// Create dummy files in the Nextflow work directory
def dummy_files = [
    "${workflow.workDir}/NO_GTF_FILE.gtf",
    "${workflow.workDir}/NO_BAM_FILE.bam",
    "${workflow.workDir}/NO_TRANSCRIPTS.fasta",
    "${workflow.workDir}/NO_R1.fastq.gz",
    "${workflow.workDir}/NO_R2.fastq.gz",
    "${workflow.workDir}/NO_PLUS_BAM_FILE.bam",
    "${workflow.workDir}/NO_MINUS_BAM_FILE.bam",
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
        f.parent.mkdirs() // Ensure parent directory exists in workDir
        f.text = ''
    }
}

// Include modules
include { TRNASCAN_SE }   from './modules/trnascan_se.nf'
include { RNASEQ_PROCESSING } from './modules/rnaseq_processing.nf'
include { FUNANNOTATE }   from './modules/funannotate.nf'
include { BRAKER_RUN }    from './modules/braker_run.nf'
include { BRAKER_POST }   from './modules/braker_post.nf'
include { COMPARE_BUSCO } from './modules/compare_busco.nf'
include { FUNCTIONAL_ANNOTATION } from './modules/func_annotate.nf'

workflow {
    def mandatoryParams = parseCsvToMap(params.mandatory_csv)
    def optionalParams = parseCsvToMap(params.optional_csv)

    // Validate mandatory inputs
    if (!file(mandatoryParams.genome_masked).exists()) exit 1, "ERROR: genome_masked missing: ${mandatoryParams.genome_masked}"
    if (!file(mandatoryParams.genome_unmasked).exists()) exit 1, "ERROR: genome_unmasked missing: ${mandatoryParams.genome_unmasked}"
    if (!file(mandatoryParams.protein_evidence).exists()) exit 1, "ERROR: protein_evidence missing: ${mandatoryParams.protein_evidence}"
    if (!file(mandatoryParams.genemark_dir).exists()) exit 1, "ERROR: genemark_dir missing: ${mandatoryParams.genemark_dir}"
    if (!file("${mandatoryParams.genemark_dir}/gm_key_64.gz").exists()) exit 1, "ERROR: GeneMark key missing: ${mandatoryParams.genemark_dir}/gm_key_64.gz"
    if (!file("${mandatoryParams.genemark_dir}/gmes_linux_64_4.tar.gz").exists()) exit 1, "ERROR: GeneMark tar missing: ${mandatoryParams.genemark_dir}/gmes_linux_64_4.tar.gz"

    // Validate optional inputs
    if (optionalParams.funanno_DB && !file(optionalParams.funanno_DB).exists()) exit 1, "ERROR: Funannotate DB missing: ${optionalParams.funanno_DB}"
    if (optionalParams.eggnog_DB && !file(optionalParams.eggnog_DB).exists()) exit 1, "ERROR: EggNOG DB missing: ${optionalParams.eggnog_DB}"
    if (optionalParams.func_tool_dir && !file(optionalParams.func_tool_dir).exists()) exit 1, "ERROR: func_tool_dir missing: ${optionalParams.func_tool_dir}"
    if (optionalParams.func_tool_dir && !file("${optionalParams.func_tool_dir}/phobius101_linux.tgz").exists()) exit 1, "ERROR: Phobius tarball missing: ${optionalParams.func_tool_dir}/phobius101_linux.tgz"
    if (optionalParams.func_tool_dir && !file("${optionalParams.func_tool_dir}/signalp-6.0h.fast.tar.gz").exists()) exit 1, "ERROR: SignalP tarball missing: ${optionalParams.func_tool_dir}/signalp-6.0h.fast.tar.gz"

    // Validate stranded
    def valid_stranded = ['no', 'forward', 'reverse']
    if (optionalParams.stranded && !(optionalParams.stranded in valid_stranded))
        exit 1, "Invalid stranded value '${optionalParams.stranded}'"

    // Set defaults
    optionalParams.stranded        = optionalParams.stranded ?: params.rnaseq_stranded
    optionalParams.nanopore_mrna  = optionalParams.nanopore_mrna ?: ''
    optionalParams.pacbio_isoseq  = optionalParams.pacbio_isoseq ?: ''
    optionalParams.gc_probability = optionalParams.gc_probability ?: ''
    optionalParams.func_tool_dir  = optionalParams.func_tool_dir ?: ''
    optionalParams.eggnog_DB      = optionalParams.eggnog_DB ?: ''
    def busco_db_fun = mandatoryParams.busco_db_fun ?: ''

    // Define metadata
    def fullMeta = [
        id: mandatoryParams.name,
        species: mandatoryParams.species.replaceAll(" ", "_"),
        organism: mandatoryParams.organism,
        busco_db: mandatoryParams.busco_db,
        busco_db_fun: busco_db_fun,
        funanno_DB: optionalParams.funanno_DB ?: '',
        eggnog_DB: optionalParams.eggnog_DB ?: '',
        stranded: optionalParams.stranded,
        use_dual_bams: optionalParams.stranded in ['forward', 'reverse'],
        gc_probability: optionalParams.gc_probability ?: ''
    ]

    // TRNASCAN_SE
    def trna_input = Channel.of([fullMeta, file(mandatoryParams.genome_masked)])
    def trnascan_scripts = Channel.of(file(params.script_dir))
    TRNASCAN_SE(trna_input, trnascan_scripts)
    def highconf_trna = TRNASCAN_SE.out.highconf // For FUNANNOTATE and BRAKER_RUN (.tbl)
    def trna_gff = TRNASCAN_SE.out.gff // For BRAKER_POST (.gff)

    // RNA-seq processing
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
            plus_bam: Channel.of([fullMeta, file("${workflow.workDir}/NO_PLUS_BAM_FILE.bam")]),
            minus_bam: Channel.of([fullMeta, file("${workflow.workDir}/NO_MINUS_BAM_FILE.bam")])
        ]
    }

    // Combine inputs for Funannotate
    def funannotate_input = highconf_trna
        .combine(rnaseq_outputs.gtf)
        .combine(rnaseq_outputs.bam)
        .combine(rnaseq_outputs.transcripts)
        .combine(rnaseq_outputs.r1)
        .combine(rnaseq_outputs.r2)
        .map { tuple ->
            def meta = tuple[0]
            def safeFile = { pathOrNull ->
                if (pathOrNull && file(pathOrNull).exists()) return file(pathOrNull)
                else return file("${workflow.workDir}/NO_${pathOrNull?.name ?: 'FILE'}")
            }

            def gtf = safeFile(tuple[3])
            def bam = safeFile(tuple[5])
            def transcripts = safeFile(tuple[7])
            def r1 = safeFile(tuple[9])
            def r2 = safeFile(tuple[11])

            if (!meta.id)       exit 1, "ERROR: 'name' is missing in mandatory.csv"
            if (!meta.species)  exit 1, "ERROR: 'species' is missing in mandatory.csv"
            if (!meta.organism) exit 1, "ERROR: 'organism' is missing in mandatory.csv"
            if (!meta.busco_db) exit 1, "ERROR: 'busco_db' is missing in mandatory.csv"

            return [
                meta,
                file(mandatoryParams.genome_masked),
                file(mandatoryParams.genome_unmasked),
                tuple[1], // trna
                file(mandatoryParams.protein_evidence),
                gtf,
                bam,
                transcripts,
                r1,
                r2,
                file("${mandatoryParams.genemark_dir}/gm_key_64.gz"),
                file("${mandatoryParams.genemark_dir}/gmes_linux_64_4.tar.gz"),
                optionalParams.nanopore_mrna ? file(optionalParams.nanopore_mrna) : '',
                optionalParams.pacbio_isoseq ? file(optionalParams.pacbio_isoseq) : ''
            ]
        }

    if (params.mode == 'both' || params.mode == 'funannotate') {
        FUNANNOTATE(funannotate_input)
    }

    // Combine inputs for BRAKER_RUN
    def braker_input = Channel.of(fullMeta)
        .combine(Channel.of(file(mandatoryParams.genome_masked)))
        .combine(Channel.of(file(mandatoryParams.protein_evidence)))
        .combine(rnaseq_outputs.bam.map { it[1] ?: file("${workflow.workDir}/NO_BAM_FILE.bam") })
        .combine(fullMeta.stranded == 'no' ? Channel.of(file("${workflow.workDir}/NO_PLUS_BAM_FILE.bam")) : rnaseq_outputs.plus_bam.map { it[1] ?: file("${workflow.workDir}/NO_PLUS_BAM_FILE.bam") })
        .combine(fullMeta.stranded == 'no' ? Channel.of(file("${workflow.workDir}/NO_MINUS_BAM_FILE.bam")) : rnaseq_outputs.minus_bam.map { it[1] ?: file("${workflow.workDir}/NO_MINUS_BAM_FILE.bam") })
        .combine(Channel.of(optionalParams.gc_probability))
        .map { tuple ->
            def meta = tuple[0]
            def genome = tuple[1]
            def protein_evidence = file(tuple[2])
            def bam = file(tuple[3])
            def plus_bam = file(tuple[4])
            def minus_bam = file(tuple[5])
            def gc_prob = tuple[6]

            return [meta, genome, protein_evidence, bam, plus_bam, minus_bam, gc_prob]
        }

    if (params.mode == 'both' || params.mode == 'braker') {
        BRAKER_RUN(braker_input)
    }

    // Combine inputs for BRAKER_POST
    if (params.mode == 'both' || params.mode == 'braker') {
        def braker_post_input = BRAKER_RUN.out.braker_gtf
            .map { meta, gtf -> [meta, gtf] }
            .combine(BRAKER_RUN.out.braker_dir.map { _, dir -> dir })
            .combine(Channel.of(file(mandatoryParams.genome_unmasked)))
            .combine(trna_gff.map { _, gff -> gff })
            .combine(Channel.of(file(params.script_dir)))
            .map { tuple ->
                def meta   = tuple[0]
                def gtf    = tuple[1]
                def dir    = tuple[2]
                def genome = tuple[3]
                def trna   = tuple[4]
                def script = tuple[5]

                if (!meta?.id) exit 1, "BRAKER_POST ERROR: meta is missing"
                if (!gtf || !file(gtf).exists()) exit 1, "BRAKER_POST ERROR: gtf file missing"
                if (!dir || !file(dir).exists()) exit 1, "BRAKER_POST ERROR: braker_dir missing"
                if (!genome || !file(genome).exists()) exit 1, "BRAKER_POST ERROR: genome_unmasked missing"
                if (!trna || !file(trna).exists()) exit 1, "BRAKER_POST ERROR: trna gff missing"
                if (!script || !file(script).exists()) exit 1, "BRAKER_POST ERROR: script_dir missing"

                return [meta, gtf, dir, genome, trna, script]
            }

        BRAKER_POST(braker_post_input)
    }

    // COMPARE_BUSCO and FUNCTIONAL_ANNOTATION
    if (params.func_annotation) {
        def no_fa_busco = file("${workflow.workDir}/NO_FA_BUSCO.txt")
        def no_br_busco = file("${workflow.workDir}/NO_BR_BUSCO.txt")
        def no_fa_proteins = file("${workflow.workDir}/NO_PROTEINS.fa")
        def no_fa_gff = file("${workflow.workDir}/NO_GFF3.gff3")
        def no_br_proteins = file("${workflow.workDir}/NO_PROTEINS.fa")
        def no_br_gff = file("${workflow.workDir}/NO_GFF3.gff3")
        def no_funanno_db = file("${workflow.workDir}/NO_FUNANNO_DB.empty")
        def no_eggnog_db = file("${workflow.workDir}/NO_EGGNOG_DB.empty")
        def no_phobius_tarball = file("${workflow.workDir}/NO_PHOBIUS_TARBALL.empty")
        def no_signalp_tarball = file("${workflow.workDir}/NO_SIGNALP_TARBALL.empty")

        // Collect Funannotate outputs with fallbacks
        def fa_busco_summary = (params.mode in ['both', 'funannotate'] ?
            FUNANNOTATE.out.busco_summary : Channel.of([fullMeta, no_fa_busco]))
            .ifEmpty([fullMeta, no_fa_busco])
        def fa_proteins = (params.mode in ['both', 'funannotate'] ?
            FUNANNOTATE.out.proteins : Channel.of([fullMeta, no_fa_proteins]))
            .ifEmpty([fullMeta, no_fa_proteins])
        def fa_gff3 = (params.mode in ['both', 'funannotate'] ?
            FUNANNOTATE.out.gff3 : Channel.of([fullMeta, no_fa_gff]))
            .ifEmpty([fullMeta, no_fa_gff])

        // Collect Braker outputs with fallbacks
        def br_busco_summary = (params.mode in ['both', 'braker'] ?
            BRAKER_POST.out.busco_summary : Channel.of([fullMeta, no_br_busco]))
            .ifEmpty([fullMeta, no_br_busco])
        def br_proteins = (params.mode in ['both', 'braker'] ?
            BRAKER_POST.out.proteins : Channel.of([fullMeta, no_br_proteins]))
            .ifEmpty([fullMeta, no_br_proteins])
        def br_gff3 = (params.mode in ['both', 'braker'] ?
            BRAKER_POST.out.gff3 : Channel.of([fullMeta, no_br_gff]))
            .ifEmpty([fullMeta, no_br_gff])

        // Combine inputs for COMPARE_BUSCO
        def compare_busco_input = fa_busco_summary
            .combine(br_busco_summary, by: 0)
            .combine(fa_proteins, by: 0)
            .combine(fa_gff3, by: 0)
            .combine(br_proteins, by: 0)
            .combine(br_gff3, by: 0)
            .map { meta, fa_busco, br_busco, fa_prot, fa_gff, br_prot, br_gff ->
                if (!meta) exit 1, "ERROR: meta is null in COMPARE_BUSCO input"
                [meta, fa_busco, br_busco, fa_prot, fa_gff, br_prot, br_gff]
            }

        // Run COMPARE_BUSCO
        if (params.mode == 'both') {
            COMPARE_BUSCO(compare_busco_input, Channel.of(file(params.script_dir)))
        }

        // Prepare FUNCTIONAL_ANNOTATION inputs
        def functional_annotation_input = fa_gff3
            .combine(fa_proteins, by: 0)
            .combine(br_gff3, by: 0)
            .combine(br_proteins, by: 0)
            .map { meta, fa_gff, fa_prot, br_gff, br_prot ->
                if (!meta) exit 1, "ERROR: meta is null in FUNCTIONAL_ANNOTATION input"

                def protein_fasta = no_br_proteins
                def gff_file = no_br_gff
                if (params.mode == 'funannotate') {
                    protein_fasta = (fa_prot && file(fa_prot).exists() && file(fa_prot).size() > 0) ? 
                        file(fa_prot) : no_fa_proteins
                    gff_file = (fa_gff && file(fa_gff).exists() && file(fa_gff).size() > 0) ? 
                        file(fa_gff) : no_fa_gff
                } else if (params.mode == 'braker' || params.mode == 'both') {
                    protein_fasta = (br_prot && file(br_prot).exists() && file(br_prot).size() > 0) ? 
                        file(br_prot) : no_br_proteins
                    gff_file = (br_gff && file(br_gff).exists() && file(br_gff).size() > 0) ? 
                        file(br_gff) : no_br_gff
                }

                [meta,
                 protein_fasta,
                 gff_file,
                 file(mandatoryParams.genome_unmasked),
                 optionalParams.funanno_DB ? file(optionalParams.funanno_DB) : no_funanno_db,
                 optionalParams.eggnog_DB ? file(optionalParams.eggnog_DB) : no_eggnog_db,
                 file("${mandatoryParams.genemark_dir}/gm_key_64.gz"),
                 file("${mandatoryParams.genemark_dir}/gmes_linux_64_4.tar.gz"),
                 optionalParams.func_tool_dir ? file("${optionalParams.func_tool_dir}/phobius101_linux.tgz") : no_phobius_tarball,
                 optionalParams.func_tool_dir ? file("${optionalParams.func_tool_dir}/signalp-6.0h.fast.tar.gz") : no_signalp_tarball]
            }

        // Run FUNCTIONAL_ANNOTATION
        FUNCTIONAL_ANNOTATION(functional_annotation_input)
    }
}

// Cleanup dummy files (optional, as they are in workDir and can be handled by Nextflow cleanup)
workflow.onComplete {
    dummy_files.each { f ->
        if (f.exists()) {
            f.delete()
        }
    }
}
