process BRAKER {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/braker", mode: 'copy', pattern: '*.braker_error.log'
    publishDir "${params.outdir}/braker", mode: 'copy', pattern: '*_braker.gff3'
    publishDir "${params.outdir}/braker", mode: 'copy', pattern: '*.braker.prot.fasta'
    publishDir "${params.outdir}/braker", mode: 'copy', pattern: '*_busco_braker.txt'
    container 'abdoallahsharaf/geneforge-braker3:2.1'

    input:
    tuple val(meta),
          path(genome_masked),
          path(protein_evidence),
          path(bam),
          path(plus_bam),
          path(minus_bam),
          val(gc_probability)
    
    path(genome_unmasked)
    path(trna_gff)
    path(script_dir)

    output:
    tuple val(meta), path("braker.gtf"),                               emit: braker_gtf
    tuple val(meta), path("braker"),                                   emit: braker_dir
    tuple val(meta), path("${meta.id}_braker.gff3"),                   emit: gff3
    tuple val(meta), path("${meta.id}.braker.prot.fasta"),             emit: proteins
    tuple val(meta), path("${meta.id}_busco_braker.txt"),              emit: busco_summary
    path "${meta.id}_braker_error.log",                               emit: error_log

    script:
    def prefix       = meta.id
    def species_base = meta.species.replaceAll(/\s+/, '_')
    def species      = "${species_base}_braker"
    def busco_db      = meta.busco_db

    // BAM placeholders to recognize as missing files
    def bam_placeholders       = ['NO_BAM_FILE.bam',      '', null]
    def plus_bam_placeholders  = ['NO_PLUS_BAM_FILE.bam', '', null]
    def minus_bam_placeholders = ['NO_MINUS_BAM_FILE.bam','', null]

    def bam_available      = (bam      != null && !bam_placeholders.contains(bam.getName()))
    def dual_bam_available = (plus_bam != null && minus_bam != null &&
                              !plus_bam_placeholders.contains(plus_bam.getName()) &&
                              !minus_bam_placeholders.contains(minus_bam.getName()))

    def bam_flag      = ''
    def stranded_flag = ''
    if (meta.stranded in ['forward', 'reverse'] && dual_bam_available) {
        bam_flag      = "--bam=${plus_bam},${minus_bam}"
        stranded_flag = "--stranded=+,-"
    } else if (bam_available) {
        bam_flag      = "--bam=${bam}"
        stranded_flag = ''
    } else {
        bam_flag      = ''
        stranded_flag = ''
    }
    def gc_flag = (gc_probability && gc_probability != '') ? "--gc_probability=${gc_probability}" : ""

    """
    #!/bin/bash
    set -euo pipefail

    # ---------------------------------------------------------------
    # Restore Clean Environment Variables (Keeps ETP/GeneMark Happy)
    # ---------------------------------------------------------------
    export PATH=/opt/BRAKER/scripts:/opt/Augustus/bin:/opt/Augustus/scripts:/opt/conda/bin:/opt/ETP/bin:/opt/ETP/bin/gmes:/opt/ETP/bin/gmes/ProtHint/bin:/opt/ETP/tools:/usr/local/bin:/usr/bin:/bin:\$PATH
    export AUGUSTUS_CONFIG_PATH=/opt/Augustus/config
    export AUGUSTUS_BIN_PATH=/opt/Augustus/bin
    export AUGUSTUS_SCRIPTS_PATH=/opt/Augustus/scripts
    export GENEMARK_PATH=/opt/ETP/bin/gmes
    export DIAMOND_PATH=/opt/ETP/tools
    export SAMTOOLS_PATH=/opt/ETP/tools
    export COMPLEASM_PATH=/opt/compleasm_kit
    export PROTHINT_PATH=/opt/ETP/bin/gmes/ProtHint/bin
    export BEDTOOLS_PATH=/opt/ETP/tools
    export BAMTOOLS_PATH=/opt/ETP/tools
    export TSEBRA_PATH=/opt/TSEBRA/bin
    export CDBTOOLS_PATH=/opt/cdbfasta

    # Define AGAT Perl library path variable (Do not export globally yet)
    AGAT_PERL5LIB=/opt/conda/envs/agat_busco/lib/perl5/5.32/site_perl:/opt/conda/envs/agat_busco/lib/perl5/5.32/vendor_perl:/opt/conda/envs/agat_busco/lib/perl5/5.32/core_perl:/opt/conda/envs/agat_busco/lib/perl5/vendor_perl:/opt/conda/envs/agat_busco/lib/perl5/site_perl:/opt/conda/envs/agat_busco/lib/perl5\${PERL5LIB:+:\$PERL5LIB}

    
    # ---------------------------------------------------------------
    # Pre-mark compleasm placement files directly in the task directory
    # ---------------------------------------------------------------
    mkdir -p mb_downloads
    ln -sf /opt/compleasm_kit/mb_downloads/placement_files mb_downloads/placement_files
    ln -sf /opt/compleasm_kit/mb_downloads/placement_files.done mb_downloads/placement_files.done
    rm -f mb_downloads/placement_files.tmp 2>/dev/null || true

    # Link compleasm.py locally to the task directory without running chmod
    ln -sf /opt/compleasm_kit/compleasm.py ./compleasm.py
    
    # Prepend the absolute path of the local task directory to PATH
    export PATH="\$PWD:\$PATH"

    # ---------------------------------------------------------------
    # Log inputs
    # ---------------------------------------------------------------
    echo '[INFO] BRAKER inputs:'                          > ${prefix}_braker_error.log
    echo '  meta.id=${meta.id}'    >> ${prefix}_braker_error.log
    echo '  meta.use_dual_bams=${meta.use_dual_bams}'    >> ${prefix}_braker_error.log
    echo '  meta.stranded=${meta.stranded}'              >> ${prefix}_braker_error.log
    echo '  bam=${bam}'                                  >> ${prefix}_braker_error.log
    echo '  plus_bam=${plus_bam}'                        >> ${prefix}_braker_error.log
    echo '  minus_bam=${minus_bam}'                      >> ${prefix}_braker_error.log
    echo '  protein_evidence=${protein_evidence}'        >> ${prefix}_braker_error.log
    echo '  gc_probability=${gc_probability}'            >> ${prefix}_braker_error.log
    echo '  bam_flag=${bam_flag}'                        >> ${prefix}_braker_error.log
    echo '  stranded_flag=${stranded_flag}'              >> ${prefix}_braker_error.log
    echo '  gc_flag=${gc_flag}'                          >> ${prefix}_braker_error.log

    # ---------------------------------------------------------------
    # Clean up any pre-existing AUGUSTUS species dir in the real home
    # ---------------------------------------------------------------
    species_dir="\$HOME/.augustus/species/${species}"
    if [ -d "\${species_dir}" ]; then
        echo "[INFO] Removing existing species directory: \${species_dir}" >> ${prefix}_braker_error.log
        rm -rf "\${species_dir}"
    fi
    echo "[INFO] Running BRAKER with species: ${species}" >> ${prefix}_braker_error.log

    # ---------------------------------------------------------------
    # BRAKER RUN (Forced isolated baseline environment execution)
    # ---------------------------------------------------------------
    /usr/bin/perl -I/usr/local/lib/perl5 -I/usr/local/share/perl5 -I/usr/lib/x86_64-linux-gnu/perl5/5.32 -I/usr/share/perl5/5.32 -I/usr/lib/x86_64-linux-gnu/perl-base -I/usr/lib/x86_64-linux-gnu/perl/5.32 -I/usr/share/perl/5.32 /opt/BRAKER/scripts/braker.pl \\
      --genome=${genome_masked} \\
      ${bam_flag} \\
      ${stranded_flag} \\
      --prot_seq=${protein_evidence} \\
      --species="${species}" \\
      ${gc_flag} \\
      --threads ${task.cpus} \\
      --busco_lineage=${busco_db} 2>> ${prefix}_braker_error.log

    mv braker/braker.gtf .

    # ---------------------------------------------------------------
    # BRAKER POST (Switch Context to Safe AGAT/Conda Env Only Now)
    # ---------------------------------------------------------------
    export PATH=/opt/conda/envs/agat_busco/bin:\$PATH
    export PERL5LIB="\${AGAT_PERL5LIB}"

    if [[ -f "braker/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff" ]]; then
        echo ">>> Adding UTRs..."
        python3.8 ${script_dir}/stringtie2utr.py \\
            -g braker.gtf \\
            -s braker/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff \\
            -o braker_with_utrs.gtf
        cat braker_with_utrs.gtf | gtf2gff.pl --gff3 -o braker.gff3
    else
        echo ">>> no RNASeq..."
        cat braker.gtf | gtf2gff.pl --gff3 -o braker.gff3
    fi

    echo ">>> Merging tRNAs..."
    agat_sp_merge_annotations.pl --gff braker.gff3 --gff ${trna_gff} --out merged.gff

    echo ">>> Post-filtering..."
    agat_sp_filter_by_ORF_size.pl -g merged.gff -s 50 -o ${prefix}_filtered.gff
    agat_sp_fix_overlaping_genes.pl -f ${prefix}_filtered_sup50.gff -o ${prefix}_braker.gff3

    echo ">>> Extracting proteins..."
    gffread ${prefix}_braker.gff3 -g ${genome_unmasked} -y ${prefix}.braker.prot.fasta

    echo ">>> Running BUSCO..."
    /opt/conda/envs/busco6_env/bin/python3 /opt/conda/envs/busco6_env/bin/busco \\
        -i ${prefix}.braker.prot.fasta \\
        -o ${prefix}_busco \\
        -m proteins \\
        -l ${busco_db} \\
        -c ${task.cpus}

    mv "${prefix}_busco/short_summary.specific.${busco_db}.${prefix}_busco.txt" ${prefix}_busco_braker.txt
    """
}
