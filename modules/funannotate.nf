nextflow.enable.dsl=2

process FUNANNOTATE {
    tag "${meta.id}"
    label 'process_high'

    conda 'bioconda::funannotate=1.8.17 bioconda::agat=1.4.0 bioconda::gffread=0.12.7 bioconda::snap bioconda::busco=5.4.7 bioconda::transdecoder=5.5.0'
 
    publishDir "${params.outdir}/funannotate", mode: 'copy', pattern: '*_funannotate.gff3'
    publishDir "${params.outdir}/funannotate", mode: 'copy', pattern: '*.funannotate.prot.fasta'
    publishDir "${params.outdir}/funannotate", mode: 'copy', pattern: '*_busco_funannotate.txt'

    input:
    tuple val(meta),
          path(genome_masked),
          path(genome_unmasked),
          val(highconf_tbl),
          path(protein_evidence),
          val(gtf),
          val(bam),
          val(transcripts),
          val(rnaseq_r1),
          val(rnaseq_r2),
          path(genemark_key),
          path(genemark_tar),
          val(nanopore_mrna),
          val(pacbio_isoseq)

    output:
    tuple val(meta), path("${meta.id}_funannotate.gff3"), emit: gff3
    tuple val(meta), path("${meta.id}.funannotate.prot.fasta"), emit: proteins
    tuple val(meta), path("${meta.id}_busco_funannotate.txt"), emit: busco_summary
    path "${meta.id}_error.log", optional: true, emit: error_log

    script:
    def prefix = meta.id
    def species = meta.species
    def organism = meta.organism
    def busco_db = meta.busco_db
    def busco_db_fun = meta.busco_db_fun
    def funanno_db = meta.funanno_DB
    def use_setup = (!funanno_db || !file(funanno_db).exists())
    def trnascan_flag = (highconf_tbl && file(highconf_tbl).exists()) ? "--trnascan ${highconf_tbl}" : ''
    def genemark_mode = (rnaseq_r1 && rnaseq_r2 && gtf && bam && transcripts) ? 'ET' : 'ES'
    
    def stranded_flag = meta.stranded && meta.stranded != 'no' ? "--stranded ${meta.stranded == 'forward' ? 'FR' : meta.stranded == 'reverse' ? 'RF' : meta.stranded}" : ''
    
    
    def bash_vars = [
        "prefix='${prefix}'",
        "gtf='${gtf ? gtf.getName() : ''}'",
        "bam='${bam ? bam.getName() : ''}'",
        "transcripts='${transcripts ? transcripts.getName() : ''}'",
        "rnaseq_r1='${rnaseq_r1 ? rnaseq_r1.getName() : ''}'",
        "rnaseq_r2='${rnaseq_r2 ? rnaseq_r2.getName() : ''}'"
    ].join('\n')
    def longread_flag = ''
    if (nanopore_mrna && file(nanopore_mrna).exists()) {
        longread_flag = "--nanopore_mrna ${nanopore_mrna}"
    } else if (pacbio_isoseq && file(pacbio_isoseq).exists()) {
        longread_flag = "--pacbio_isoseq ${pacbio_isoseq}"
    }

    return """
    #!/bin/bash
    set -euo pipefail

    ${bash_vars}

    mkdir -p funannotate_\${prefix} gmes

    if [ -f "${genemark_key}" ]; then
        export HOME=\$(pwd)
        gunzip -c "${genemark_key}" > "\$HOME/.gm_key"
    fi

    if [ -f "${genemark_tar}" ]; then
        tar -xzf ${genemark_tar} -C gmes
        subdir=\$(find gmes -mindepth 1 -maxdepth 1 -type d | head -n1)
        if [ -n "\$subdir" ]; then
            mv "\$subdir"/* gmes/
            rmdir "\$subdir"
        fi
        chmod +x gmes/gmes_petap.pl
        (cd gmes && perl change_path_in_perl_scripts.pl \$(which perl))
    fi

    export GENEMARK_PATH=\$(realpath gmes)
    export PATH=\$GENEMARK_PATH:\$PATH

    if ${use_setup}; then
        mkdir -p ./funannotate_db
        funannotate setup --install all -b ${busco_db_fun} -f --database ./funannotate_db
        export FUNANNOTATE_DB=\$(realpath ./funannotate_db)
    else
        export FUNANNOTATE_DB="${funanno_db}"
    fi

    # Fix permissions for aux_scripts
    SITE_PACKAGES=\$(python3 -c 'import site; print(site.getsitepackages()[0])' 2>/dev/null || echo "")
    echo "SITE_PACKAGES: \$SITE_PACKAGES" >> ${prefix}_error.log
    if [ -n "\$SITE_PACKAGES" ]; then
        AUX_SCRIPTS_DIR="\$SITE_PACKAGES/funannotate/aux_scripts"
        echo "AUX_SCRIPTS_DIR: \$AUX_SCRIPTS_DIR" >> ${prefix}_error.log
        if [ -d "\$AUX_SCRIPTS_DIR" ]; then
            chmod -R +x "\$AUX_SCRIPTS_DIR" 2>> ${prefix}_error.log || {
                echo "WARNING: Failed to set executable permissions for \$AUX_SCRIPTS_DIR" >> ${prefix}_error.log
            }
        else
            echo "WARNING: Directory \$AUX_SCRIPTS_DIR not found" >> ${prefix}_error.log
        fi
    else
        echo "WARNING: Could not determine site-packages directory for Python" >> ${prefix}_error.log
    fi

    # Build RNA-seq flags
    rna_bam_flag=""
    stringtie_flag=""
    transcript_flag=""

    if [[ -f "${bam}" && -s "${bam}" ]]; then
        rna_bam_flag="--rna_bam ${bam}"
    fi
    if [[ -f "${gtf}" && -s "${gtf}" ]]; then
        stringtie_flag="--stringtie ${gtf}"
    fi
    if [[ -f "${transcripts}" && -s "${transcripts}" ]]; then
        transcript_flag="--transcript_evidence ${transcripts}"
    fi

    echo "DEBUG: RNA-seq inputs: r1=${rnaseq_r1} r2=${rnaseq_r2} gtf=${gtf} bam=${bam} transcripts=${transcripts}" >> \${prefix}_error.log

    # Short-read training FIRST (if all short-read inputs are valid)
    if [[ ! -z "${rnaseq_r1}" && -s "${rnaseq_r1}" && \
          ! -z "${rnaseq_r2}" && -s "${rnaseq_r2}" && \
          -f "${gtf}" && -s "${gtf}" && \
          -f "${bam}" && -s "${bam}" && \
          -f "${transcripts}" && -s "${transcripts}" ]]; then
        echo ">>> Running short-read training..." >> \${prefix}_error.log
        funannotate train --species "${species}" -i ${genome_unmasked} -o funannotate_\${prefix} \
            -l ${rnaseq_r1} -r ${rnaseq_r2} --no_trimmomatic ${stranded_flag} \
            --cpus ${task.cpus} 2>> \${prefix}_error.log
    fi

    # THEN long-read training (if longread exists)
    if [[ -n "${longread_flag}" ]]; then
        echo ">>> Running long-read training..." >> \${prefix}_error.log
        funannotate train --species "${species}" -i ${genome_unmasked} -o funannotate_\${prefix} \
            ${longread_flag} --cpus ${task.cpus} 2>> \${prefix}_error.log
    fi

    

    echo ">>> Running predict step..." >> \${prefix}_error.log
    funannotate predict --species "${species}" -i ${genome_masked} -o funannotate_\${prefix} --name "${prefix}" \\
        \${rna_bam_flag} \\
        \${stringtie_flag} \\
        \${transcript_flag} \\
        ${protein_evidence ? "--protein_evidence ${protein_evidence}" : ''} \\
        ${trnascan_flag} \\
        --organism ${organism} \\
        --database \$FUNANNOTATE_DB \\
        --busco_db ${busco_db_fun} \\
        --genemark_mode ${genemark_mode} \\
        --GENEMARK_PATH \$GENEMARK_PATH \\
        --cpus ${task.cpus} 2>> \${prefix}_error.log

    funannotate update -i funannotate_\${prefix} --species "${species}" --cpus ${task.cpus} 2>> \${prefix}_error.log

    find funannotate_\${prefix}/update_results/ -name '*.gff3' | head -n1 > gff_path.txt
    GFF_FILE=\$(cat gff_path.txt)
    agat_sp_filter_by_ORF_size.pl -g \$GFF_FILE -s 50 -o \${prefix}_filtered.gff
    agat_sp_fix_overlaping_genes.pl -f \${prefix}_filtered_sup50.gff -o \${prefix}_funannotate.gff3
    gffread \${prefix}_funannotate.gff3 -g ${genome_unmasked} -y \${prefix}.funannotate.prot.fasta

    busco -i \${prefix}.funannotate.prot.fasta -o \${prefix}_busco -m proteins -l ${busco_db} -c ${task.cpus}
    cp "${prefix}_busco/short_summary.specific.${busco_db}.${prefix}_busco.txt" ${prefix}_busco_funannotate.txt
    """
}
return FUNANNOTATE

