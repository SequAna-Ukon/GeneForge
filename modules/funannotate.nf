process FUNANNOTATE {
    tag "${meta.id}"
    label 'process_high'

    container 'docker://abdoallahsharaf/geneforge-funannotate:2.1'

    publishDir "${params.outdir}/funannotate", mode: 'copy', pattern: '*_funannotate.gff3'
    publishDir "${params.outdir}/funannotate", mode: 'copy', pattern: '*.funannotate.prot.fasta'
    publishDir "${params.outdir}/funannotate", mode: 'copy', pattern: '*_busco_funannotate.txt'

    input:
    tuple val(meta),
          path(genome_masked),
          path(genome_unmasked),
          val(highconf_tbl),
          path(protein_evidence),
          path(gtf),
          path(bam),
          path(transcripts),
          path(rnaseq_r1),
          path(rnaseq_r2),
          path(genemark_key),
          path(genemark_tar),
          path(nanopore_mrna),
          path(pacbio_isoseq)

    output:
    tuple val(meta), path("${meta.id}_funannotate.gff3"),       emit: gff3
    tuple val(meta), path("${meta.id}.funannotate.prot.fasta"), emit: proteins
    tuple val(meta), path("${meta.id}_busco_funannotate.txt"),  emit: busco_summary
    path "${meta.id}_error.log", optional: true,                 emit: error_log

    script:
    def prefix        = meta.id
    def species       = meta.species
    def species_db    = meta.species.replaceAll(' ', '_')
    def organism      = meta.organism
    def busco_db      = meta.busco_db
    def busco_db_fun  = meta.busco_db_fun
    def funanno_db    = meta.funanno_DB ?: "${launchDir}/work/funannotate_db"
    def trnascan_flag = (highconf_tbl && file(highconf_tbl).exists()) ? "--trnascan ${highconf_tbl}" : ''

    def isReal = { f -> f && f.toString() != "" && !f.toString().contains('NO_') && !f.toString().endsWith('.empty') }
    def hasShortReads = isReal(rnaseq_r1) && isReal(rnaseq_r2)
    def hasBAM        = isReal(bam)
    def hasNano       = isReal(nanopore_mrna)
    def hasPacBio     = isReal(pacbio_isoseq)
    def genemark_mode = hasBAM ? 'ET' : 'ES'
    def stranded_flag = (meta.stranded && meta.stranded != 'no')
        ? "--stranded ${meta.stranded == 'forward' ? 'FR' : meta.stranded == 'reverse' ? 'RF' : meta.stranded}"
        : ''

    def runProteinEvidence = isReal(protein_evidence) ? "true" : "false"
    def runBam             = hasBAM              ? "true" : "false"
    def runGtf             = isReal(gtf)         ? "true" : "false"
    def runTranscripts     = isReal(transcripts) ? "true" : "false"

    """
    #!/bin/bash
    set -euo pipefail

    export USER=\$(whoami)
    
    # -------------------------------------------------------------------------
    # 0. Restore PATH — funannotate env first, busco_env for busco only
    # -------------------------------------------------------------------------
    export PATH=/opt/conda/envs/funannotate/bin:/opt/conda/envs/busco_env/bin:/usr/local/bin:/usr/bin:/bin:\$PATH
    export CONDA_PREFIX=/opt/conda/envs/funannotate
    export EVM_HOME=/opt/conda/envs/funannotate/opt/evidencemodeler-1.1.1
    export AUGUSTUS_CONFIG_PATH=/opt/conda/envs/funannotate/config
    export PERL5LIB=/opt/conda/envs/funannotate/lib/perl5/5.32/site_perl:/opt/conda/envs/funannotate/lib/perl5/site_perl:/opt/conda/envs/funannotate/lib/perl5/5.32/vendor_perl:/opt/conda/envs/funannotate/lib/perl5/vendor_perl:/opt/conda/envs/funannotate/lib/perl5/5.32/core_perl:/opt/conda/envs/funannotate/lib/perl5/core_perl
    export TRINITYHOME=/opt/conda/envs/funannotate/opt/trinity-2.8.5
    export TRINITY_HOME=/opt/conda/envs/funannotate/opt/trinity-2.8.5
    export PATH="\$TRINITYHOME:\$PATH"
    export QUARRY_PATH=/opt/conda/envs/funannotate/bin


    mkdir -p funannotate_${prefix} gmes

    # -------------------------------------------------------------------------
    # 1. Environment
    # -------------------------------------------------------------------------
    ENV_PREFIX=\$(dirname \$(dirname \$(which funannotate)))
    export PASAHOME="\$ENV_PREFIX/opt/pasa-2.5.3"
    export PERL5LIB="\$PASAHOME/SAMPLE_HOOKS:\$PASAHOME/PerlLib:\${PERL5LIB:+:\$PERL5LIB}"
    export PATH="\$PASAHOME/bin:\$PATH"
    echo "Using ENV_PREFIX: \$ENV_PREFIX" >> ${prefix}_error.log

    # -------------------------------------------------------------------------
    # 3. GeneMark
    # -------------------------------------------------------------------------
    if [ -f "${genemark_key}" ] && [[ "${genemark_key}" != *.empty ]]; then
        export HOME=\$(pwd)
        gunzip -c "${genemark_key}" > "\$HOME/.gm_key"
    fi
    if [ -f "${genemark_tar}" ] && [[ "${genemark_tar}" != *.empty ]]; then
        tar -xzf ${genemark_tar} -C gmes
        subdir=\$(find gmes -mindepth 1 -maxdepth 1 -type d | head -n1)
        if [ -n "\$subdir" ]; then
            mv "\$subdir"/* gmes/ && rmdir "\$subdir"
        fi
        chmod +x gmes/gmes_petap.pl
        (cd gmes && perl change_path_in_perl_scripts.pl \$(which perl))
    fi
    export GENEMARK_PATH=\$(realpath gmes)
    export PATH=\$GENEMARK_PATH:\$PATH

    # -------------------------------------------------------------------------
    # 4. Funannotate database
    # -------------------------------------------------------------------------
    export FUNANNOTATE_DB="${funanno_db}"
    if [[ ! -d "\$FUNANNOTATE_DB" ]]; then
        if [[ "\$FUNANNOTATE_DB" == *"/work/funannotate_db" ]]; then
            mkdir -p "\$FUNANNOTATE_DB"
        else
            echo "ERROR: Provided DB path \$FUNANNOTATE_DB does not exist." >> ${prefix}_error.log
            echo "       Please run: funannotate setup --install all -b ${busco_db_fun} --wget -f --database \$FUNANNOTATE_DB" >> ${prefix}_error.log
            exit 1
        fi
    fi

    if [[ -f "\$FUNANNOTATE_DB/funannotate-db-info.txt" ]]; then
        echo ">>> Using existing funannotate database at \$FUNANNOTATE_DB" >> ${prefix}_error.log
    else
        (
            flock -x 200
            if [[ ! -f "\$FUNANNOTATE_DB/funannotate-db-info.txt" ]]; then
                echo ">>> Installing funannotate database to \$FUNANNOTATE_DB" >> ${prefix}_error.log
                funannotate setup --install all -b ${busco_db_fun} --wget -f --database "\$FUNANNOTATE_DB"
            else
                echo ">>> Reusing existing funannotate database at \$FUNANNOTATE_DB" >> ${prefix}_error.log
            fi
        ) 200>"\$FUNANNOTATE_DB/.install.lock"
    fi

    # -------------------------------------------------------------------------
    # 4b. Pre-process long reads: FASTQ->FASTA, RNA->DNA
    #     Filters reads <200bp to prevent seqclean IndexError in funannotate
    # -------------------------------------------------------------------------
    TRAIN_LR_FLAG=""

    preprocess_lr() {
        local infile="\$1"
        local outfile="\$2"
        local tmpfile="\${outfile}.tmp"
        local filtered="\${outfile}.filtered"
        zcat -f "\$infile" > "\$tmpfile" || cp "\$infile" "\$tmpfile"
        local first_char
        first_char=\$(head -c1 "\$tmpfile")
        if [[ "\$first_char" == "@" ]]; then
            local h s skip
            while IFS= read -r h && IFS= read -r s && IFS= read -r skip && IFS= read -r skip; do
                s=\${s//U/T}
                s=\${s//u/t}
                h=\${h:1}
                h=\${h%% *}
                [[ -n "\$s" ]] && printf '>%s\n%s\n' "\$h" "\$s"
            done < "\$tmpfile" > "\$filtered"
        else
            sed '/^>/ s/ .*//' "\$tmpfile" | tr 'Uu' 'Tt' > "\$filtered"
        fi
        rm -f "\$tmpfile"

        # Filter sequences <200bp to prevent seqclean IndexError
        awk '/^>/{header=\$0; next} length(\$0)>=200{print header; print}' \
            "\$filtered" > "\$outfile"
        rm -f "\$filtered"

        local nseqs
        nseqs=\$(grep -c "^>" "\$outfile" 2>/dev/null || echo 0)
        echo "    Preprocessed \$nseqs sequences (>=200bp) -> \$outfile" >> ${prefix}_error.log

        if [[ "\$nseqs" -eq 0 ]]; then
            echo "WARNING: No long reads passed length filter, skipping long read training" >> ${prefix}_error.log
            return 1
        fi
        return 0
    }

    if [[ "${hasNano}" == "true" ]]; then
        echo ">>> Pre-processing ONT reads..." >> ${prefix}_error.log
        if preprocess_lr "${nanopore_mrna}" lr_processed.fa; then
            TRAIN_LR_FLAG="--nanopore_mrna \$(pwd)/lr_processed.fa"
        else
            TRAIN_LR_FLAG=""
        fi
    elif [[ "${hasPacBio}" == "true" ]]; then
        echo ">>> Pre-processing PacBio reads..." >> ${prefix}_error.log
        if preprocess_lr "${pacbio_isoseq}" lr_processed.fa; then
            TRAIN_LR_FLAG="--pacbio_isoseq \$(pwd)/lr_processed.fa"
        else
            TRAIN_LR_FLAG=""
        fi
    fi

    # -------------------------------------------------------------------------
    # 5. Training
    # -------------------------------------------------------------------------
    TRAIN_SR_FLAG=""
    [[ "${hasShortReads}" == "true" ]] && \
        TRAIN_SR_FLAG="-l ${rnaseq_r1} -r ${rnaseq_r2} --no_trimmomatic ${stranded_flag}"

    if [[ -n "\$TRAIN_SR_FLAG" || -n "\$TRAIN_LR_FLAG" ]]; then
        echo ">>> Running training..." >> ${prefix}_error.log
        funannotate train --species "${species}" -i ${genome_unmasked} -o funannotate_${prefix} \
            \$TRAIN_SR_FLAG \$TRAIN_LR_FLAG \
            --max_intronlen 50000 \
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G 2>> ${prefix}_error.log
    else
        echo ">>> No RNA evidence, skipping train." >> ${prefix}_error.log
    fi

    # Materialise symlinks left by train
    for f in funannotate_${prefix}/training/funannotate_train.coordSorted.bam \
             funannotate_${prefix}/training/funannotate_train.transcripts.gff3 \
             funannotate_${prefix}/training/funannotate_long-reads.fasta \
             funannotate_${prefix}/training/funannotate_train.trinity-GG.fasta; do
        if [[ -L "\$f" && -e "\$f" ]]; then
            cp --remove-destination \$(realpath "\$f") "\$f"
        fi
    done
    if [[ -f funannotate_${prefix}/training/funannotate_train.coordSorted.bam ]]; then
        samtools index funannotate_${prefix}/training/funannotate_train.coordSorted.bam
    fi

    # -------------------------------------------------------------------------
    # 5b. Filter protein evidence using PASA TransDecoder proteins
    # -------------------------------------------------------------------------
    PROTEIN_EVIDENCE=""

    if [[ "${runProteinEvidence}" == "true" ]]; then
        if [[ -d "funannotate_${prefix}/training/pasa" ]]; then
            PASA_PEP=\$(find funannotate_${prefix}/training/pasa -name '*.assemblies.fasta.transdecoder.pep' | head -n1)
        else
            PASA_PEP=""
        fi

        if [[ -n "\$PASA_PEP" ]]; then
            echo ">>> Filtering protein evidence using PASA TransDecoder proteins..." >> ${prefix}_error.log
            echo "    Original DB size: \$(grep -c '>' ${protein_evidence}) proteins" >> ${prefix}_error.log
            diamond blastp \\
                -p ${task.cpus} \\
                -q "\$PASA_PEP" \\
                -d ${protein_evidence} \\
                --max-target-seqs 5 \\
                --outfmt 6 qseqid sseqid \\
                -o diamond_protein_hits.txt 2>> ${prefix}_error.log
            cut -f2 diamond_protein_hits.txt | sort -u > protein_hit_ids.txt
            seqtk subseq ${protein_evidence} protein_hit_ids.txt > filtered_proteins.fa
            echo "    Filtered DB size: \$(grep -c '>' filtered_proteins.fa) proteins" >> ${prefix}_error.log
            PROTEIN_EVIDENCE=\$(realpath filtered_proteins.fa)
        else
            echo ">>> WARNING: PASA training directory or .pep not found, using original protein evidence" >> ${prefix}_error.log
            PROTEIN_EVIDENCE=\$(realpath ${protein_evidence})
        fi
    fi

    # -------------------------------------------------------------------------
    # 6. Predict
    # -------------------------------------------------------------------------
    RNA_FLAGS=""
    [[ "${runBam}" == "true" ]]         && RNA_FLAGS+=" --rna_bam ${bam}"
    [[ "${runGtf}" == "true" ]]         && RNA_FLAGS+=" --stringtie ${gtf}"
    [[ "${runTranscripts}" == "true" ]] && RNA_FLAGS+=" --transcript_evidence ${transcripts}"

    run_predict() {
        local gm_mode="\$1"
        local prot_flag=""
        [[ -n "\$PROTEIN_EVIDENCE" ]] && prot_flag="--protein_evidence \$PROTEIN_EVIDENCE"

        funannotate predict \\
            --species "${species}" \\
            -i ${genome_masked} \\
            -o funannotate_${prefix} \\
            --name "${prefix}" \\
            \$RNA_FLAGS \\
            \$prot_flag \\
            ${trnascan_flag} \\
            --organism "${organism}" \\
            --database \$FUNANNOTATE_DB \\
            --busco_db ${busco_db_fun} \\
            --genemark_mode "\$gm_mode" \\
            --GENEMARK_PATH \$GENEMARK_PATH \\
            --cpus ${task.cpus} 2>> ${prefix}_error.log
    }

    if ! run_predict "${genemark_mode}"; then
        echo "WARNING: GeneMark-${genemark_mode} failed; retrying with ES" >> ${prefix}_error.log
        rm -rf funannotate_${prefix}/predict_misc/genemark* || true
        run_predict "ES"
    fi

    # -------------------------------------------------------------------------
    # 7. Update
    # -------------------------------------------------------------------------
    funannotate update -i funannotate_${prefix} --species "${species}" \\
     --cpus ${task.cpus} 2>> ${prefix}_error.log

    # -------------------------------------------------------------------------
    # 8. Post-process
    # -------------------------------------------------------------------------
    GFF_FILE=\$(find funannotate_${prefix}/update_results/ -name '*.gff3' | head -n1)
    if [[ -f "\$GFF_FILE" ]]; then
        agat_sp_filter_by_ORF_size.pl -g \$GFF_FILE -s 50 -o ${prefix}_filtered.gff
        agat_sp_fix_overlaping_genes.pl -f ${prefix}_filtered_sup50.gff -o ${prefix}_funannotate.gff3
        gffread ${prefix}_funannotate.gff3 -g ${genome_unmasked} -y ${prefix}.funannotate.prot.fasta
        /opt/conda/envs/busco_env/bin/python /opt/conda/envs/busco_env/bin/busco \\
            -i ${prefix}.funannotate.prot.fasta \\
            -o ${prefix}_busco \\
            -m proteins \\
            -l ${busco_db} \\
            -c ${task.cpus}
        cp "${prefix}_busco/short_summary.specific.${busco_db}.${prefix}_busco.txt" \\
           ${prefix}_busco_funannotate.txt || true
    else
        echo "ERROR: No GFF3 found after update." >> ${prefix}_error.log
        exit 1
    fi
    """
}

