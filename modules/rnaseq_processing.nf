process RNASEQ_PROCESSING {
    tag "$meta.id"
    label 'process_medium'
    container 'docker://abdoallahsharaf/geneforge-rnaseq:2.1'

    publishDir "${params.outdir}/RNASeq", mode: 'copy', pattern: "${meta.id}_RNASeqAll.Stringtie.gtf"
    publishDir "${params.outdir}/RNASeq", mode: 'copy', pattern: "${meta.id}_RNASeqAll.STAR.bam"
    publishDir "${params.outdir}/RNASeq", mode: 'copy', pattern: "${meta.id}_RNASeqAll.transcripts.fasta"
    publishDir "${params.outdir}/RNASeq", mode: 'copy', pattern: "${meta.id}_plus_strand.bam"
    publishDir "${params.outdir}/RNASeq", mode: 'copy', pattern: "${meta.id}_minus_strand.bam"

    input:
      val(meta)
      path(genome)
      val(reads_dir)
      path(script_dir)
      val(stranded)

    output:
    tuple val(meta), path("${meta.id}_RNASeqAll.Stringtie.gtf"), emit: gtf
    tuple val(meta), path("${meta.id}_RNASeqAll.STAR.bam"), emit: bam
    tuple val(meta), path("${meta.id}_RNASeqAll.transcripts.fasta"), emit: transcripts
    tuple val(meta), path("${meta.id}_gtf_summary.txt"), emit: gtf_summary, optional: true
    tuple val(meta), path("${meta.id}_plus_strand.bam"), emit: plus_strand, optional: true
    tuple val(meta), path("${meta.id}_minus_strand.bam"), emit: minus_strand, optional: true
    path "${meta.id}_error.log", emit: error
    tuple val(meta), path("trimmed/${meta.id}_trimmed_1P.fastq.gz"), emit: trimmed_r1
    tuple val(meta), path("trimmed/${meta.id}_trimmed_2P.fastq.gz"), emit: trimmed_r2

    script:
    def prefix = "${meta.id}"
  
    def star_strand = stranded == 'forward' ? '--twopassMode Basic' : 
                      stranded == 'reverse' ? '--twopassMode Basic' : 
                      '--outFilterIntronMotifs RemoveNoncanonical'
    def stringtie_strand = stranded == 'forward' ? '--fr' : stranded == 'reverse' ? '--rf' : ''
    def read_cmd = stranded in ['forward', 'reverse'] ? 'zcat' : 'gunzip -c'
   
    """
    #!/bin/bash
    set -euo pipefail

    # Initialize empty required output structures to avoid Nextflow tracking crashes on premature errors
    mkdir -p trimmed index
    touch ${prefix}_RNASeqAll.Stringtie.gtf ${prefix}_RNASeqAll.STAR.bam ${prefix}_RNASeqAll.transcripts.fasta
    touch trimmed/${prefix}_trimmed_1P.fastq.gz trimmed/${prefix}_trimmed_2P.fastq.gz

    # ---------------------------------------------------------------
    # Restore container PATH environments
    # ---------------------------------------------------------------
    export PATH=/opt/conda/bin:/opt/conda/envs/rnaseq/bin:/usr/local/bin:/usr/bin:/bin:\$PATH

    # Log start
    echo "Starting RNA-Seq processing at \$(date)" > ${prefix}_error.log

    # Validate genome
    if ! grep -q '^>' ${genome} || ! grep -q '[ACGTN]' ${genome}; then
        echo "Invalid FASTA format: missing headers or sequences" >> ${prefix}_error.log
        exit 1
    fi

    # Set directory permissions 
    chmod -R u+rwX trimmed index || { echo "Failed to set directory permissions" >> ${prefix}_error.log; exit 1; }

    # Validate and process FASTQ reads
    r1_files=(\$(find "${reads_dir}" -maxdepth 1 -type f \\( -iname "*_R1*.fastq.gz" -o -iname "*_R1*.fq.gz" -o -iname "*_1*.fastq.gz" -o -iname "*_1*.fq.gz" \\)))
    r2_files=(\$(find "${reads_dir}" -maxdepth 1 -type f \\( -iname "*_R2*.fastq.gz" -o -iname "*_R2*.fq.gz" -o -iname "*_2*.fastq.gz" -o -iname "*_2*.fq.gz" \\)))
    
    if [ \${#r1_files[@]} -eq 0 ] || [ \${#r2_files[@]} -eq 0 ]; then
          echo "No paired-end FASTQ files found in ${reads_dir}" >> ${prefix}_error.log
          exit 1
    fi

    if [ \${#r1_files[@]} -eq 1 ]; then
         # Single dataset: copy files
         cp \${r1_files[0]} ${prefix}_merged_R1.fastq.gz
         cp \${r2_files[0]} ${prefix}_merged_R2.fastq.gz
    else
        # Multiple datasets: merge
        cat \${r1_files[@]} > ${prefix}_merged_R1.fastq.gz
        cat \${r2_files[@]} > ${prefix}_merged_R2.fastq.gz
    fi

    # Validate merged FASTQ files
    if [ ! -s ${prefix}_merged_R1.fastq.gz ] || [ ! -s ${prefix}_merged_R2.fastq.gz ]; then
        echo "Merged FASTQ files are empty" >> ${prefix}_error.log
        exit 1
    fi

    # Determine Trimmomatic call style (Direct wrapper vs Jar execution)
    if command -v trimmomatic &> /dev/null; then
        TRIM_CMD="trimmomatic"
    elif [ -f "/opt/conda/share/trimmomatic/trimmomatic.jar" ]; then
        TRIM_CMD="java -jar /opt/conda/share/trimmomatic/trimmomatic.jar"
    else
        echo "Could not find trimmomatic executable or jar package inside image." >> ${prefix}_error.log
        exit 1
    fi

    # Trim reads with Trimmomatic
    \$TRIM_CMD PE -threads ${task.cpus} -phred33 \\
        ${prefix}_merged_R1.fastq.gz \\
        ${prefix}_merged_R2.fastq.gz \\
        trimmed/${prefix}_trimmed_1P.fastq.gz trimmed/${prefix}_trimmed_1U.fastq.gz \\
        trimmed/${prefix}_trimmed_2P.fastq.gz trimmed/${prefix}_trimmed_2U.fastq.gz \\
        ILLUMINACLIP:${script_dir}/adapt.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 2>> ${prefix}_error.log || \\
        { echo "Trimmomatic failed" >> ${prefix}_error.log; exit 1; }

    # Index genome with STAR
    STAR --runThreadN ${task.cpus} \\
        --runMode genomeGenerate \\
        --genomeDir index \\
        --genomeFastaFiles ${genome} \\
        --genomeSAindexNbases 10 \\
        --limitGenomeGenerateRAM ${task.memory.toBytes()} 2>> ${prefix}_error.log

    # Map reads with STAR
    STAR --runThreadN ${task.cpus} \\
        --genomeDir index \\
        --readFilesIn trimmed/${prefix}_trimmed_1P.fastq.gz trimmed/${prefix}_trimmed_2P.fastq.gz \\
        --readFilesCommand ${read_cmd} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMstrandField intronMotif \\
        ${star_strand} \\
        --outFileNamePrefix ${prefix}_  --limitBAMsortRAM ${task.memory.toBytes()}  2>> ${prefix}_error.log

    # Run StringTie
    stringtie -p ${task.cpus} \\
        -o ${prefix}_RNASeqAll.Stringtie.gtf \\
        ${stringtie_strand} \\
        ${prefix}_Aligned.sortedByCoord.out.bam 2>> ${prefix}_error.log

    # Optional GTF summary
    grep -v "#" ${prefix}_RNASeqAll.Stringtie.gtf | cut -f3 | sort | uniq -c > ${prefix}_gtf_summary.txt

    # Generate transcript FASTA
    gtf_genome_to_cdna_fasta.pl \\
        ${prefix}_RNASeqAll.Stringtie.gtf \\
        ${genome} > ${prefix}_RNASeqAll.transcripts.fasta 2>> ${prefix}_error.log

    # Split BAM by strand if stranded
    if [ "${stranded}" != "no" ]; then
        # Plus strand
        samtools view -@ ${task.cpus} -h ${prefix}_Aligned.sortedByCoord.out.bam | \\
            awk 'BEGIN{OFS="\\t"} /^@/ || (\$2==99 || \$2==83)' | \\
            samtools view -@${task.cpus} -b -o ${prefix}_plus_strand.bam 2>> ${prefix}_error.log

        # Minus strand
        samtools view -@ ${task.cpus} -h ${prefix}_Aligned.sortedByCoord.out.bam | \\
            awk 'BEGIN{OFS="\\t"} /^@/ || (\$2==147 || \$2==163)' | \\
            samtools view -@ ${task.cpus} -b -o ${prefix}_minus_strand.bam 2>> ${prefix}_error.log
    fi

    # Rename BAM for output
    mv ${prefix}_Aligned.sortedByCoord.out.bam ${prefix}_RNASeqAll.STAR.bam
    """
}
