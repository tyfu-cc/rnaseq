process CHECKSTRANDEDNESS {
    tag "${meta.id}"
    label "process_medium"

    container "ghcr.io/tyfu-cc/checkstrandedness:1.0.1"

    input:
    tuple val(meta), path(reads)
    path  gtf
    path  transcripts

    output:
    tuple val(meta), path("*.strandedness.txt"), emit: txt
    path "versions.yml"                        , emit: versions, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_reads = meta.single_end 
        ? "--reads_1 ${reads[0]}" 
        : "--reads_1 ${reads[0]} --reads_2 ${reads[1]}"
    """
    check_strandedness \\
        --gtf ${gtf} \\
        --transcripts ${transcripts} \\
        ${input_reads} \\
        > ${prefix}.check_strandedness.output.txt

    last_line=`tail -n 1 ${prefix}.check_strandedness.output.txt`

    if [[ \$last_line == *"FR/fr-stranded"* || \$last_line == *"FR/fr-secondstrand"* ]]; then
        echo "forward" > ${prefix}.strandedness.txt
    elif [[ \$last_line == *"RF/rf-stranded" || \$last_line == *"RF/fr-firststrand"* ]]; then
        echo "reverse" > ${prefix}.strandedness.txt
    else
        echo "unstranded" > ${prefix}.strandedness.txt
    fi

    # cat <<-END_VERSIONS > versions.yml
    # "${task.process}":
    #     cutadapt: \$(cutadapt --version)
    # END_VERSIONS
    """

    //stub:
    //def prefix = task.ext.prefix ?: "${meta.id}"
    //def trimmed = meta.single_end ? "${prefix}.trimmed.fastq.gz" : "${prefix}_1.trimmed.fastq.gz ${prefix}_2.trimmed.fastq.gz"
    //"""
    //touch ${prefix}.cutadapt.log
    //touch ${trimmed}

    //cat <<-END_VERSIONS > versions.yml
    //"${task.process}":
    //    cutadapt: \$(cutadapt --version)
    //END_VERSIONS
    //"""
}
