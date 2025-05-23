process RSEM_CALCULATEEXPRESSION {
    tag "${meta.id}"
    label "process_high"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:64aad4a4e144878400649e71f42105311be7ed87-0' :
        'biocontainers/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:64aad4a4e144878400649e71f42105311be7ed87-0' }"

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("*.genes.results")         , emit: counts_gene
    tuple val(meta), path("*.isoforms.results")      , emit: counts_transcript
    tuple val(meta), path("*.stat")                  , emit: stat
    tuple val(meta), path("*.log")                   , emit: logs
    path  "versions.yml"                             , emit: versions

    tuple val(meta), path("*.STAR.genome.bam")       , emit: bam_star      , optional: true
    tuple val(meta), path("${prefix}.genome.bam")    , emit: bam_genome    , optional: true
    tuple val(meta), path("${prefix}.transcript.bam"), emit: bam_transcript, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    def strandedness = ""
    if (meta.strandedness == "forward") {
        strandedness = "--strandedness forward"
    } else if (meta.strandedness == "reverse") {
        strandedness = "--strandedness reverse"
    }
    def paired_end = meta.single_end ? "" : "--paired-end"
    """
    INDEX=`find -L ./ -name "*.grp" | sed 's/\\.grp\$//'`
    rsem-calculate-expression \\
        --num-threads ${task.cpus} \\
        --temporary-folder ./tmp/ \\
        ${strandedness} \\
        ${paired_end} \\
        ${args} \\
        ${reads} \\
        \$INDEX \\
        ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rsem: \$(rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g")
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.genes.results
    touch ${prefix}.isoforms.results
    touch ${prefix}.stat
    touch ${prefix}.log
    touch ${prefix}.STAR.genome.bam
    touch ${prefix}.genome.bam
    touch ${prefix}.transcript.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rsem: \$(rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g")
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
