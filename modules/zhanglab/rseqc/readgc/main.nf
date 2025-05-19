process RSEQC_READGC {
    tag "${meta.id}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:5.0.3--py39hf95cd2a_0' :
        'biocontainers/rseqc:5.0.3--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.GC.xls")     , emit: xls
    tuple val(meta), path("*.GC_plot.r")  , emit: rscript
    tuple val(meta), path("*.GC_plot.pdf"), emit: pdf
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    read_GC.py \\
        -i ${bam} \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(read_GC.py --version | sed -e "s/read_GC.py //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.GC.xls
    touch ${prefix}.GC_plot.r
    touch ${prefix}.GC_plot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(read_GC.py --version | sed -e "s/read_GC.py //g")
    END_VERSIONS
    """
}
