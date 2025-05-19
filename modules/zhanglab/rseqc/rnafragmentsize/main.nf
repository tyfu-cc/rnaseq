process RSEQC_RNAFRAGMENTSIZE {
    tag "${meta.id}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:5.0.3--py39hf95cd2a_0' :
        'biocontainers/rseqc:5.0.3--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  bed

    output:
    tuple val(meta), path("*.fragSize.txt"), emit: txt
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    RNA_fragment_size.py \\
        -i ${bam} \\
        -r ${bed} \\
        > ${prefix}.fragSize.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(RNA_fragment_size.py --version | sed -e "s/RNA_fragment_size.py //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fragSize.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(RNA_fragment_size.py --version | sed -e "s/RNA_fragment_size.py //g")
    END_VERSIONS
    """
}
