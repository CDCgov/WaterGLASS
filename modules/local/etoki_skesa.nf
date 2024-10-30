process ETOKI_MLST_SKESA {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(scaffolds)
    path(reference_alleles)
    path(etoki_md5sum)

    output:
    tuple val(meta), path('*.etoki.fasta')    , emit: etoki_alleles_fasta
    path 'versions.yml'                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = scaffolds.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? scaffolds.getBaseName() : scaffolds

    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${scaffolds} > ${fasta_name}
    fi
    singularity exec ${projectDir}/bin/etoki_latest.sif EToKi.py MLSType \\
        -i ${fasta_name} \\
        -r $reference_alleles \\
        -k ${prefix} \\
        -d $etoki_md5sum \\
        -o ${prefix}.etoki.fasta \\
        $args \\
    &> etoki_mlst_main.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        etoki_mlst_main: \$(etoki_mlst_main.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}

