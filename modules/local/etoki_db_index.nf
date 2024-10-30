process ETOKI_DB_INDEX {
    // tag "$fasta"
    label 'process_medium'

    input:
    path(alleles_reformat)
    path(reference_alleles)
    
    output:
    path "*.csv"            , emit: etoki_md5sum 
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
       
    singularity exec ${projectDir}/bin/etoki_05142024.sif EToKi.py MLSTdb \\
        -i $alleles_reformat \\
        -r $reference_alleles \\
        -d etoki_index_md5hash.csv \\
        $args \\
        &> etoki_index.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        etoki_index: \$(etoki_index_indexversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
