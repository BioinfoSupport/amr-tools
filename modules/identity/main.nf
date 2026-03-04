process IDENTITY {
    memory '2 GB'
    cpus 1
    time '1h'
    executor 'local'
    input:
    	tuple val(meta), path(file)
    output:
    	tuple val(meta), path(file)
    script:
    """
    """
}



