
process CHECKM2_DB {
    container 'quay.io/biocontainers/checkm2:1.1.0--pyh7e72e81_1'
    memory '4 GB'
    cpus 1
    time '1h'
    output:
		    path('checkm2_db', type: 'dir')
    script:
    """
    checkm2 database --download --path checkm2_db
    """
}

