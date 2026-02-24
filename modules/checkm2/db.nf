
process CHECKM2_DB {
    container 'community.wave.seqera.io/library/checkm2:1.1.0--60f287bc25d7a10d'
    memory '4 GB'
    cpus 1
    time '1h'
    output:
		    path('checkm2_db', type: 'dir')
    script:
    """
    export CHECKM2_DB_PATH="$PWD/checkm2_db"
    checkm2 database --download
    """
}

