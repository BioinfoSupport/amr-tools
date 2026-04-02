
process FASTANI {
    container 'docker.io/staphb/fastani:1.34'
    memory '20 GB'
    cpus 8
    time '15 min'
    input:
        tuple val(meta), path('query/*'), path('ref/*')
    output:
    		tuple val(meta), path('fastani.tsv')
    script:
		    """
		    find ref/ -type f -or -type l > refList.txt
		    find query/ -type f -or -type l > queryList.txt
		    fastANI --threads ${task.cpus} ${task.ext.args?:''} --refList refList.txt --queryList queryList.txt --output fastani.tsv
		    """
		stub:
				"""
				touch fastani.tsv
				"""
}

