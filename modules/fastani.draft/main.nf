
process FASTANI {
    container 'registry.gitlab.unige.ch/amr-genomics/orgfinder:v0.2'
    memory '20 GB'
    cpus 8
    time '15 min'
    input:
        tuple val(meta), path(query_fa), path(ref_fa)
    output:
    		tuple val(meta), path('fastani.tsv')
    script:
    		file("refList.txt").text = ref_fa.collect({it.toString()}).join('\n') + '\n'
    		file("queryList.txt").text = query_fa.collect({it.toString()}).join('\n') + '\n'
		    """
		    /app/fastANI --threads ${task.cpus} --refList 'refList.txt' --queryList 'queryList.txt' --output 'fastani.tsv'
		    """
		stub:
				"""
				touch fastani.tsv
				"""
}

