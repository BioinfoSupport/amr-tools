
process PROKKA {
	  container "docker.io/staphb/prokka:1.14.6"
    memory '12 GB'
    cpus 4
    time '1h'
    input:
        tuple val(meta), path('assembly.fasta')
    output:
				tuple val(meta), path("prokka/", type: 'dir')
    script:
		    """
				prokka ${task.ext.args?:''} --outdir 'prokka/' --cpus ${task.cpus} --prefix prokka assembly.fasta
		    """    
}

