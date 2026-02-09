process MEDAKA_CONSENSUS {
    container 'docker.io/ontresearch/medaka:shac4e11bfa4e65668b28739ba32edc3af12baf7574-amd64'
    memory '10 GB'
    cpus 4
    time '2h'
    input:
        tuple val(meta), path('assembly.fasta'), path('long_reads.fastq.gz'), val(args)
    output:
        tuple val(meta), path('medaka',type:'dir')
    script:
		    """
		    medaka_consensus \\
		      ${args?:''} \\
		      ${task.ext.args?:'--bacteria'} \\
			    -f -t ${task.cpus} \\
			    -d assembly.fasta \\
			    -i long_reads.fastq.gz \\
			    -o medaka
		    """
		stub:		    
		    """
		    mkdir -p medaka/
		    touch medaka/consensus.fasta
		    """
}
