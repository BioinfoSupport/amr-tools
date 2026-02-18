
process FQ_SUBSAMPLE {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory '2 GB'
    cpus 4
    time '30 min'
    input:
    		tuple val(meta), path(reads)
    output:
    		tuple val(meta), path("subsampled_reads_*.fastq.gz")
    script:
    		reads = reads instanceof List?reads:[reads]
    		def bp_limit = task.ext.bp_limit?:-1
    		if (bp_limit<0) {
    			reads.withIndex().collect({x,i -> "ln -s '${x}' 'subsampled_reads_${i+1}.fastq.gz'"}).join("\n")
    		} else {
					"""
					set +o pipefail
					awk '
						BEGIN{bp=0;for(i=1;i<length(ARGV);i++){print ARGV[i]}}
						{for(i=1;i<length(ARGV);i++){"gzip -dc " ARGV[i] | getline line[i "_" NR%4]}}
						NR%4==0 {
					    for(i=1;i<length(ARGV);i++) {
					      if (line[i "_" 1]!~/^@/) {print "Expected @ at line " NR " of file " ARGV[i] > "/dev/stderr";exit 1}
					      if (line[i "_" 3]!~/^[+]/) {print "Expected + at line " NR " of file " ARGV[i] > "/dev/stderr";exit 1}
					      bp += length(line[i "_" 2])
					      ocmd = "gzip > subsampled_reads_" i ".fastq.gz"
					      print line[i "_" 1] | ocmd
					      print line[i "_" 2] | ocmd
					      print line[i "_" 3] | ocmd
					      print line[i "_" 0] | ocmd
					    }
					    if (bp>=${bp_limit}){exit 0}
					  }
					' ${reads.join(" ")}
					"""
    		}
    stub:
    	reads.withIndex().collect({x,i -> "touch 'subsampled_reads_${i+1}.fastq.gz'"}).join("\n")
}

