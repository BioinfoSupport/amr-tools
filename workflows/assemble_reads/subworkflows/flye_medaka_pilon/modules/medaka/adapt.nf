process MEDAKA_ADAPT {
  input:
      tuple val(meta), path('medaka')
  output:
      tuple val(meta), path('consensus.fasta'), emit:fasta
	script:
	"""
		ln -s medaka/consensus.fasta .
	"""
} 