process FLYE_ADAPT {
  input:
      tuple val(meta), path('flye')
  output:
      tuple val(meta), path('assembly.fasta'), emit:fasta
	script:
	"""
			ln -s flye/assembly.fasta .
	"""
} 