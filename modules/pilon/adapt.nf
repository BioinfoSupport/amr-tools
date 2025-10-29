process PILON_ADAPT {
  input:
      tuple val(meta), path('pilon')
  output:
      tuple val(meta), path('pilon.fasta'), emit:fasta
	script:
	"""
			ln -s pilon/pilon.fasta ./
	"""
}
