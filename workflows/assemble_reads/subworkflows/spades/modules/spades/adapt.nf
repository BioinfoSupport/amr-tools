process SPADES_ADAPT {
  input:
      tuple val(meta), path('spades')
  output:
      tuple val(meta), path('assembly.fasta'), emit:fasta
	script:
	"""
		if [ -f spades/scaffolds.fasta ]; then
			ln -s spades/scaffolds.fasta assembly.fasta
		else
			ln -s spades/contigs.fasta assembly.fasta
		fi
	"""
} 