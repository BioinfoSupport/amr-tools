process HYBRACTER_ADAPT {
  input:
      tuple val(meta), path('hybracter')
  output:
      tuple val(meta), path('assembly.fasta'), emit:fasta
	script:
	"""
	if [ -f hybracter/FINAL_OUTPUT/complete/sample_final.fasta ]; then
			ln -s hybracter/FINAL_OUTPUT/complete/sample_final.fasta assembly.fasta
	fi
	"""
} 