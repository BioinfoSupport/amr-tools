

// Simple variant calling from a BAM between 2 assembled genomes

process BCFTOOLS_ASM5_MPILEUP {
    container 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'
    memory '10 GB'
    cpus 4
    time '1h'
    input:
	    tuple val(meta), path('ref.fasta'), path('alignment.bam')
    output:
	    tuple val(meta), path("mutations.vcf.gz"), emit: vcf
	    tuple val(meta), path("mutations.vcf.gz.csi"), emit: vcf_csi
	    tuple val(meta), path("mutations.txt"), emit: txt
    script:
	    """
			bcftools mpileup -Ou \
			    --threads ${task.cpus} \
					${task.ext.args?:''} \
			  	--per-sample-mF \
			    -a FORMAT/AD \
			  	--fasta-ref="ref.fasta" \
			  	alignment.bam \
				| bcftools norm -Ou -m- -f "ref.fasta" \
			  | bcftools filter -Oz -i '(FORMAT/AD[:1] >= 1)' \
			  > mutations.vcf.gz
			bcftools index mutations.vcf.gz
	    bcftools query -f "%CHROM\t%POS\t%REF>%ALT" mutations.vcf.gz > mutations.txt
	    """
		stub:
			"""
			touch mutations.vcf.gz mutations.vcf.gz.csi mutations.txt
			"""
}






