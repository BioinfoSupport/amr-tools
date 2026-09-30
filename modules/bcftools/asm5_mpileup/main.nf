

// Simple variant calling from a BAM between 2 assembled genomes

process BCFTOOLS_ASM5_MPILEUP {
    container 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'
    memory '10 GB'
    cpus 4
    time '1h'
    input:
	    tuple val(meta), path('ref.fasta'), path('alignment.bam'), path(ref_gff)
    output:
	    tuple val(meta), path("mutations.vcf.gz"), emit: vcf
	    tuple val(meta), path("mutations.vcf.gz.csi"), emit: vcf_csi
	    tuple val(meta), path("mutations.txt"), emit: txt
    script:
      def csq_cmd = ref_gff ? "| bcftools csq -Oz --force --local-csq --fasta-ref='ref.fasta' --gff-annot='${ref_gff}'" : ""
      def csq_qry = ref_gff ? "\t%INFO/BCSQ" : ""
	    """
			bcftools mpileup -Ou \
			    --threads ${task.cpus} \
					${task.ext.args?:''} \
			  	--per-sample-mF \
			    -a FORMAT/AD \
			  	--fasta-ref="ref.fasta" \
			  	alignment.bam \
				| bcftools norm -Ou -m- -f "ref.fasta" \
			  | bcftools filter -Ou -i '(FORMAT/AD[:1] >= 1)' \
			  | bcftools +setGT -Oz -- -t a -n c:1 \
			  ${csq_cmd} \
			  > mutations.vcf.gz
			bcftools index mutations.vcf.gz
	    bcftools query -f "%CHROM\t%POS\t%REF>%ALT${csq_qry}" mutations.vcf.gz > mutations.txt
	    """
		stub:
			"""
			touch mutations.vcf.gz mutations.vcf.gz.csi mutations.txt
			"""
}






