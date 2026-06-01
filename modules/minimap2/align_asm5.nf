process MINIMAP2_ALIGN_ASM5 {
    container 'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c'
    memory '10 GB'
    cpus 4
    time '1h'
    input:
	    tuple val(meta), path('ref.fasta'), path('asm.fasta')
    output:
	    tuple val(meta), path("out.bam"), emit: bam
	    tuple val(meta), path("out.bam.bai"), emit: bai
    script:
	    """
	    minimap2 -x asm5 ${task.ext.args?:''} -t ${task.cpus} ref.fasta asm.fasta -a \
	    | samtools sort -@ ${task.cpus} --write-index -O BAM -o out.bam##idx##out.bam.bai
	    """
		stub:
			"""
			touch out.bam
			touch out.bam.bai
			"""
}




/*
NCBI_FASTA=paper/bc12_illumina_polished.fasta
minimap2 -cx asm5 -z1000000 -p 0.1 -Y --cs -a "${NCBI_FASTA}" "${REF_FASTA}" \
  | samtools sort - \
  | bcftools mpileup -Ou --no-BAQ --indel-size 1000 \
  	--ff=UNMAP --ambig-reads=drop --gap-frac=0 --min-ireads=0 --per-sample-mF \
    -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR \
  	--fasta-ref="${NCBI_FASTA}" - \
	| bcftools norm -Ou -m- -f "${NCBI_FASTA}" \
  | bcftools filter -Ou -i '(FORMAT/AD[:1] >= 1)' \
  | bcftools query -f "%CHROM\t%POS\t%REF>%ALT" \
  > paper/bc12_vs_bc01.txt
#  | bcftools view
  
*/



