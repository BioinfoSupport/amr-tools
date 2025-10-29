

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



