

include { MINIMAP2_ALIGN_ASM5 } from './modules/minimap2/align_asm5' 
include { BCFTOOLS_ASM5_MPILEUP } from './modules/bcftools/asm5_mpileup' 
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'


workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		MINIMAP2_ALIGN_ASM5([[sample_id:'input'],params.ref,params.query])
		
		MINIMAP2_ALIGN_ASM5.out.bam
		| map({meta,bam -> [meta,params.ref,bam]})
		| BCFTOOLS_ASM5_MPILEUP

	publish:
		bam = MINIMAP2_ALIGN_ASM5.out.bam
		bai = MINIMAP2_ALIGN_ASM5.out.bai
		mut_vcf = BCFTOOLS_ASM5_MPILEUP.out.vcf
		mut_vcf_csi = BCFTOOLS_ASM5_MPILEUP.out.vcf_csi
		mut_txt = BCFTOOLS_ASM5_MPILEUP.out.txt
		html = Channel.empty()
}

output {

	bam {
		path { m,x -> x >> "${m.sample_id}.bam"}
	}
	bai {
		path { m,x -> x >> "${m.sample_id}.bai"}
	}
	mut_vcf {
		path { m,x -> x >> "${m.sample_id}.vcf.gz"}
	}
	mut_vcf_csi {
		path { m,x -> x >> "${m.sample_id}.vcf.gz.csi"}
	}
	mut_txt {
		path { m,x -> x >> "${m.sample_id}.txt"}
	}	
	html {
		path { "reports/"}
	}
	
}




