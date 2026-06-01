

include { MINIMAP2_ALIGN_ASM5 } from './modules/minimap2/align_asm5' 
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'


workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		MINIMAP2_ALIGN_ASM5([[sample_id:'input'],params.ref,params.fasta])

	publish:
		bam = MINIMAP2_ALIGN_ASM5.out.bam
		bai = MINIMAP2_ALIGN_ASM5.out.bai
		vcf = Channel.empty()
		html = Channel.empty()
}

output {

	bam {
		path { m,x -> x >> "minimap2/${m.sample_id}.bam"}
	}
	bai {
		path { m,x -> x >> "minimap2/${m.sample_id}.bai"}
	}
	vcf {
		path { "/"}
	}
	html {
		path { "reports/"}
	}
	
}




