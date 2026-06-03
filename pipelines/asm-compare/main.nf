

include { MINIMAP2_ALIGN_ASM5   } from './modules/minimap2/align_asm5'
include { BCFTOOLS_ASM5_MPILEUP } from './modules/bcftools/asm5_mpileup'
include { RMD_RENDER            } from './modules/rmd/render'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'


workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		def query_ch = Channel.fromPath(params.query).map({
				def id = it.name.replaceAll(/\.(fasta|fna|fa)$/,'')
				[[sample_id:id],params.ref,it]
		})

		query_ch 
		| MINIMAP2_ALIGN_ASM5
		
		MINIMAP2_ALIGN_ASM5.out.bam
		| map({meta,bam -> [meta,params.ref,bam]})
		| BCFTOOLS_ASM5_MPILEUP

		RMD_RENDER(
			Channel.empty(),
			file("${projectDir}/assets/report.qmd"),
			file("${projectDir}/assets/report.qmd")
		)

	publish:
		bam = MINIMAP2_ALIGN_ASM5.out.bam
		bai = MINIMAP2_ALIGN_ASM5.out.bai
		mut_vcf = BCFTOOLS_ASM5_MPILEUP.out.vcf
		mut_vcf_csi = BCFTOOLS_ASM5_MPILEUP.out.vcf_csi
		mut_txt = BCFTOOLS_ASM5_MPILEUP.out.txt
		html = RMD_RENDER.out.html
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




