
include { MINIMAP2_ALIGN_ONT } from './modules/minimap2/align_ont'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_LONG  } from './modules/samtools/stats'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_SHORT } from './modules/samtools/stats'
include { BWA_MEM            } from './modules/bwa/mem'
include { BWA_INDEX          } from './modules/bwa/index'
include { ORGANIZE_FILES     } from './modules/organize_files'



process COMBINE_QC_STATS {
  container "registry.gitlab.unige.ch/amr-genomics/rscript:main"
  memory '8 GB'
  cpus 2
  time '30 min'
  input:
		tuple val(meta), path("stats")
		path("assets") 
  output:
    tuple val(meta), path('qc.rds'), emit: 'rds'
  script:
		def builder = new groovy.json.JsonBuilder(meta)
	  """
		#!/usr/bin/env Rscript
		source("assets/lib_assembly_stats.R")
		stats <- read_assembly_stats("./stats")
		stats[["meta"]] <- '${builder.toString()}'
		saveRDS(stats,"qc.rds")
		"""
}

process ASSEMBLIES_MULTIQC {
  container "registry.gitlab.unige.ch/amr-genomics/rscript:main"
  memory '8 GB'
  cpus 2
  time '30 min'
  input:
		path("stats/qc*.rds")
		path("assets") 
  output:
    path('multiqc.txt')
  script:
		"""
		#!/usr/bin/env Rscript
		source("assets/lib_assembly_stats.R")
		rlang::inject(aggregate_assembly_stats(!!!(fs::dir_ls("stats",glob="*.rds") |> map(readRDS) |> unname()))) |>
		write.table(file="multiqc.txt",sep="\\t",quote=FALSE,row.names=FALSE)
		"""
}

workflow ASSEMBLIES_QC {
	take:
		fa_ch
		fqs_ch
		fql_ch
	main:
	
		// Short reads alignment and statistics
		BWA_MEM(BWA_INDEX(fa_ch).join(fqs_ch))
		SAMTOOLS_STATS_SHORT(BWA_MEM.out.bam)
		
		// Long reads alignment and statistics
		MINIMAP2_ALIGN_ONT(fa_ch.join(fql_ch))
		SAMTOOLS_STATS_LONG(MINIMAP2_ALIGN_ONT.out.bam)
		
		//TODO: RUN VCF_LONG
		//TODO: RUN VCF_SHORT
		//TODO: HTML_REPORT()
		//TODO: CHARACTERIZE_UNMAPPED_READS
		fa_ch
			.join(SAMTOOLS_STATS_LONG.out,remainder:true)
			.join(SAMTOOLS_STATS_SHORT.out,remainder:true)
			.map({meta,x1,x2,x3 -> [meta,[[x1,"assembly.fasta"],[x2,"long_reads.bam.stats"],[x3,"short_reads.bam.stats"]].findAll({x,y -> x})]})
			| ORGANIZE_FILES
		
		COMBINE_QC_STATS(ORGANIZE_FILES.out,"${moduleDir}/assets")
		ASSEMBLIES_MULTIQC(COMBINE_QC_STATS.out.rds.map({m,x -> x}).collect(),"${moduleDir}/assets")

	emit:
		long_bam        = MINIMAP2_ALIGN_ONT.out.bam
		long_bai        = MINIMAP2_ALIGN_ONT.out.bai
		long_bam_stats  = SAMTOOLS_STATS_LONG.out
		long_vcf        = Channel.empty()
		
		short_bam       = BWA_MEM.out.bam
		short_bai       = BWA_MEM.out.bai
		short_bam_stats = SAMTOOLS_STATS_SHORT.out
		short_vcf       = Channel.empty()
		
		assembly_qc          = COMBINE_QC_STATS.out.rds
		assembly_multiqc_txt = ASSEMBLIES_MULTIQC.out
}








/* Assembly stats
Short Reads
  - % properly paired
  - Number chimeric pair
  - Number of identified mutation in the VCF
Long Reads
  - % mapped
  - Number of identified mutation in the VCF
*/

/* Contigs stats
Short Reads
  - % mapped
  - % properly paired
  - Avg coverage
Long Reads
  - % mapped
  - Avg coverage
*/

/* Report
Short Reads
   - inter contigs links
*/


/*
process IGV_SCRIPT {
	input:
	  val(meta)
  output:
    tuple val(meta), path("load_in_igv.sh")
	script:
	"""
	file("${moduleDir}/assets/load_in_igv.sh")
	"""
}
*/