

//include { FQ_SUBSAMPLE      } from './modules/fq_subsample'
include { SAMTOOLS_FASTQ    } from './modules/samtools/fastq'
include { SEQUENCING_QC     } from './workflows/sequencing_qc'
include { ASSEMBLE_READS    } from './workflows/assemble_reads'
include { ASSEMBLY_QC       } from './workflows/assembly_qc'
include { ANNOTATE_ASSEMBLY } from './workflows/annotate_assembly'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'


def get_samplesheet() {
	ss = [
		lr_ch: Channel.empty(),
		sr_ch: Channel.empty()
	]
	if (params.samplesheet) {
		SS = Channel.fromList(samplesheetToList(params.samplesheet, "assets/schema_samplesheet.json"))
			.multiMap({x ->
				lr_ch: [[sample_id:x[0].sample_id],x[0].long_reads]
				sr_ch: [[sample_id:x[0].sample_id],[x[0].short_reads_1,x[0].short_reads_2]]
			})
		ss.lr_ch = SS.lr_ch
		ss.sr_ch = SS.sr_ch
	} else {
		if (params.long_reads) {
			ss.lr_ch = Channel.fromPath(params.long_reads)
					.map({x -> tuple(["sample_id":x.name.replaceAll(/\.(fastq\.gz|fq\.gz|bam|cram)$/,'')],x)})
		}
		if (params.short_reads) {
			ss.sr_ch = Channel
					.fromFilePairs(params.short_reads,size:-1) { file -> file.name.replaceAll(/_(R?[12])(_001)?\.(fq|fastq)\.gz$/, '') }
					.map({id,x -> [["sample_id":id],x]})
		}	
	}
	// Filter missing values
	ss.sr_ch = ss.sr_ch.map({x,y -> [x,y.findAll({v->v})]}).filter({x,y -> y})
	ss.lr_ch = ss.lr_ch.filter({x,y -> y})
	return ss	
}


workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		// Prepare SampleSheet
		ss = get_samplesheet()
		
		// CONVERT long_reads given in BAM/CRAM format into FASTQ format
		ss.lr_ch = ss.lr_ch.branch({meta,f -> 
			bam: f.name =~ /\.(bam|cram)$/
			fq: true
		})
		ss.lr_ch = ss.lr_ch.fq.mix(SAMTOOLS_FASTQ(ss.lr_ch.bam))

		// Reduce FASTQ size if needed
		//ss.lr_ch = FQ_SUBSAMPLE(ss.lr_ch)

		// Reads Quality Controls, get a multiQC
		SEQUENCING_QC(
			ss.lr_ch.filter({!params.skip_seq_qc}),
			ss.sr_ch.filter({!params.skip_seq_qc})
		)

		// de novo read assembly
    //ASSEMBLE_READS(params.assembler,ss.sr_ch,ss.lr_ch)
    ASSEMBLE_READS = [out: [fasta: Channel.empty(),dir: Channel.empty()]]
    
		// Assembly QC
		/*
		meta_map = ASSEMBLE_READS.out.fasta.map({m1,m2,fa -> [m1,m1+m2]})
		ASSEMBLY_QC(
			ASSEMBLE_READS.out.fasta.map({m1,m2,fa -> [m1+m2,fa]}),
			meta_map.combine(ss.sr_ch,by:0).map({m1,m2,x -> [m2,x]}),
			meta_map.combine(ss.lr_ch,by:0).map({m1,m2,x -> [m2,x]})
		)
		*/
		ASSEMBLY_QC = [out: [assembly_qc: Channel.empty(),long_bam: Channel.empty(),long_bai: Channel.empty(),long_bam_stats: Channel.empty(),short_bam: Channel.empty(),short_bai: Channel.empty(),short_bam_stats: Channel.empty(),assembly_multiqc: Channel.empty()]]
		//ANNOTATE_ASSEMBLY(params.annot,Channel.empty(),ss.lr_ch,ss.sr_ch)	
		

	publish:
		asm_fasta = ASSEMBLE_READS.out.fasta
		asm_dir = ASSEMBLE_READS.out.dir
		
		asm_qc  = ASSEMBLY_QC.out.assembly_qc
		asm_qc_long_bam = ASSEMBLY_QC.out.long_bam
		asm_qc_long_bai = ASSEMBLY_QC.out.long_bai
		asm_qc_long_bam_stats = ASSEMBLY_QC.out.long_bam_stats		
		asm_qc_short_bam = ASSEMBLY_QC.out.short_bam
		asm_qc_short_bai = ASSEMBLY_QC.out.short_bai
		asm_qc_short_bam_stats = ASSEMBLY_QC.out.short_bam_stats
		
		seq_multiqc = Channel.empty() //SEQUENCING_QC.out.multiqc_html
		asm_multiqc = ASSEMBLY_QC.out.assembly_multiqc
}



output {
	asm_fasta {
		path { m1,m2,x -> x >> "samples/${m1.sample_id}/assemblies/${m2.assembly_name}/assembly.fasta"}
	}
	asm_dir {
		path { m1,m2,x -> x >> "samples/${m1.sample_id}/assemblies/${m2.assembly_name}/assembler_output"}
	}
	asm_qc_long_bam {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/long_reads.bam"}
	}
	asm_qc_long_bai {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/long_reads.bam.bai"}
	}
	asm_qc_long_bam_stats {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/long_reads.bam.stats"}
	}
	asm_qc_short_bam {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/short_reads.bam"}
	}
	asm_qc_short_bai {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/short_reads.bam.bai"}
	}
	asm_qc_short_bam_stats {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/short_reads.bam.stats"}
	}
	asm_qc {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/assembly_qc.rds"}
	}
	
	seq_multiqc {
		path { it >> "./multiqc_sequencing.html" }
	}
	asm_multiqc {
		path { it >> "./multiqc_assembly.txt"}
	}
}
