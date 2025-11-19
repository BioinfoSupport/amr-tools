
import nextflow.Channel

class AmrUtils {

	static def parse_generic_params(opts,samplesheetToList) {
		def ss = [
			fa_ch : Channel.empty(),
			fql_ch: Channel.empty(),
			fqs_ch: Channel.empty(),
			ss_ch : Channel.empty()
		]
		if (opts.samplesheet) {
				ss.ss_ch = Channel.fromList(samplesheetToList(opts.samplesheet))
				ss.ss_ch.view()
				def SS = ss.ss_ch
					.multiMap({x ->
						meta = [sample_id:x[0].sample_id]
						fa_ch: [meta,x[0].assembly_fasta]
						fql_ch: [meta,x[0].long_reads]
						fqs_ch: [meta,[x[0].short_reads_1,x[0].short_reads_2]]
					})
				ss.fa_ch = SS.fa_ch
				ss.fql_ch = SS.fql_ch
				ss.fqs_ch = SS.fqs_ch
		} else {
			if (opts.fasta) {
				ss.fa_ch = Channel.fromPath(opts.fasta)
						.map({x -> [[sample_id:x.name.replaceAll(/\.(fasta|fa|fna)(\.gz)?$/,'')],x]})
			}
			if (opts.long_reads) {
				ss.fql_ch = Channel.fromPath(opts.long_reads)
						.map({x -> [[sample_id:x.name.replaceAll(/\.(fastq\.gz|fq\.gz|bam|cram)$/,'')],x]})
			}
			if (opts.short_reads) {
				ss.fqs_ch = Channel
						.fromFilePairs(opts.short_reads,size:-1) { file -> file.name.replaceAll(/_(R?[12])(_001)?\.(fq|fastq)\.gz$/, '') }
						.map({id,x -> [[sample_id:id],x]})
			}
			def meta_ch = Channel.empty()
					.mix(ss.fa_ch.map({[it[0]]}))
					.mix(ss.fql_ch.map({[it[0]]}))
					.mix(ss.fqs_ch.map({[it[0]]}))
					.unique()
				// DEBUGGING
				/*
				meta_ch
					.combine(ss.fa_ch)
					.view()
					.combine(ss.fql_ch)
					.join(ss.fqs_ch,remainder:true)
				*/
				//.map({m,fa,lr,sr -> m + [assembly_fasta:fa,long_reads:lr,short_reads_1:(sr && sr.size()>0?sr[0]:null),short_reads_2:(sr && sr.size()>1?sr[1]:null)]})
				
		}
		
		// Filter out missing files from the channels
		ss.fa_ch  = ss.fa_ch.filter({x,y -> y})
		ss.fqs_ch = ss.fqs_ch.map({x,y -> [x,y.findAll({v->v})]}).filter({x,y -> y})
		ss.fql_ch = ss.fql_ch.filter({x,y -> y})

		return ss
	}

}


