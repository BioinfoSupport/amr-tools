
import nextflow.Channel

class AmrUtils {
	
	static def get_readsets(opts,samplesheetToList) {
		def readsets = [
      long_reads : Channel.empty(),
      short_reads: Channel.empty(),
      ss         : Channel.empty(),
      csv        : Channel.empty()
		]
		
		if (opts.csv) {
				readsets.ss = Channel.fromList(samplesheetToList(opts.csv,"assets/schema/sheets/readsets.json"))
					.map({it[0].sample_id = it[0].sample_id?:it[0].readset_id;it})
				readsets.long_reads  = readsets.ss.map({[it[0],it[1]]}).filter({it[1]})
				readsets.short_reads = readsets.ss.map({[it[0],[it[2],it[3]]]}).filter({it[1].findAll({it})})
		} else {
	  	if (opts.long_reads) {
	  		readsets.long_reads = Channel.fromPath(opts.long_reads)
	  			.map({
	  				def id = it.name.replaceAll(/\.(fastq\.gz|fq\.gz|bam|cram)$/,'')
	  				[[sample_id:id,readset_id:id],it]
	  			})
	  	}
	  	if (opts.short_reads) {
				readsets.short_reads = Channel.fromFilePairs(opts.short_reads,size:-1) { 
						file -> file.name.replaceAll(/_(R?[12])(_001)?\.(fq|fastq)\.gz$/, '') 
					}
					.map({id,x -> [[sample_id:id,readset_id:id],x]})
	  	}
			def meta_ch = readsets.long_reads.map({it[0]}).mix(readsets.short_reads.map({it[0]})).unique()
			readsets.ss = meta_ch
				.join(readsets.long_reads.ifEmpty(['__sentinel__',null]),failOnDuplicate:true,remainder:true)
				.join(readsets.short_reads.ifEmpty(['__sentinel__',null]),failOnDuplicate:true,remainder:true)
				.map({[it[0],it[1],it[2]?it[2][0]:null,it[2]?it[2][1]:null]})
				.filter({it[0]!='__sentinel__'})
		}
		
		readsets.csv = readsets.ss.map({[
			readset_id: it[0].readset_id,
			sample_id:  it[0].sample_id,
			long_reads: it[1]?it[1]:null,
			short_reads_1: it[2]?it[2]:null,
			short_reads_2: it[3]?it[3]:null
		]})

		return readsets
	}

	static def get_assemblies(opts,samplesheetToList) {
		def assemblies = [
      fasta : Channel.empty(),
      ss    : Channel.empty(),
      csv   : Channel.empty()
		]
		
		if (opts.csv) {
				assemblies.ss = Channel.fromList(samplesheetToList(opts.csv,"assets/schema/sheets/assemblies.json"))
					.map({it[0].sample_id = it[0].sample_id?:it[0].assembly_id;it})
				assemblies.fasta  = assemblies.ss.map({[it[0],it[1]]}).filter({it[1]})
		} else {
	  	if (opts.fasta) {
	  		assemblies.fasta = Channel.fromPath(opts.fasta)
	  			.map({
	  				def id = it.name.replaceAll(/\.(fasta|fa)$/,'')
	  				[[sample_id:id,assembly_id:id],it]
	  			})
	  	}
			def meta_ch = assemblies.fasta.map({it[0]}).unique()
			assemblies.ss = meta_ch
				.join(assemblies.fasta.ifEmpty(['__sentinel__',null]),failOnDuplicate:true,remainder:true)
				.map({[it[0],it[1]]})
				.filter({it[0]!='__sentinel__'})
		}
		
		assemblies.csv = assemblies.ss.map({[
			assembly_id: it[0].assembly_id,
			sample_id:  it[0].sample_id,
			fasta: it[1]?it[1]:null,
			assembler_name: it[0].assembler_name?:null,
			assembler_readset_id: it[0].assembler_readset_id?:null,
			assembler_parameters: it[0].assembler_parameters?:null
		]})

		return assemblies
	}
	
	

	static def parse_generic_params(opts,samplesheetToList) {
		def ss = [
			fa_ch : Channel.empty(),
			fql_ch: Channel.empty(),
			fqs_ch: Channel.empty(),
			ss_ch : Channel.empty()
		]
		if (opts.samplesheet) {
				ss.ss_ch = Channel.fromList(samplesheetToList(opts.samplesheet))
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


