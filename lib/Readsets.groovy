
import nextflow.Channel
import groovy.transform.TupleConstructor



class Readsets {
	def meta_ch
	def short_reads
	def long_reads
	
	Readsets() {
		this.meta_ch     = Channel.empty()
		this.long_reads  = Channel.empty()
		this.short_reads = Channel.empty()
	}
	
	static def fromChannels(short_reads,long_reads) {
	  	def x = new Readsets()
	  	x.short_reads = short_reads
	  	x.long_reads  = long_reads
	  	x.meta_ch     = x.long_reads.map({it[0]}).mix(x.short_reads.map({it[0]})).unique()
	  	return x
	}
	
	static def fromCSV(csv_path,samplesheetToList) {
			def x = new Readsets()
			def ss = Channel.fromList(samplesheetToList(csv_path,"assets/schema/sheets/readsets.json"))
				.map({it[0].sample_id = it[0].sample_id?:it[0].readset_id;it})
			x.meta_ch     = ss.map({it[0]})
			x.short_reads = ss.map({[it[0],it[1]]}).filter({it[1]})
			x.long_reads  = ss.map({[it[0],[it[2],it[3]]]}).filter({it[1].findAll({it})})
			return x
	}

	static def fromParams(opts,samplesheetToList) {
		if (opts.csv) {
			return fromCSV(opts.csv,samplesheetToList)
		} else {
			def sr_ch = Channel.empty()			
			def lr_ch = Channel.empty()
	  	if (opts.long_reads) {
	  		lr_ch = Channel.fromPath(opts.long_reads).map({
	  				def id = it.name.replaceAll(/\.(fastq\.gz|fq\.gz|bam|cram)$/,'')
	  				[[sample_id:id,readset_id:id],it]
	  		})
	  	}
	  	if (opts.short_reads) {
				sr_ch = Channel.fromFilePairs(opts.short_reads,size:-1) { 
						file -> file.name.replaceAll(/_(R?[12])(_001)?\.(fq|fastq)\.gz$/, '') 
				}
				.map({id,x -> [[sample_id:id,readset_id:id],x]})
	  	}
	  	return fromChannels(sr_ch,lr_ch)
		}
	}
	
	
	def flat_csv() {
		def ss = this.meta_ch
			.join(this.long_reads.ifEmpty(['__sentinel__',null]),failOnDuplicate:true,remainder:true)
			.join(this.short_reads.ifEmpty(['__sentinel__',null]),failOnDuplicate:true,remainder:true)
			.map({[it[0],it[1],it[2]?it[2][0]:null,it[2]?it[2][1]:null]})
			.filter({it[0]!='__sentinel__'})
		ss.map({[
			readset_id: it[0].readset_id,
			sample_id:  it[0].sample_id,
			long_reads: it[1]?it[1].toString():null,
			short_reads_1: it[2]?it[2].toString():null,
			short_reads_2: it[3]?it[3].toString():null
		]})
	}

}

