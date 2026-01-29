
import nextflow.Channel

class Readsets {
	def meta_ch
	def long_reads
	def short_reads
	
	Readsets(opts,samplesheetToList) {
		this.meta_ch     = Channel.empty()
		this.long_reads  = Channel.empty()
		this.short_reads = Channel.empty()
		if (opts.csv) {
			// Directly load the sheet from CSV file
			def ss = Channel.fromList(samplesheetToList(opts.csv,"assets/schema/sheets/readsets.json"))
				.map({it[0].sample_id = it[0].sample_id?:it[0].readset_id;it})
			this.meta_ch     = ss.map({it[0]})
			this.long_reads  = ss.map({[it[0],it[1]]}).filter({it[1]})
			this.short_reads = ss.map({[it[0],[it[2],it[3]]]}).filter({it[1].findAll({it})})
		} else {
			// Generate the sheet from parameters
	  	if (opts.long_reads) {
	  		this.long_reads = Channel.fromPath(opts.long_reads).map({
	  				def id = it.name.replaceAll(/\.(fastq\.gz|fq\.gz|bam|cram)$/,'')
	  				[[sample_id:id,readset_id:id],it]
	  			})
	  	}
	  	if (opts.short_reads) {
				this.short_reads = Channel.fromFilePairs(opts.short_reads,size:-1) { 
						file -> file.name.replaceAll(/_(R?[12])(_001)?\.(fq|fastq)\.gz$/, '') 
				}
				.map({id,x -> [[sample_id:id,readset_id:id],x]})
	  	}
	  	this.meta_ch = this.long_reads.map({it[0]}).mix(this.short_reads.map({it[0]})).unique()
		}
	}
	
	def flat_csv_channel() {
		def ss = this.meta_ch
			.join(this.long_reads.ifEmpty(['__sentinel__',null]),failOnDuplicate:true,remainder:true)
			.join(this.short_reads.ifEmpty(['__sentinel__',null]),failOnDuplicate:true,remainder:true)
			.map({[it[0],it[1],it[2]?it[2][0]:null,it[2]?it[2][1]:null]})
			.filter({it[0]!='__sentinel__'})
		ss.map({[
			readset_id: it[0].readset_id,
			sample_id:  it[0].sample_id,
			long_reads: it[1]?it[1]:null,
			short_reads_1: it[2]?it[2]:null,
			short_reads_2: it[3]?it[3]:null
		]})
	}

}

