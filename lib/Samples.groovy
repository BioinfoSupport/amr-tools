
import nextflow.Channel

class Samples {
	
	static def fromCSV(csv_path,samplesheetToList) {
		return Channel.fromList(samplesheetToList(csv_path,"assets/schema/samples_sheets_schema.json"))
	}

	static def fromChannels(sr_ch,lr_ch) {
		return lr_ch.ifEmpty(['__sentinel__',null])
			.join(sr_ch.ifEmpty(['__sentinel__',null]),failOnDuplicate:true,remainder:true)
			.filter({it[0]!='__sentinel__'})
			.map({it[0] + [long_reads:it[1],short_reads_1:it[2]?it[2][0]:null,short_reads_2:it[2]?it[2][1]:null]})
			.map({it + [sample_id:it.sample_id?:it.readset_id]})
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
	  				[[readset_id:id],it]
	  		})
	  	}
	  	if (opts.short_reads) {
				sr_ch = Channel.fromFilePairs(opts.short_reads,size:-1) { 
						file -> file.name.replaceAll(/_(R?[12])(_001)?\.(fq|fastq)\.gz$/, '') 
				}
				.map({id,x -> [[readset_id:id],x]})
	  	}
			return fromChannels(sr_ch,lr_ch)			
		}
	}
	
	static def toCSV(readsets) {
		return readsets
	}
	
	static def long_reads_channel(readsets) {
		return readsets
			.filter({it.long_reads})
			.map({[it.subMap('sample_id'),it.long_reads]})
	}

	static def short_reads_channel(readsets) {
		return readsets
			.filter({it.short_reads_1})
			.map({[it.subMap('sample_id'),[it.short_reads_1,it.short_reads_2].findAll({it})]})
	}
}


