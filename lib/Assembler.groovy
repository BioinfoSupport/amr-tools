
import nextflow.Channel



class Assembler {
	static Map fromParams(opts) {
		return [name:opts.name?:"long_flye_medaka",args:opts.args?:[:]]
	}
}

/*
class Assemblies {
	def meta_ch
	def fasta
	
	Assemblies() {
		this.meta_ch = Channel.empty()
		this.fasta   = Channel.empty()
	}
	
	static def fromChannels(fasta) {
	  	def x = new Assemblies()
  		x.fasta = fasta
	  	x.meta_ch = fasta.map({it[0]}).unique()
	  	return x
	}

	static def fromCSV(csv_path,samplesheetToList) {
			def x = new Assemblies()
			def ss = Channel.fromList(samplesheetToList(csv_path,"assets/schema/sheets/assemblies.json"))
				.map({it[0].sample_id = it[0].sample_id?:it[0].assembly_id;it})
			x.meta_ch = ss.map({it[0]})
			x.fasta   = ss.map({[it[0],it[1]]}).filter({it[1]})
			return x
	}
	
	static def fromParams(opts,samplesheetToList) {
		this.meta_ch = Channel.empty()
		this.fasta   = Channel.empty()
		if (opts.csv) {
			return fromCSV(opts.csv,samplesheetToList)
		} else {
			def fa_ch = Channel.empty()
	  	if (opts.fasta) {
	  		fa_ch = Channel.fromPath(opts.fasta).map({
	  				def id = it.name.replaceAll(/\.(fasta|fa)$/,'')
	  				[[sample_id:id,assembly_id:id],it]
  			})
	  	}
	  	return fromChannels(fa_ch)
		}
	}

	
	def flat_csv() {
		ss = this.meta_ch
				.join(this.fasta.ifEmpty(['__sentinel__',null]),failOnDuplicate:true,remainder:true)
				.map({[it[0],it[1]]})
				.filter({it[0]!='__sentinel__'})
		ss.map({[
			assembly_id: it[0].assembly_id,
			sample_id:  it[0].sample_id,
			fasta: it[1]?it[1]:null,
			assembler_name: it[0].assembler_name?:null,
			assembler_readset_id: it[0].assembler_readset_id?:null,
			assembler_args: it[0].assembler_args?:null
		]})
	}

}


*/