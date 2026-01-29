
import nextflow.Channel

class Assemblies {
	def meta_ch
	def fasta
	
	Assemblies(opts,samplesheetToList) {
		this.meta_ch = Channel.empty()
		this.fasta   = Channel.empty()
		if (opts.csv) {
				def ss = Channel.fromList(samplesheetToList(opts.csv,"assets/schema/sheets/assemblies.json"))
					.map({it[0].sample_id = it[0].sample_id?:it[0].assembly_id;it})
				this.meta_ch = ss.map({it[0]})
				this.fasta   = ss.map({[it[0],it[1]]}).filter({it[1]})
		} else {
	  	if (opts.fasta) {
	  		this.fasta = Channel.fromPath(opts.fasta)
	  			.map({
	  				def id = it.name.replaceAll(/\.(fasta|fa)$/,'')
	  				[[sample_id:id,assembly_id:id],it]
	  			})
	  	}
			this.meta_ch = assemblies.fasta.map({it[0]}).unique()
		}
	}

	
	def flat_csv() {
		ss = meta_ch
				.join(assemblies.fasta.ifEmpty(['__sentinel__',null]),failOnDuplicate:true,remainder:true)
				.map({[it[0],it[1]]})
				.filter({it[0]!='__sentinel__'})
		ss.map({[
			assembly_id: it[0].assembly_id,
			sample_id:  it[0].sample_id,
			fasta: it[1]?it[1]:null,
			assembler_name: it[0].assembler_name?:null,
			assembler_readset_id: it[0].assembler_readset_id?:null,
			assembler_parameters: it[0].assembler_parameters?:null
		]})
	}

}
