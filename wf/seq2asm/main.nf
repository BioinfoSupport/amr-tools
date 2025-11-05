

include { UNICYCLER as LONG_UNICYCLER                    } from './subworkflows/unicycler'
include { UNICYCLER as SHORT_UNICYCLER                   } from './subworkflows/unicycler'
include { UNICYCLER as HYBRID_UNICYCLER                  } from './subworkflows/unicycler'
include { HYBRACTER as LONG_HYBRACTER                    } from './subworkflows/hybracter'
include { HYBRACTER as HYBRID_HYBRACTER                  } from './subworkflows/hybracter'
include { SPADES    as SHORT_SPADES                      } from './subworkflows/spades'
include { FLYE_MEDAKA_PILON as HYBRID_FLYE_MEDAKA_PILON  } from './subworkflows/flye_medaka_pilon'
include { FLYE_MEDAKA_PILON as LONG_FLYE_MEDAKA          } from './subworkflows/flye_medaka_pilon'

workflow SEQ2ASM {
	take:
		opts
		fqs_ch    // channel: [ val(meta), path(short_reads) ]
		fql_ch    // channel: [ val(meta), path(long_reads) ]
	main:

		// Short reads only assemblies
		SHORT_SPADES(fqs_ch.filter({opts.short_spades}),Channel.empty())
		SHORT_UNICYCLER(fqs_ch.filter({opts.short_unicycler}),Channel.empty())
		
		// Long reads only assemblies
		LONG_HYBRACTER(Channel.empty(),fql_ch.filter({opts.long_hybracter}))
		LONG_UNICYCLER(Channel.empty(),fql_ch.filter({opts.long_unicycler}))
		LONG_FLYE_MEDAKA(Channel.empty(),fql_ch.filter({opts.long_flye_medaka}))

		// Hybrid assemblies
		HYBRID_HYBRACTER(fqs_ch.filter({opts.hybrid_hybracter}),fql_ch.filter({opts.hybrid_hybracter}))
		HYBRID_UNICYCLER(fqs_ch.filter({opts.hybrid_unicycler}),fql_ch.filter({opts.hybrid_unicycler}))
		HYBRID_FLYE_MEDAKA_PILON(fqs_ch.filter({opts.hybrid_flye_medaka_pilon}),fql_ch.filter({opts.hybrid_flye_medaka_pilon}))
		
	emit:
		fasta = Channel.empty().mix(
			SHORT_SPADES.out.fasta.map({meta,x -> [meta,[assembly_name:'short_spades'],x]}),
			SHORT_UNICYCLER.out.fasta.map({meta,x -> [meta,[assembly_name:'short_unicycler'],x]}),
			LONG_FLYE_MEDAKA.out.fasta.map({meta,x -> [meta,[assembly_name:'long_flye_medaka'],x]}),
			LONG_UNICYCLER.out.fasta.map({meta,x -> [meta,[assembly_name:'long_unicycler'],x]}),
			LONG_HYBRACTER.out.fasta.map({meta,x -> [meta,[assembly_name:'long_hybracter'],x]}),
		  HYBRID_UNICYCLER.out.fasta.map({meta,x -> [meta,[assembly_name:'hybrid_unicycler'],x]}),
		  HYBRID_HYBRACTER.out.fasta.map({meta,x -> [meta,[assembly_name:'hybrid_hybracter'],x]}),
		  HYBRID_FLYE_MEDAKA_PILON.out.fasta.map({meta,x -> [meta,[assembly_name:'hybrid_flye_medaka_pilon'],x]})
		)
		dir = Channel.empty().mix(
			SHORT_SPADES.out.dir.map({meta,x -> [meta,[assembly_name:'short_spades'],x]}),
			SHORT_UNICYCLER.out.dir.map({meta,x -> [meta,[assembly_name:'short_unicycler'],x]}),
			LONG_FLYE_MEDAKA.out.dir.map({meta,x -> [meta,[assembly_name:'long_flye_medaka'],x]}),
			LONG_UNICYCLER.out.dir.map({meta,x -> [meta,[assembly_name:'long_unicycler'],x]}),
			LONG_HYBRACTER.out.dir.map({meta,x -> [meta,[assembly_name:'long_hybracter'],x]}),
		  HYBRID_UNICYCLER.out.dir.map({meta,x -> [meta,[assembly_name:'hybrid_unicycler'],x]}),
		  HYBRID_HYBRACTER.out.dir.map({meta,x -> [meta,[assembly_name:'hybrid_hybracter'],x]}),
		  HYBRID_FLYE_MEDAKA_PILON.out.dir.map({meta,x -> [meta,[assembly_name:'hybrid_flye_medaka_pilon'],x]})
		)
}




// ------------------------------------------------------------------
// Main entry point when running the pipeline from command line
// ------------------------------------------------------------------
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

params.samplesheet = null
params.long_reads  = []
params.short_reads = []
params.long_unicycler           = false,
params.long_hybracter           = false,
params.long_flye_medaka         = true,
params.short_spades             = false,
params.short_unicycler          = false,
params.hybrid_unicycler         = false,
params.hybrid_hybracter         = false,
params.hybrid_flye_medaka_pilon = false

workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		// Extract Long_read and Short_read channels from params 
		fa_ch = Channel.empty()
		lr_ch = Channel.empty()
		sr_ch = Channel.empty()
		if (params.samplesheet) {
				SS = Channel.fromList(samplesheetToList(params.samplesheet, "assets/schema_samplesheet.json"))
					.multiMap({x ->
						meta = [sample_id:x[0].sample_id,assembly_name:x[0].assembly_name]
						fa_ch: [meta,x[0].assembly_fasta]
						lr_ch: [meta,x[0].long_reads]
						sr_ch: [meta,[x[0].short_reads_1,x[0].short_reads_2]]
					})
				fa_ch = SS.fa_ch
				lr_ch = SS.lr_ch
				sr_ch = SS.sr_ch
		} else {
			fa_ch = Channel.fromPath(params.fasta)
				.map({x -> [[sample_id:x.name.replaceAll(/\.(fasta|fa|fna)(\.gz)?$/,''),assembly_name:'fasta'],x]})
		}
		
		// Filter out missing values
		fa_ch = fa_ch.filter({x,y -> y})
		

		SEQ2ASM(params,fqs_ch,fql_ch)

	publish:
		orgfinder           = AMR_ANNOT.out.orgfinder
		amrfinderplus       = AMR_ANNOT.out.amrfinderplus
		resfinder           = AMR_ANNOT.out.resfinder
		mobtyper            = AMR_ANNOT.out.mobtyper
		plasmidfinder       = AMR_ANNOT.out.plasmidfinder
		cgemlst             = AMR_ANNOT.out.cgemlst
		MLST                = AMR_ANNOT.out.MLST
		prokka              = AMR_ANNOT.out.prokka
		resfinder_long      = AMR_ANNOT.out.resfinder_long
		resfinder_short     = AMR_ANNOT.out.resfinder_short
		plasmidfinder_long  = AMR_ANNOT.out.plasmidfinder_long
		plasmidfinder_short = AMR_ANNOT.out.plasmidfinder_short
		
  	multireport_html    = AMR_ANNOT.out.html_report
  	multireport_xlsx    = AMR_ANNOT.out.xlsx_report
}

output {
	orgfinder {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/orgfinder"}
	}
	amrfinderplus {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/amrfinderplus"}
	}
	resfinder {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/resfinder"}
	}
	mobtyper {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/mobtyper"}
	}
	plasmidfinder {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/plasmidfinder"}
	}
	cgemlst {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/cgemlst"}
	}
	MLST {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/MLST"}
	}
	prokka {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/prokka"}
	}
	
	resfinder_long {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/resfinder_long"}
	}
	resfinder_short {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/resfinder_short"}
	}
	plasmidfinder_long {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/plasmidfinder_long"}
	}
	plasmidfinder_short {
		path { m,x -> x >> "samples/${m.sample_id}/assemblies/${m.assembly_name}/plasmidfinder_short"}
	}
	
	multireport_html {
		path { it >> "./amr_annot.html" }
	}
	multireport_xlsx {
		path { it >> "./amr_annot.xlsx" }
	}
}













