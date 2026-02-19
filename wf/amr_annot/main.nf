#!/usr/bin/env nextflow

include { GZIP_DECOMPRESS as GZIP_DECOMPRESS_FASTA } from './modules/gzip'
include { RESFINDER                                } from './modules/cgetools/resfinder'
include { RESFINDER     as RESFINDER_LONGREAD      } from './modules/cgetools/resfinder'
include { RESFINDER     as RESFINDER_SHORTREAD     } from './modules/cgetools/resfinder'
include { PLASMIDFINDER                            } from './modules/cgetools/plasmidfinder'
include { PLASMIDFINDER as PLASMIDFINDER_LONGREAD  } from './modules/cgetools/plasmidfinder'
include { PLASMIDFINDER as PLASMIDFINDER_SHORTREAD } from './modules/cgetools/plasmidfinder'
include { ORGFINDER_DETECT                         } from './modules/orgfinder/detect'
include { SPECIATOR                                } from './modules/speciator'
include { AMRFINDERPLUS_UPDATE                     } from './modules/amrfinderplus/update'
include { AMRFINDERPLUS_RUN                        } from './modules/amrfinderplus/run'
include { PROKKA                                   } from './modules/tseemann/prokka'
include { MLST                                     } from './modules/tseemann/mlst'
include { MLST as CGEMLST                          } from './modules/cgetools/mlst'
include { MOBTYPER                                 } from './modules/mobsuite/mobtyper.nf'
include { SAMTOOLS_FAIDX                           } from './modules/samtools/faidx'
include { MULTIREPORT                              } from './subworkflows/multireport'



defaultOrgArgs = new groovy.json.JsonSlurper().parseText(file("${moduleDir}/assets/default_org_args.json").text)

// Get name of the organism to use for given sample
def org_name(meta,opts) {
	if (opts.org_name!=null) return opts.org_name
  if (meta.containsKey('org_name')) return meta['org_name']
	return null
}

// Retreive arguments to use for the given tool on the given sample
def tool_args(tool_name,meta,opts,org_name=null) {
	def key = tool_name + "_args"
	def default_args_key = 'default_' + key
	if (opts.defaults.containsKey(key)) return opts[key]
  if (meta[1].containsKey(key)) return meta[key]
  if (org_name==null) return opts[default_args_key]
  def org_args = defaultOrgArgs.containsKey(org_name)?defaultOrgArgs[org_name] : [:]
  if (org_args.containsKey(key)) return org_args[key]
  return opts[default_args_key]
}



workflow AMR_ANNOT_ASSEMBLY_FOR_ORG {
	take:
		opts
		fa_ch                // channel: [ val(meta), path(assembly_fna) ]
		detected_orgname_ch  // channel: [ val(meta), val(orgname) ]
	main:
		// Update fa_ch with appropriate org_name
		fa_org_ch = fa_ch
			.join(detected_orgname_ch,remainder:true)
			.map({meta,fa,detected_org_name -> [meta,fa,org_name(meta,opts)?:detected_org_name]})

		// MLST typing
		cgemlst_ch = fa_org_ch
			.filter({!opts.skip.cgemlst})
			.map({meta,fa,org_name -> [meta, fa, tool_args('cgemlst',meta,opts,org_name)]})
			.filter({meta,fasta,args -> args!=null})
			| CGEMLST

	emit:
		cgemlst = cgemlst_ch
		org_name = fa_org_ch.map({meta,fa,org_name -> [meta,org_name]})
}



workflow AMR_ANNOT_ASSEMBLY {
	take:
		opts
		asm_ch
	main:

		// CGE - RESFINDER
		def resfinder_ch = RESFINDER(asm_ch.filter({!opts.skip.resfinder}),'fasta')

		// Plasmid typing
		def plasmidfinder_ch = asm_ch.filter({!opts.skip.plasmidfinder}) | PLASMIDFINDER

		// NCBI AMRfinder+
		def amrfinderplus_ch = Channel.empty()
		if (!opts.skip.amrfinderplus) {
			def amrfinderplus_db = AMRFINDERPLUS_UPDATE()
			amrfinderplus_ch = AMRFINDERPLUS_RUN(asm_ch,amrfinderplus_db)
		}

		// MOBsuite - MOBtyper
		def mobtyper_ch = Channel.empty()
		if (!opts.skip.mobtyper) {
			mobtyper_ch = asm_ch | MOBTYPER
		}

		// PROKKA annotations
		def prokka_ch = Channel.empty()
		if (!opts.skip.prokka) {
			prokka_ch = asm_ch | PROKKA
		}

		// Speciator
		def speciator_ch = Channel.empty()
		if (!opts.skip.speciator) {
			speciator_ch = asm_ch | SPECIATOR	
		}
		
		def MLST_ch = Channel.empty()
		if (!opts.skip.MLST) {
			MLST_ch = asm_ch | MLST
		}

    // Run orgfinder to auto detect organism
    def orgfinder_dir_ch = Channel.empty()
    if (!opts.skip.orgfinder) {
    	orgfinder_dir_ch = ORGFINDER_DETECT(asm_ch).orgfinder	
    }
			
		// ---------------------------------------------------------------------
		// Organism specific tools
		// ---------------------------------------------------------------------
		//ann_ch = AMR_ANNOT_ASSEMBLY_FOR_ORG(opts,asm_ch,orgfinder_ch.org_name)

	emit:
			orgfinder           = orgfinder_dir_ch
			amrfinderplus       = amrfinderplus_ch
			resfinder           = resfinder_ch
			mobtyper            = mobtyper_ch
			plasmidfinder       = plasmidfinder_ch
			prokka              = prokka_ch
			org_name            = Channel.empty() //ann_ch.org_name			
			cgemlst             = Channel.empty() //ann_ch.cgemlst
			MLST                = MLST_ch
			
}



workflow AMR_ANNOT_READS {
	take:
		opts
		fqs_ch
		fql_ch
	main:
		PLASMIDFINDER_LONGREAD(fql_ch.filter({!opts.skip.plasmidfinder_longread}))
		RESFINDER_LONGREAD(fql_ch.filter({!opts.skip.resfinder_longread}),'nanopore')
		PLASMIDFINDER_SHORTREAD(fqs_ch.filter({!opts.skip.plasmidfinder_shortread}))
		RESFINDER_SHORTREAD(fqs_ch.filter({!opts.skip.resfinder_shortread}),'illumina')
	emit:
			resfinder_long      = RESFINDER_LONGREAD.out
			resfinder_short     = RESFINDER_SHORTREAD.out
			plasmidfinder_long  = PLASMIDFINDER_LONGREAD.out
			plasmidfinder_short = PLASMIDFINDER_SHORTREAD.out
}




workflow AMR_ANNOT {
	take:
		opts
		asm_ch
		fqs_ch
		fql_ch
	main:
		fai_ch = SAMTOOLS_FAIDX(asm_ch)
		AMR_ANNOT_ASSEMBLY(opts,asm_ch)
		AMR_ANNOT_READS(opts,fqs_ch,fql_ch)
		/*
		MULTIREPORT(
			asm_ch,
			fai_ch,
			AMR_ANNOT_ASSEMBLY.out.org_name,
			AMR_ANNOT_ASSEMBLY.out.orgfinder,
			AMR_ANNOT_ASSEMBLY.out.amrfinderplus,
			AMR_ANNOT_ASSEMBLY.out.resfinder,
			AMR_ANNOT_ASSEMBLY.out.mobtyper,
			AMR_ANNOT_ASSEMBLY.out.plasmidfinder,
			AMR_ANNOT_ASSEMBLY.out.cgemlst,
			AMR_ANNOT_ASSEMBLY.out.MLST,
			AMR_ANNOT_ASSEMBLY.out.prokka,
			AMR_ANNOT_READS.out.resfinder_long,
			AMR_ANNOT_READS.out.plasmidfinder_long,
			AMR_ANNOT_READS.out.resfinder_short,
			AMR_ANNOT_READS.out.plasmidfinder_short
		)
		*/
	emit:
			orgfinder           = AMR_ANNOT_ASSEMBLY.out.orgfinder
			amrfinderplus       = AMR_ANNOT_ASSEMBLY.out.amrfinderplus
			resfinder           = AMR_ANNOT_ASSEMBLY.out.resfinder
			mobtyper            = AMR_ANNOT_ASSEMBLY.out.mobtyper
			plasmidfinder       = AMR_ANNOT_ASSEMBLY.out.plasmidfinder
			cgemlst             = AMR_ANNOT_ASSEMBLY.out.cgemlst
			MLST                = AMR_ANNOT_ASSEMBLY.out.MLST
			prokka              = AMR_ANNOT_ASSEMBLY.out.prokka
			resfinder_long      = AMR_ANNOT_READS.out.resfinder_long
			resfinder_short     = AMR_ANNOT_READS.out.resfinder_short
			plasmidfinder_long  = AMR_ANNOT_READS.out.plasmidfinder_long
			plasmidfinder_short = AMR_ANNOT_READS.out.plasmidfinder_short
			
	    multireport_folder = Channel.empty() //MULTIREPORT.out.folder
    	multireport_html   = Channel.empty() //MULTIREPORT.out.html.map({m,x -> x})
    	multireport_xlsx   = Channel.empty() //MULTIREPORT.out.xlsx.map({m,x -> x})
}














