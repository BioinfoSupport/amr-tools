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
include { PROKKA_RUN                               } from './modules/tseemann/prokka'
include { MLST                                     } from './modules/tseemann/mlst'
include { MLST as CGEMLST                          } from './modules/cgetools/mlst'
include { MOBTYPER_RUN                             } from './modules/mobsuite/mobtyper.nf'
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
	if (opts.containsKey(key)) return opts[key]
  if (meta.containsKey(key)) return meta[key]
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
			.filter({!opts.skip_cgemlst})
			.map({meta,fa,org_name -> [meta, fa, tool_args('cgemlst',meta,opts,org_name)]})
			.filter({meta,fasta,args -> args!=null})
			| CGEMLST
		MLST_ch = fa_org_ch
			.filter({!opts.skip_MLST})
			.map({meta,fa,org_name -> [meta, fa, tool_args('MLST',meta,opts,org_name)]})
			.filter({meta,fasta,args -> args!=null})
			| MLST

		// PROKKA annotations
		prokka_ch = fa_org_ch
			.filter({!opts.skip_prokka})
			.map({meta,fa,org_name -> [meta, fa, tool_args('prokka',meta,opts,org_name)]})
			.filter({meta,fasta,args -> args!=null})
			| PROKKA_RUN

	emit:
		cgemlst = cgemlst_ch
		MLST = MLST_ch
		prokka = prokka_ch
		org_name = fa_org_ch.map({meta,fa,org_name -> [meta,org_name]})
}



workflow AMR_ANNOT_ASSEMBLY {
	take:
		opts
		asm_ch
	main:
		// CGE - RESFINDER
		resfinder_ch = RESFINDER(asm_ch.filter({!opts.skip_resfinder}),'fasta')

		// Plasmid typing
		plasmidfinder_ch = asm_ch.filter({!opts.skip_plasmidfinder}) | PLASMIDFINDER
		
		// NCBI AMRfinder+
		if (opts.skip_amrfinderplus) {
			amrfinderplus_ch = Channel.empty()
		} else {
			amrfinderplus_db = AMRFINDERPLUS_UPDATE()
			amrfinderplus_ch = AMRFINDERPLUS_RUN(
					asm_ch
						.map({meta,fasta -> [meta,fasta,tool_args('amrfinderplus',meta,opts)]})
						.filter({meta,fasta,args -> args!=null}),
					amrfinderplus_db
			)
		}
		
		// MOBsuite - MOBtyper
		mobtyper_ch = asm_ch
			.filter({!opts.skip_mobtyper})
			.map({meta,fasta -> [meta,fasta,tool_args('mobtyper',meta,opts)]})
			.filter({meta,fasta,args -> args!=null})
      | MOBTYPER_RUN

    // Run orgfinder to auto detect organism
		orgfinder_ch = asm_ch.filter({!opts.skip_orgfinder}) | ORGFINDER_DETECT
			
		// Speciator
		speciator_ch = asm_ch.filter({!opts.skip_speciator}) | SPECIATOR

			
		// ---------------------------------------------------------------------
		// Organism specific tools
		// ---------------------------------------------------------------------
		ann_ch = AMR_ANNOT_ASSEMBLY_FOR_ORG(opts,asm_ch,orgfinder_ch.org_name)

	emit:
			orgfinder           = orgfinder_ch.orgfinder
			amrfinderplus       = amrfinderplus_ch
			resfinder           = resfinder_ch
			mobtyper            = mobtyper_ch
			plasmidfinder       = plasmidfinder_ch
			org_name            = ann_ch.org_name			
			cgemlst             = ann_ch.cgemlst
			MLST                = ann_ch.MLST
			prokka              = ann_ch.prokka
}



workflow AMR_ANNOT_READS {
	take:
		opts
		fqs_ch
		fql_ch
	main:
		PLASMIDFINDER_LONGREAD(fql_ch.filter({!opts.skip_plasmidfinder_longread}))
		RESFINDER_LONGREAD(fql_ch.filter({!opts.skip_resfinder_longread}),'nanopore')
		PLASMIDFINDER_SHORTREAD(fqs_ch.filter({!opts.skip_plasmidfinder_shortread}))
		RESFINDER_SHORTREAD(fql_ch.filter({!opts.skip_resfinder_shortread}),'illumina')
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
			
	    multireport_folder = MULTIREPORT.out.folder
    	multireport_html   = MULTIREPORT.out.html.map({m,x -> x})
    	multireport_xlsx   = MULTIREPORT.out.xlsx.map({m,x -> x})
}














