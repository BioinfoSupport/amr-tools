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



workflow ANNOTATE_ASSEMBLY_FOR_ORG {
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



workflow ANNOTATE_ASSEMBLY {
	take:
		opts
		asm_ch
		lr_ch
		sr_ch
	main:
		// ------------------------------------------------------------------			
		// Convert input files formats
		// ------------------------------------------------------------------
		// Uncompress .fasta.gz files when needed
		asm_ch = asm_ch.branch({meta,f -> 
			gz: f.name =~ /\.gz$/
			fa: true
		})
		asm_ch = asm_ch.fa.mix(GZIP_DECOMPRESS_FASTA(asm_ch.gz))

		// -------------------
		// Reads processing
		// -------------------
		PLASMIDFINDER_LONGREAD(lr_ch.filter({!opts.skip_plasmidfinder_longread}))
		RESFINDER_LONGREAD(lr_ch.filter({!opts.skip_resfinder_longread}),'nanopore')
		PLASMIDFINDER_SHORTREAD(sr_ch.filter({!opts.skip_plasmidfinder_shortread}))
		RESFINDER_SHORTREAD(lr_ch.filter({!opts.skip_resfinder_shortread}),'illumina')
			
		// ---------------------------------------------------------------------
		// Tools that can run directly on a FASTA without specifying an organism
		// ---------------------------------------------------------------------
		fai_ch = SAMTOOLS_FAIDX(asm_ch)

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
		ann_ch = ANNOTATE_ASSEMBLY_FOR_ORG(opts,asm_ch,orgfinder_ch.org_name)


		// Results aggregation and reporting
		MULTIREPORT(
			asm_ch,
			fai_ch,
			ann_ch.org_name,
			orgfinder_ch.orgfinder,
			amrfinderplus_ch,
			resfinder_ch,
			mobtyper_ch,
			plasmidfinder_ch,
			ann_ch.cgemlst,
			ann_ch.MLST,
			ann_ch.prokka,
			RESFINDER_LONGREAD.out,
			PLASMIDFINDER_LONGREAD.out,
			RESFINDER_SHORTREAD.out,
			PLASMIDFINDER_SHORTREAD.out
		)
		
	emit:
			orgfinder           = orgfinder_ch.orgfinder
			amrfinderplus       = amrfinderplus_ch
			resfinder           = resfinder_ch
			mobtyper            = mobtyper_ch
			plasmidfinder       = plasmidfinder_ch
			cgemlst             = ann_ch.cgemlst
			MLST                = ann_ch.MLST
			prokka              = ann_ch.prokka
			resfinder_long      = RESFINDER_LONGREAD.out
			resfinder_short     = RESFINDER_SHORTREAD.out
			plasmidfinder_long  = PLASMIDFINDER_LONGREAD.out
			plasmidfinder_short = PLASMIDFINDER_SHORTREAD.out
			
	    multireport_folder = MULTIREPORT.out.folder
    	html_report        = MULTIREPORT.out.html
    	xlsx_report        = MULTIREPORT.out.xlsx
}
