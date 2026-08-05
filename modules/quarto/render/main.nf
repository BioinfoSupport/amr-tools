
process QUARTO_RENDER {
	  container "registry.gitlab.unige.ch/amr-genomics/rscript:v2"
    memory '8 GB'
    cpus 2
    time '30 min'
    input:
    		tuple(val(meta),path(files),val(params))
    		each path('report_template.qmd')
    		each path(extra)
    output:
        tuple(val(meta),path('report_template.html'),emit:'html')
    script:
				"""
				HOME=/tmp quarto render ${task.ext.args?:''} 'report_template.qmd' ${params}
				"""
		stub:
				"""
				touch report_template.html
				"""
}
