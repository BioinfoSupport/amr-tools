


workflow {
	main:
		println('Hello world')

	publish:
		bam = Channel.empty()
		bai = Channel.empty()
		vcf = Channel.empty()
		html = Channel.empty()
}

output {

	bam {
		path { "/"}
	}
	bai {
		path { "/"}
	}
	vcf {
		path { "/"}
	}
	html {
		path { "reports/"}
	}
	
}




