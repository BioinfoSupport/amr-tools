
import nextflow.Channel

class Assembler {
	static Map fromParams(opts) {
		return [name:opts.name?:"long_flye_medaka",args:opts.args?:[:]]
	}
}

