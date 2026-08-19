


class Utils {

	static void logFailed(t) {
		if (t.workDir) {
    	def errFile = t.workDir.resolve('.command.err')
	    System.err.println """=== FAILED: ${t.name} (exit ${t.exitStatus}) ===
	      Work dir: ${t.workDir}
	      ${errFile.exists() && errFile.size() > 0 ? "  Stderr:\n${errFile.text.takeRight(2000)}" : "  (no stderr)"}
	    """.stripIndent()
		}
	}
	
	static String errorStrategyRetryOnce(t) {
		if (t.attempt<2) return 'retry'
		logFailed(t)
		return 'ignore'
	}
	
}


