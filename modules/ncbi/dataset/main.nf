

process NCBI_DATASET_DOWNLOAD_GENOME {
    container 'docker.io/staphb/ncbi-datasets:18.18.0'
    memory '8 GB'
    cpus 1
    time '30 min'
    input:
    	val args
    output:
  		path('dataset')
    script:
	    """
	    datasets download genome ${args} --filename dataset.zip \
	    	&& unzip dataset.zip -d dataset
	    """
}