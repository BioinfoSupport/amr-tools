


all:nextflow_schema.json

nextflow_schema.json:assets/schema/params/readsets.json assets/schema/nextflow_schema_template.json
	jq \
	  --slurpfile readsets assets/schema/params/readsets.json \
	  --slurpfile assembly assets/schema/params/assembly.json \
	  '(.properties.readsets=$$readsets[0]) | (.properties.assembly=$$assembly[0])' \
	  assets/schema/nextflow_schema_template.json \
	> $@

clean:
	rm -f nextflow_schema.json
