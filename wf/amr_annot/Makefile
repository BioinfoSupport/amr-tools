


all:nextflow_schema.json

nextflow_schema.json:assets/schema/nextflow_schema_template.json
	jq \
	  --slurpfile assemblies assets/schema/params/assemblies.json \
	  '(.properties.assemblies=$$assemblies[0])' \
	  assets/schema/nextflow_schema_template.json \
	> $@

clean:
	rm -f nextflow_schema.json
