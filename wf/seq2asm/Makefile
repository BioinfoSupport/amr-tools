


all:nextflow_schema.json

nextflow_schema.json:assets/schema/nextflow_schema_template.json assets/schema/params/readsets.json assets/schema/params/assembler.json
	jq \
	  --slurpfile readsets assets/schema/params/readsets.json \
    --slurpfile assembler assets/schema/params/assembler.json \
    '(.properties.assembler=$$assembler[0]) | (.properties.readsets=$$readsets[0])' \
	  $< \
	> $@

clean:
	rm -f nextflow_schema.json
