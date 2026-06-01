
# asm-compare

`asm-compare` is a nextflow pipeline to process compare two assemblies and generate a VCF and an HTML report.

The pipeline steps include:

  1) Map the assemblies with `minimap2 -x asm5` and produce an indexed `bam` file

  2) Generate a VCF listing the mutation with `bcftools`
  
  3) 


# Installation

The pipeline depends on nextflow, that can be installed with: 

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run BioinfoSupport/amr-tools/pipelines/asm-compare --ref=ref_assembly.fasta --fasta=assembly.fasta
```

# Test the pipeline

```bash
nextflow run BioinfoSupport/amr-tools/pipelines/asm-compare -profile standard,test
```
