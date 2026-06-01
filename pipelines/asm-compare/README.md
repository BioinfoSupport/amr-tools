
# asm-compare

`asm-compare` a simple nextflow pipeline to compare two assemblies and generate a VCF.

The pipeline steps include:

  1) Map the two assemblies with `minimap2 -x asm5` and produce an indexed `bam` file

  2) Generate a VCF listing the mutations with `bcftools mpileup`
  
  3) Generate a simple 3 columns table list of mutations from the VCF


# Installation

The pipeline depends on nextflow, that can be installed with: 

```bash
curl -s https://get.nextflow.io | bash
```

# Usage

```bash
nextflow run BioinfoSupport/amr-tools/pipelines/asm-compare --ref=ref_assembly.fasta --query=query_assembly.fasta
```

# Test the pipeline

```bash
nextflow run BioinfoSupport/amr-tools/pipelines/asm-compare -profile standard,test
```

