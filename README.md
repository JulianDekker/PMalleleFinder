# TR_diversity

## Table of contents
* [General info](#general-info)
* [Requirements](#Requirements)
* [Setup](#setup)
* [Command line options](#Command-line-options)

## General info
This pipeline will allow you to extract allele sequences and population information for *TR* genes.
	
## Requirements
* python=3.7.1
* biopython=1.76
* BCBio.gff=0.6.6
* pandas=1.0.1
* numpy=1.18.1
* R=3.6.1
* samtools=1.6
* tabix=1.6
	
## Setup
Clone the project structure and run `main_pipe.py`

Example usage:
```
$ python3 main_pipe.py --gff=gencode.v32.chr_patch_hapl_scaff.annotation.gff3 --vcf="ALL.chr7.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz;ALL.chr14.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz" --popfile=pop-info.tsv --ref=GRCh38_full_analysis_set_plus_decoy_hla.fa
```

## Command line options
* `--gff=gencodefile` **required**
This option required a .gff or .gff3 input file to read gene locations from
* `--vcf=vcffile.vcf` **required**
This option requires a .vcf.gz and .vcf.tbi file to extract variation
* `--popfile=populationfile` **required**
This option requires a tab seperated population file in the format of `Sample  pop Superpop`
* `--ref=reference.fasta` **required**
This option requires the reference genome fasta and an index for the fasta in the form of a .fai file. 
* `--filter=filterfile`
This option allows the specification of genes to filter on. Every new filter should be on a new line in this file. If not specified no filter will be used.
* `--threshold=integer`
this option allows the specification of the threshold used to discard sequences. If not specified the threshold is set to 4.
* `--rss=TRUE|FALSE` 
This option if set to TRUE will also generate sequences containing the RSS for *V,D,J* genes. 60bp sequence will be cut.