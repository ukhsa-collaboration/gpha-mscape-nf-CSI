# CSI: Clinical Sample Investigation

Pipeline for reporting investigations into clinical metagenomic samples. Currently
the pipeline only includes a feature to analyse unclassified ONT reads.

## Installation

Currently as a workaround, a module requires the installation of a conda environment.

Clone repo and create environment:

`git clone git@github.com:ukhsa-collaboration/gpha-mscape-nf-CSI.git`

`conda env create –f assets/environment.yaml`

`conda activate taxaplease`

## Inputs

To run the pipeline you require a samplesheet specifying the path to the input fastq (already
downloaded from S3). An example can be found here: assets/samplesheet_example.csv


Parameters can be specified on the CLI or within nextflow.config

| Parameter | Required | Description |
| -------- | ------- | ------- |
| --input | Yes | Path to samplesheet |
| --outdir | No | Path to output dir  |
| --ref_fasta | Yes | Path to references used for mapping. For unclassified samples used the tobacco genome |
| --collapse_rank | No | Taxonomic rank used to group Kraken2 results. Defaults to family  |
| --pluspf | Yes | Path to PlusPF version in production  |
| --pluspf_recent | Yes | Path to most recent PlusPF database |
| --viper | Yes | Path to viral database |

## Usage


```
nextflow run main.nf --input samplesheet.csv --ref_fasta assets.tobacco_refs.fasta.gz --pluspf path/to/pluspfdb --pluspf_recent path/to/recent/pluspfdb --viper path/to/viraldb
```

## Outputs

| Output Dir | Output | Description |
| -------- | ------- | ------- |
| Mapping | id_mapping_stats.csv | Mapping information |
| Mapping | id_mapping_stats_text.txt | Descriptive text for reporting |
| Flye | assembly_graph.gfa | Flye output |
| Flye | assembly_info.txt | Flye output |
| Flye | assembly.fasta | Flye assembled contigs |
| Flye | id_de_novo_assembly_stats_text.txt | Descriptive text for reporting |
| Flye | id_read2contig_mapping_stats.csv | Reads mapped to contigs |
| Flye | id_de_novo_assembly_stats.html | Contig summary plot |
| Kraken | id.sampletype.database.kraken_out.txt | Kraken2 stdout |
| Kraken | id_sampletype_kraken2_database_rank.csv | Summary of Kraken2 classifications. 1 per contig/read per database |
| Sankey | id_sampletype_sankey.html | Kraken2 classification plot 1 per read/contig |
| Sankey | id_sampletype_sankey_text.txt | Descriptive text for reporting |

## Downloading reference genome

You can download reference genomes using the NCBI datasets command line tool. Installation instructions can be found here: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/

To download the tobacco genome use the following commands:
```
datasets download genome taxon 4097 --reference --filename tobacco_refs.zip
unzip -o tobacco_refs.zip
cat  ncbi_dataset/data/*/*.fna > tobacco_refs.fasta
gzip tobacco_refs.fasta.gz
```
