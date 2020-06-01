# augur-nf
This repository contains a Nextflow implementation of Nexstrain's Augur pipeline (https://github.com/nextstrain/augur). From Augur's github page: 

"Augur is the bioinformatics toolkit to track evolution from sequence and serological data. It provides a collection of commands which are designed to be composable into larger processing pipelines."

The processes in this pipeline depend on a [Docker-ized version of Augur](https://github.com/StaPH-B/docker-builds/tree/master/augur/7.0.2).

#### Using the pipeline
The pipeline is designed to start from assembled genomes in fasta format. All assemblies must be in the same fasta file. Additionally, a reference genome in Genbank format, a tsv of metadata, a tsv containing a hex color code for each country of origin in your sample, a tsv containing latitude and longitude coordinates for each country of origin in your sample and an Auspice config file are required. Please see this [tutorial](https://nextstrain-augur.readthedocs.io/en/stable/tutorials/zika_tutorial.html) for more information on the file format requirements for each input.

Start the pipeline using `nextflow augur.nf` and follow the options for running the pipeline.

```
usage: nextflow augur.nf [--output] [--filter] [--min_date] [--traits]
             --sequences <path> --refernce <path> --metadata <path> --colors <path> 
             --lat_long <path> --config <path>

required arguments:
  sequences                Path to fasta file of assemblies
  reference                Path to reference genome in Genbank format
  metadata                 Path to metadata file
  colors                   Path to tsv file containing hex color codes and country of origin
  lat_long                 Path to tsv file containing longitude and latitude coordinates of countries of origin
  config                   Path to Auspice config file
  
optional arguments:
  --output  <output_path>  Path to ouput directory, default "augur_results"
  --threads                Number of threads to use during alingment, default 1
  --filter  <path_to_file> Filter out a list of sequences
  --group_by               What metadata column to group sequences by, default "country year month" (use with --filter)
  --seq_per_group          Number of sequences per group (use with --filter and --seq_per_group)
  --traits                 Perform ancestral reconstruction of a trait, default true
  --columns                What metadata column(s) to use for ancestral reconstruction, default "region country" (use with --trait)
  --<process_name>_addn    Any additional parameters for each Augur process
```
