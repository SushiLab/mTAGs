<p align="center">
<img src="img/mtags_logo.png" width="500"  />
</p>

# mTAGs: taxonomic profiling using degenerate consensus reference sequences of ribosomal RNA genes






  * [Installation](#installation)
    + [Installation using conda](#installation-using-conda)
    + [Manual installation](#manual-installation)
  * [Usage](#usage)
    + [PROFILE](#profile)
    + [MERGE](#merge)
  * [Output Files](#output-files)
    + [Profile file](#profile-file)
    + [Single insert annotation](#single-insert-annotation)
  * [Database](#database)
  * [References](#references)

mTAGs is a tool for the taxonomic profiling of metagenomes. It detects sequencing reads belonging to the small subunit of the ribosomal RNA (SSU-rRNA) gene and annotates them through the alignment to full-length degenerate consensus SSU-rRNA reference sequences. The tool is capable of processing single-end and pair-end metagenomic reads, takes advantage of the information contained in any region of the SSU-rRNA gene and provides relative abundance profiles at multiple taxonomic ranks (Domain, Phylum, Class, Order, Family, Genus and OTUs defined at a 97% sequence identity cutoff).


The tool is developed by Hans-Joachim Ruscheweyh and Guillem Salazar and distributed under the [![License GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html).


If you use mTAGs, please cite:


> Salazar G\*, Ruscheweyh H-J\*, Hildebrand F, Acinas S and Sunagawa S. **mTAGs: taxonomic profiling using degenerate
consensus reference sequences of ribosomal RNA gene.** *Bioinformatics*, 2021.

Analyses in the publication were executed using version 1.0.0

Questions/Comments? Write a github issue.

## Installation


`mTAGs` is written in `python` and has the following dependencies:


- Python>=v3.7
- [vsearch](https://github.com/torognes/vsearch/) (tested: v2.15.0)
- [hmmer](http://hmmer.org/) (tested: v3.3)



### Installation using conda

The easiest way to install mTAGs is to use the conda package manager, which will automatically create an environment with dependencies installed in the correct version.

```bash
$ conda create -n mtags python=3.7 hmmer vsearch
$ source activate mtags
# or
$ conda activate mtags

$ git clone https://github.com/SushiLab/mTAGs.git
$ cd mTAGs
$ pip install -r requirements.txt -e . 

# Download mTAGs database
$ mtags download

2021-06-21 12:02:40,294 INFO: Starting mTAGs
2021-06-21 12:02:40,294 INFO: Start downloading the mTAGs database. ~600MB
2021-06-21 12:05:17,883 INFO: Finished downloading the mTAGs database.
2021-06-21 12:05:17,883 INFO: Finishing mTAGs

$ mtags
```
<details><summary>mTAGs output</summary>
<p>


```bash
Program:    mTAGs - taxonomic profiling using degenerate consensus reference
            sequences of ribosomal RNA gene
Version:    1.0.0
Reference:  Salazar, Ruscheweyh, et al. mTAGs: taxonomic profiling using
            degenerate consensus reference sequences of ribosomal RNA
            gene. Bioinformatics (2021)

Usage: mtags <command> [options]
Command:

-- General
    profile     Extract and taxonomically annotate rRNA reads in metagenomic samples
    merge       Merge profiles

-- Expert
    extract     Extract rRNA reads in metagenomic samples
    annotate    Annotate and quantify rRNA reads

-- Installation
    download    Download the mTAGs database - Once after download of the tool

The database needs to be downloaded in the last step of the installation. This needs to be done once and before the first metagenomic samples can be processed:
```

</p>
</details>


### Manual installation

Manual installation is possible but not recommended. Install via pip after installation of dependencies:

```bash
$ git clone https://github.com/SushiLab/mTAGs.git
$ cd mTAGs
$ pip install -r requirements.txt -e . 

# Download mTAGs database
$ mtags download

2021-06-21 12:02:40,294 INFO: Starting mTAGs
2021-06-21 12:02:40,294 INFO: Start downloading the mTAGs database. ~600MB
2021-06-21 12:05:17,883 INFO: Finished downloading the mTAGs database.
2021-06-21 12:05:17,883 INFO: Finishing mTAGs

$ mtags
```
<details><summary>mTAGs output</summary>
<p>


```bash
Program:    mTAGs - taxonomic profiling using degenerate consensus reference
            sequences of ribosomal RNA gene
Version:    1.0.0
Reference:  Salazar, Ruscheweyh, et al. mTAGs: taxonomic profiling using
            degenerate consensus reference sequences of ribosomal RNA
            gene. Bioinformatics (2021)

Usage: mtags <command> [options]
Command:

-- General
    profile     Extract and taxonomically annotate rRNA reads in metagenomic samples
    merge       Merge profiles

-- Expert
    extract     Extract rRNA reads in metagenomic samples
    annotate    Annotate and quantify rRNA reads

-- Installation
    download    Download the mTAGs database - Once after download of the tool

The database needs to be downloaded in the last step of the installation. This needs to be done once and before the first metagenomic samples can be processed:
```

</p>
</details>


## Usage

The tool is split in to two steps: profiling and merging. The first step (`mtags profile [options]`) uses HMM models to extract potential rRNA sequences from metagenomic data and annotates them taxonomically through the alignment of these sequences against a modified Silva database. The second step (`mtags merge [options]`) is a function that merges taxonomic profiles from different metagenomic samples. The steps for extraction and annotation of rRNA sequences are grouped into a single command but can also be run independently (`mtags extract [options]` and `mtags annotate [options]`).



### PROFILE

This step uses precomputed HMM models to extract rRNA sequences from a metagenomic sample. The rRNA sequences are then aligned against a clustered rRNA database to annotate sequences and profile samples. `mTAGs` takes as input fasta/fastq files with quality controlled sequencing data.



```bash

$ mtags profile
Program:    mTAGs - taxonomic profiling using degenerate consensus reference
            sequences of ribosomal RNA gene
Version:    1.0.0
Reference:  Salazar, Ruscheweyh, et al. mTAGs: taxonomic profiling using
            degenerate consensus reference sequences of ribosomal RNA
            gene. Bioinformatics (2021)

Usage: mtags profile [options]

Input options:
    -f  FILE [FILE ...]   Forward reads file. Can be fasta/fastq and gzipped.
    -r  FILE [FILE ...]   Reverse reads file. Can be fasta/fastq and gzipped.
    -s  FILE [FILE ...]   Single/merge reads file. Can be fasta/fastq and gzipped.

Output options:
    -o  DIR               Output folder [Required]

Other options:
    -n  STR               Samplename [Required]
    -t  INT               Number of threads. [4]
    -ma INT               Maxaccepts, vsearch parameter. Larger
                          numbers increase sensitivity and runtime. [1000]
    -mr INT               Maxrejects, vsearch parameter. Larger
                          numbers increase sensitivity and runtime. [1000]

# Example usage of the mTAGs profile routine
$ mtags profile -f sample.1.fq.gz -r sample.2.fq.gz -s sample.s.fq.gz sample.m.fq.gz -o output -t 4 -n sample -ma 1000 -mr 1000
```
<details><summary>mTAGs log</summary>
<p>


```bash
2021-06-21 09:04:48,644 INFO: Starting mTAGs
2021-06-21 09:04:48,646 INFO: Extracting FastA and revcomp FastA from input/sample.1.fq.gz
2021-06-21 09:04:59,536 INFO: Processed reads:	824523
2021-06-21 09:04:59,536 INFO: Finished extracting. Found 824523 sequences.
2021-06-21 09:04:59,536 INFO: Start detecting rRNA sequences in FastA files
2021-06-21 09:04:59,536 INFO: Start detecting rRNA sequences for molecule=ssu
2021-06-21 09:04:59,536 INFO: Executing:	hmmsearch --cpu 4 -o sample/sample.1.fq.gz_fw.fasta_ssu.hmmer --domtblout sample/sample.1.fq.gz_fw.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm sample/sample.1.fq.gz_fw.fasta
2021-06-21 09:05:06,982 INFO: Finished hmmsearch
2021-06-21 09:05:06,988 INFO: Executing:	hmmsearch --cpu 4 -o sample/sample.1.fq.gz_rev.fasta_ssu.hmmer --domtblout sample/sample.1.fq.gz_rev.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm sample/sample.1.fq.gz_rev.fasta
2021-06-21 09:05:14,442 INFO: Finished hmmsearch
2021-06-21 09:05:14,449 INFO: Finished detecting rRNA sequences for molecule=ssu
2021-06-21 09:05:14,449 INFO: Start detecting rRNA sequences for molecule=lsu
2021-06-21 09:05:14,450 INFO: Executing:	hmmsearch --cpu 4 -o sample/sample.1.fq.gz_fw.fasta_lsu.hmmer --domtblout sample/sample.1.fq.gz_fw.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm sample/sample.1.fq.gz_fw.fasta
2021-06-21 09:05:35,255 INFO: Finished hmmsearch
2021-06-21 09:05:35,266 INFO: Executing:	hmmsearch --cpu 4 -o sample/sample.1.fq.gz_rev.fasta_lsu.hmmer --domtblout sample/sample.1.fq.gz_rev.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm sample/sample.1.fq.gz_rev.fasta
2021-06-21 09:05:55,845 INFO: Finished hmmsearch
2021-06-21 09:05:55,859 INFO: Finished detecting rRNA sequences for molecule=lsu
2021-06-21 09:05:55,859 INFO: Found 4143 potential rRNA sequences.
2021-06-21 09:05:55,859 INFO: Finished detecting rRNA sequences from FastA files.
2021-06-21 09:05:55,859 INFO: Finding best molecule for each read
2021-06-21 09:05:55,867 INFO: Finished finding best molecule for each read
2021-06-21 09:05:55,867 INFO: Start extracting reads/writing output
2021-06-21 09:05:58,911 INFO: Processed reads:	824523
2021-06-21 09:05:58,912 INFO: Finished extracting reads/writing output
2021-06-21 09:05:58,912 INFO: euk_lsu	1571
2021-06-21 09:05:58,912 INFO: bac_lsu	1114
2021-06-21 09:05:58,912 INFO: euk_ssu	865
2021-06-21 09:05:58,912 INFO: bac_ssu	583
2021-06-21 09:05:58,912 INFO: arc_lsu	9
2021-06-21 09:05:58,912 INFO: arc_ssu	1
2021-06-21 09:05:58,927 INFO: Extracting FastA and revcomp FastA from input/sample.2.fq.gz
2021-06-21 09:06:10,001 INFO: Processed reads:	824523
2021-06-21 09:06:10,002 INFO: Finished extracting. Found 824523 sequences.
2021-06-21 09:06:10,002 INFO: Start detecting rRNA sequences in FastA files
2021-06-21 09:06:10,002 INFO: Start detecting rRNA sequences for molecule=ssu
2021-06-21 09:06:10,002 INFO: Executing:	hmmsearch --cpu 4 -o sample/sample.2.fq.gz_fw.fasta_ssu.hmmer --domtblout sample/sample.2.fq.gz_fw.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm sample/sample.2.fq.gz_fw.fasta
2021-06-21 09:06:17,074 INFO: Finished hmmsearch
2021-06-21 09:06:17,084 INFO: Executing:	hmmsearch --cpu 4 -o sample/sample.2.fq.gz_rev.fasta_ssu.hmmer --domtblout sample/sample.2.fq.gz_rev.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm sample/sample.2.fq.gz_rev.fasta
2021-06-21 09:06:24,016 INFO: Finished hmmsearch
2021-06-21 09:06:24,022 INFO: Finished detecting rRNA sequences for molecule=ssu
2021-06-21 09:06:24,022 INFO: Start detecting rRNA sequences for molecule=lsu
2021-06-21 09:06:24,022 INFO: Executing:	hmmsearch --cpu 4 -o sample/sample.2.fq.gz_fw.fasta_lsu.hmmer --domtblout sample/sample.2.fq.gz_fw.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm sample/sample.2.fq.gz_fw.fasta
2021-06-21 09:06:44,641 INFO: Finished hmmsearch
2021-06-21 09:06:44,653 INFO: Executing:	hmmsearch --cpu 4 -o sample/sample.2.fq.gz_rev.fasta_lsu.hmmer --domtblout sample/sample.2.fq.gz_rev.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm sample/sample.2.fq.gz_rev.fasta
2021-06-21 09:07:04,995 INFO: Finished hmmsearch
2021-06-21 09:07:05,005 INFO: Finished detecting rRNA sequences for molecule=lsu
2021-06-21 09:07:05,005 INFO: Found 4144 potential rRNA sequences.
2021-06-21 09:07:05,005 INFO: Finished detecting rRNA sequences from FastA files.
2021-06-21 09:07:05,005 INFO: Finding best molecule for each read
2021-06-21 09:07:05,012 INFO: Finished finding best molecule for each read
2021-06-21 09:07:05,012 INFO: Start extracting reads/writing output
2021-06-21 09:07:08,038 INFO: Processed reads:	824523
2021-06-21 09:07:08,039 INFO: Finished extracting reads/writing output
2021-06-21 09:07:08,039 INFO: euk_lsu	1607
2021-06-21 09:07:08,039 INFO: bac_lsu	1121
2021-06-21 09:07:08,039 INFO: euk_ssu	833
2021-06-21 09:07:08,039 INFO: bac_ssu	577
2021-06-21 09:07:08,039 INFO: arc_lsu	5
2021-06-21 09:07:08,039 INFO: arc_ssu	1
2021-06-21 09:07:08,055 INFO: Extracting FastA and revcomp FastA from input/sample.s.fq.gz
2021-06-21 09:07:14,410 INFO: Processed reads:	472552
2021-06-21 09:07:14,410 INFO: Finished extracting. Found 472552 sequences.
2021-06-21 09:07:14,410 INFO: Start detecting rRNA sequences in FastA files
2021-06-21 09:07:14,410 INFO: Start detecting rRNA sequences for molecule=ssu
2021-06-21 09:07:14,411 INFO: Executing:	hmmsearch --cpu 4 -o sample/sample.s.fq.gz_fw.fasta_ssu.hmmer --domtblout sample/sample.s.fq.gz_fw.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm sample/sample.s.fq.gz_fw.fasta
2021-06-21 09:07:18,525 INFO: Finished hmmsearch
2021-06-21 09:07:18,529 INFO: Executing:	hmmsearch --cpu 4 -o sample/sample.s.fq.gz_rev.fasta_ssu.hmmer --domtblout sample/sample.s.fq.gz_rev.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm sample/sample.s.fq.gz_rev.fasta
2021-06-21 09:07:22,690 INFO: Finished hmmsearch
2021-06-21 09:07:22,696 INFO: Finished detecting rRNA sequences for molecule=ssu
2021-06-21 09:07:22,697 INFO: Start detecting rRNA sequences for molecule=lsu
2021-06-21 09:07:22,697 INFO: Executing:	hmmsearch --cpu 4 -o sample/sample.s.fq.gz_fw.fasta_lsu.hmmer --domtblout sample/sample.s.fq.gz_fw.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm sample/sample.s.fq.gz_fw.fasta
2021-06-21 09:07:34,547 INFO: Finished hmmsearch
2021-06-21 09:07:34,554 INFO: Executing:	hmmsearch --cpu 4 -o sample/sample.s.fq.gz_rev.fasta_lsu.hmmer --domtblout sample/sample.s.fq.gz_rev.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm sample/sample.s.fq.gz_rev.fasta
2021-06-21 09:07:46,666 INFO: Finished hmmsearch
2021-06-21 09:07:46,675 INFO: Finished detecting rRNA sequences for molecule=lsu
2021-06-21 09:07:46,675 INFO: Found 2236 potential rRNA sequences.
2021-06-21 09:07:46,676 INFO: Finished detecting rRNA sequences from FastA files.
2021-06-21 09:07:46,676 INFO: Finding best molecule for each read
2021-06-21 09:07:46,679 INFO: Finished finding best molecule for each read
2021-06-21 09:07:46,680 INFO: Start extracting reads/writing output
2021-06-21 09:07:48,487 INFO: Processed reads:	472552
2021-06-21 09:07:48,488 INFO: Finished extracting reads/writing output
2021-06-21 09:07:48,488 INFO: euk_lsu	909
2021-06-21 09:07:48,489 INFO: bac_lsu	602
2021-06-21 09:07:48,489 INFO: euk_ssu	406
2021-06-21 09:07:48,489 INFO: bac_ssu	313
2021-06-21 09:07:48,489 INFO: arc_lsu	6
2021-06-21 09:07:48,498 INFO: Extracting FastA and revcomp FastA from input/sample.m.fq.gz
2021-06-21 09:08:04,216 INFO: Processed reads:	1000000
2021-06-21 09:08:16,996 INFO: Processed reads:	1888111
2021-06-21 09:08:16,997 INFO: Finished extracting. Found 1888111 sequences.
2021-06-21 09:08:16,997 INFO: Start detecting rRNA sequences in FastA files
2021-06-21 09:08:16,997 INFO: Start detecting rRNA sequences for molecule=ssu
2021-06-21 09:08:16,997 INFO: Executing:	hmmsearch --cpu 4 -o sample/sample.m.fq.gz_fw.fasta_ssu.hmmer --domtblout sample/sample.m.fq.gz_fw.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm sample/sample.m.fq.gz_fw.fasta
2021-06-21 09:08:39,040 INFO: Finished hmmsearch
2021-06-21 09:08:39,051 INFO: Executing:	hmmsearch --cpu 4 -o sample/sample.m.fq.gz_rev.fasta_ssu.hmmer --domtblout sample/sample.m.fq.gz_rev.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm sample/sample.m.fq.gz_rev.fasta
2021-06-21 09:09:00,018 INFO: Finished hmmsearch
2021-06-21 09:09:00,036 INFO: Finished detecting rRNA sequences for molecule=ssu
2021-06-21 09:09:00,036 INFO: Start detecting rRNA sequences for molecule=lsu
2021-06-21 09:09:00,037 INFO: Executing:	hmmsearch --cpu 4 -o sample/sample.m.fq.gz_fw.fasta_lsu.hmmer --domtblout sample/sample.m.fq.gz_fw.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm sample/sample.m.fq.gz_fw.fasta
2021-06-21 09:10:08,160 INFO: Finished hmmsearch
2021-06-21 09:10:08,184 INFO: Executing:	hmmsearch --cpu 4 -o sample/sample.m.fq.gz_rev.fasta_lsu.hmmer --domtblout sample/sample.m.fq.gz_rev.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm sample/sample.m.fq.gz_rev.fasta
2021-06-21 09:11:12,726 INFO: Finished hmmsearch
2021-06-21 09:11:12,762 INFO: Finished detecting rRNA sequences for molecule=lsu
2021-06-21 09:11:12,762 INFO: Found 10237 potential rRNA sequences.
2021-06-21 09:11:12,762 INFO: Finished detecting rRNA sequences from FastA files.
2021-06-21 09:11:12,762 INFO: Finding best molecule for each read
2021-06-21 09:11:12,784 INFO: Finished finding best molecule for each read
2021-06-21 09:11:12,784 INFO: Start extracting reads/writing output
2021-06-21 09:11:16,805 INFO: Processed reads:	1000000
2021-06-21 09:11:20,280 INFO: Processed reads:	1888111
2021-06-21 09:11:20,282 INFO: Finished extracting reads/writing output
2021-06-21 09:11:20,282 INFO: euk_lsu	3713
2021-06-21 09:11:20,282 INFO: bac_lsu	3017
2021-06-21 09:11:20,282 INFO: euk_ssu	1902
2021-06-21 09:11:20,282 INFO: bac_ssu	1582
2021-06-21 09:11:20,282 INFO: arc_lsu	16
2021-06-21 09:11:20,282 INFO: arc_ssu	7
False
2021-06-21 09:11:20,334 INFO: Input files:
2021-06-21 09:11:20,334 INFO: 	Forward readfiles:
2021-06-21 09:11:20,334 INFO: 		sample/sample.1.fq.gz_bac_ssu.fasta
2021-06-21 09:11:20,334 INFO: 		sample/sample.1.fq.gz_euk_ssu.fasta
2021-06-21 09:11:20,334 INFO: 		sample/sample.1.fq.gz_arc_ssu.fasta
2021-06-21 09:11:20,334 INFO: 	Reverse readfiles:
2021-06-21 09:11:20,334 INFO: 		sample/sample.2.fq.gz_euk_ssu.fasta
2021-06-21 09:11:20,334 INFO: 		sample/sample.2.fq.gz_bac_ssu.fasta
2021-06-21 09:11:20,334 INFO: 		sample/sample.2.fq.gz_arc_ssu.fasta
2021-06-21 09:11:20,334 INFO: 	Single/Merged readfiles:
2021-06-21 09:11:20,334 INFO: 		sample/sample.s.fq.gz_bac_ssu.fasta
2021-06-21 09:11:20,334 INFO: 		sample/sample.s.fq.gz_euk_ssu.fasta
2021-06-21 09:11:20,334 INFO: 		sample/sample.m.fq.gz_euk_ssu.fasta
2021-06-21 09:11:20,334 INFO: 		sample/sample.m.fq.gz_bac_ssu.fasta
2021-06-21 09:11:20,334 INFO: 		sample/sample.m.fq.gz_arc_ssu.fasta
2021-06-21 09:11:20,334 INFO: Output pattern:
2021-06-21 09:11:20,334 INFO: 	sample/sample
2021-06-21 09:11:20,334 INFO: Reading and testing FastA files
2021-06-21 09:11:20,355 INFO: Paired files stats:
2021-06-21 09:11:20,355 INFO: 	Found 1682 inserts in paired files of which 1178 are actually paired and 504 are singletons
2021-06-21 09:11:20,355 INFO: Singleton/Merged files stats:
2021-06-21 09:11:20,355 INFO: 	Found 4210 inserts in singleton files.
2021-06-21 09:11:20,355 INFO: Total stats:
2021-06-21 09:11:20,355 INFO: 	Paired inserts:	1178
2021-06-21 09:11:20,355 INFO: 	Singleton inserts:	4714
2021-06-21 09:11:20,355 INFO: Finished reading and testing FastA files
2021-06-21 09:11:20,356 INFO: Writing paired sequences to FastA:	sample/sample.pe.fasta
2021-06-21 09:11:20,358 INFO: Writing finished
2021-06-21 09:11:20,358 INFO: Writing singleton sequences to FastA:	sample/sample.se.fasta
2021-06-21 09:11:20,363 INFO: Writing finished
2021-06-21 09:11:20,363 INFO: Aligning paired sequence file:	vsearch --usearch_global sample/sample.pe.fasta --db SILVA-138_NR-97_complink_cons.vsearch.udb --id 0.97 --maxaccepts 1000 --maxrejects 1000 --strand both --userout sample/sample.pe.m8 --userfields query+target+qcov+id+qilo+qihi+tilo+tihi+alnlen+ids+mism+gaps+qstrand --query_cov 0.5 --threads 4 --output_no_hits
vsearch v2.15.0_macos_x86_64, 16.0GB RAM, 4 cores
https://github.com/torognes/vsearch

Reading UDB file SILVA-138_NR-97_complink_cons.vsearch.udb 100%
Reorganizing data in memory 100%
Creating bitmaps 100%
Parsing abundances 100%
341368570 nt in 233519 seqs, min 900, max 3718, avg 1462
Searching 100%
Matching unique query sequences: 1890 of 2356 (80.22%)
2021-06-21 09:13:26,382 INFO: Finished alignment
2021-06-21 09:13:26,383 INFO: Aligning singleton sequence file:	vsearch --usearch_global sample/sample.se.fasta --db SILVA-138_NR-97_complink_cons.vsearch.udb --id 0.97 --maxaccepts 1000 --maxrejects 1000 --strand both --userout sample/sample.se.m8 --userfields query+target+qcov+id+qilo+qihi+tilo+tihi+alnlen+ids+mism+gaps+qstrand --query_cov 0.5 --threads 4 --output_no_hits
vsearch v2.15.0_macos_x86_64, 16.0GB RAM, 4 cores
https://github.com/torognes/vsearch

Reading UDB file SILVA-138_NR-97_complink_cons.vsearch.udb 100%
Reorganizing data in memory 100%
Creating bitmaps 100%
Parsing abundances 100%
341368570 nt in 233519 seqs, min 900, max 3718, avg 1462
Searching 100%
Matching unique query sequences: 3644 of 4714 (77.30%)
2021-06-21 09:19:54,957 INFO: Finished alignment
2021-06-21 09:19:54,958 INFO: Parsing alignment files
2021-06-21 09:20:09,033 INFO: Finished parsing alignment files
2021-06-21 09:20:09,033 INFO: Filtering alignments by score and length
2021-06-21 09:20:10,883 INFO: Finished filtering
2021-06-21 09:20:10,883 INFO: Performing LCA calculations
2021-06-21 09:20:14,174 INFO: Finished LCA calculations
2021-06-21 09:20:14,174 INFO: Start writing output files
2021-06-21 09:20:14,174 INFO: Start writing:	sample/sample.bins
2021-06-21 09:20:14,631 INFO: Start writing:	sample/sample.otu.tsv
2021-06-21 09:20:14,759 INFO: Start writing:	sample/sample.genus.tsv
2021-06-21 09:20:14,763 INFO: Start writing:	sample/sample.family.tsv
2021-06-21 09:20:14,765 INFO: Start writing:	sample/sample.order.tsv
2021-06-21 09:20:14,766 INFO: Start writing:	sample/sample.class.tsv
2021-06-21 09:20:14,766 INFO: Start writing:	sample/sample.phylum.tsv
2021-06-21 09:20:14,767 INFO: Start writing:	sample/sample.domain.tsv
2021-06-21 09:20:14,767 INFO: Start writing:	sample/sample.root.tsv
2021-06-21 09:20:14,768 INFO: Finished writing output files
2021-06-21 09:20:15,178 INFO: Finished annotate command
2021-06-21 09:20:15,178 INFO: mTAGs finished successfully
```

</p>
</details>

 


### MERGE

The previous step produces taxonomic profiles for different samples. It is useful for downstream analysis to have multiple profiles in one file. The merge function combines multiple profiles into a tab-separated file with each column representating a sample:

```bash
$  mtags merge

2021-06-21 10:01:33,628 INFO: Starting mTAGs
Program:    mTAGs - taxonomic profiling using degenerate consensus reference
            sequences of ribosomal RNA gene
Version:    1.0.0
Reference:  Salazar, Ruscheweyh, et al. mTAGs: taxonomic profiling using
            degenerate consensus reference sequences of ribosomal RNA
            gene. Bioinformatics (2021)

Usage: mtags merge [options]

Input options:
    -i  FILE [FILE ...]   List of mTAGs bin files [Required]

Output options:
    -o  STR               Output prefix [Required]
                      
$ mtags merge -i *bins -o merged_profile

# Example usage of mTAGs merge:

$ mtags merge -i mtags2/*bins -o out
2021-06-21  16:36:06,642 INFO: Starting mTAGs
2021-06-21  16:36:06,642 INFO: Executing command merge
2021-06-21  16:36:06,676 INFO: Finished reading sample1.mTAG.bins. Found 35405 inserts.
2021-06-21  16:36:06,691 INFO: Finished reading sample2.mTAG.bins. Found 16265 inserts.
2021-06-21  16:36:06,718 INFO: Finished reading sample3.mTAG.bins. Found 25680 inserts.
2021-06-21  16:36:06,749 INFO: Finished reading sample4.mTAG.bins. Found 34954 inserts.

```



## Output Files

The tool produces two types of output, the taxonomic profiles (`*.tsv`) and the annotation of every insert (`*.bins`).

### Profile file

The tool produces profiles for 8 taxonomic ranks:

- otu
- genus
- family
- order
- class
- phylum
- domain
- root

All levels but the otu and the root levels are standard Silva level names. The root level was added to group all domains. The otu level entries represent taxonomic units that were generated by complete linkage clustering of distances between all sequences inside the respective genus.

The otu level is comprised of a complete taxonomic path plus the name of the otu (e.g. `otu__silva_138_complink_cons_otu_215825`) and the number of inserts that map to this otu. The output of a file looks the following:

```
root__Root;domain__Archaea;phylum__Euryarchaeota;class__Methanobacteria;order__Methanobacteriales;family__Methanobacteriaceae;genus__Methanobrevibacter;otu__silva_138_complink_cons_otu_128716	1
root__Root;domain__Bacteria;phylum__Actinobacteriota;class__Actinobacteria;order__Corynebacteriales;family__Mycobacteriaceae;genus__Mycobacterium;otu__silva_138_complink_cons_otu_99206	1
root__Root;domain__Bacteria;phylum__Proteobacteria;class__Alphaproteobacteria;order__Rickettsiales;family__Mitochondria;genus__unknown;otu__silva_138_complink_cons_otu_102401	1
root__Root;domain__Bacteria;phylum__Proteobacteria;class__Alphaproteobacteria;order__Rickettsiales;family__Mitochondria;genus__unknown;otu__silva_138_complink_cons_otu_133706	1
root__Root;domain__Bacteria;phylum__Proteobacteria;class__Gammaproteobacteria;order__Burkholderiales;family__Comamonadaceae;genus__Comamonas;otu__silva_138_complink_cons_otu_29404	1
root__Root;domain__Eukaryota;phylum__Basidiomycota;class__Agaricomycetes;order__unknown;family__unknown;genus__unknown;otu__silva_138_complink_cons_otu_187209	4
root__Root;domain__Eukaryota;phylum__Basidiomycota;class__Agaricomycetes;order__unknown;family__unknown;genus__unknown;otu__silva_138_complink_cons_otu_215825	16
Unaligned	427
Unassigned	738
```

The two last rows represent all sequences that could not be aligned or that could be aligned but could not be assigned at this taxonomic level. For almost all downstream analyses, it is better to remove this value, since it does not represent a single taxa but it needs to be taken into account when computing relative abundances.

The output for the phylum level then looks the following:

```
root__Root;domain__Archaea;phylum__Euryarchaeota	1
root__Root;domain__Bacteria;phylum__Actinobacteriota	2
root__Root;domain__Bacteria;phylum__Firmicutes	1
root__Root;domain__Bacteria;phylum__Proteobacteria	5
root__Root;domain__Eukaryota;phylum__Apicomplexa	1
root__Root;domain__Eukaryota;phylum__Ascomycota	57
root__Root;domain__Eukaryota;phylum__Basidiomycota	1310
root__Root;domain__Eukaryota;phylum__Phragmoplastophyta	1
root__Root;domain__Eukaryota;phylum__uncultured	1
root__Root;domain__Eukaryota;phylum__unknown	6
Unaligned	427
Unassigned	47
```

The file with merged profiles has one column per sample.

### Single insert annotation

The `*.bins` file with the taxonomic annotation of each insert looks the following:

`INSERTNAME	RANK	PATH`

```
HISEQ:384:HCVGNBCXY:2:2212:4002:35447   domain  root__Root;domain__Unaligned
HISEQ:384:HCVGNBCXY:1:1102:14918:19710  order   root__Root;domain__Eukaryota;phylum__Basidiomycota;class__Agaricomycetes;order__Polyporales
HISEQ:384:HCVGNBCXY:1:1102:9303:55752   domain  root__Root;domain__Unaligned
HISEQ:384:HCVGNBCXY:1:1102:17349:43609  class   root__Root;domain__Eukaryota;phylum__Basidiomycota;class__Agaricomycetes
HISEQ:384:HCVGNBCXY:1:1102:14716:53354  phylum  root__Root;domain__Eukaryota;phylum__Basidiomycota
HISEQ:384:HCVGNBCXY:1:1103:7884:25295   otu     root__Root;domain__Eukaryota;phylum__Basidiomycota;class__Agaricomycetes;order__Polyporales;family__unknown;genus__Postia;otu__silva_138_complink_cons_otu_15163
HISEQ:384:HCVGNBCXY:1:1103:19448:91267  otu     root__Root;domain__Eukaryota;phylum__Basidiomycota;class__Agaricomycetes;order__Agaricales;family__Tricholomataceae;genus__Leucopaxillus;otu__silva_138_complink_cons_otu_14867
HISEQ:384:HCVGNBCXY:1:1103:2033:80107   domain  root__Root;domain__Unaligned
HISEQ:384:HCVGNBCXY:1:1103:14200:48593  family  root__Root;domain__Eukaryota;phylum__Ascomycota;class__Sordariomycetes;order__Xylariales;family__unknown
```

## Database

The default database of mTAGs is a modified version of the Silva release 138. We cluster the database at genus level using a complete linkage algorithm on pairwise distances and then build one consensus sequence per cluster. We also generated a database from Silva release 128 using the same strategy for reproducibility reasons. However, **we encourage users to perform their analysis with most recent database version**, which is the one installed by default through the `mtags download` command.

If users need to use the Silva release 128 (NOT RECOMMENDED) they can do it by using the following commands:

```bash
# Use Silva release 128

$ cd mTAGs

# move default database 
$ mv db db_r138
# create database folder for release 128
$ mkdir db
$ cd db
$ wget https://sunagawalab.ethz.ch/share/MTAGS_DB/silva/s128/SILVA-128_NR-97_complink_cons.vsearch.udb.gz
$ wget https://sunagawalab.ethz.ch/share/MTAGS_DB/silva/s128/SILVA-128_NR-97_complink_cons.taxmap.gz
$ gunzip *gz
$ touch download.done
```
Using mTAGs will now use release version 128 of the database.


## References


Silva:

**Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO** (2013) *The SILVA ribosomal RNA gene database project: improved data processing and web-based tools.* Nucl. Acids Res. 41 (D1): D590-D596.

Vsearch:

**Rognes T, Flouri T, Nichols B, Quince C, Mahé F.** (2016) *VSEARCH: a versatile open source tool for metagenomics.* PeerJ 4:e2584. doi: 10.7717/peerj.2584


Hmmer:

**S. R. Eddy.** *Accelerated profile HMM searches.* PLOS Comp. Biol., 7:e1002195, 2011.

Bioconda:

**Grüning, Björn, Ryan Dale, Andreas Sjödin, Brad A. Chapman, Jillian Rowe, Christopher H. Tomkins-Tinch, Renan Valieris, the Bioconda Team, and Johannes Köster.** 2018. *Bioconda: Sustainable and Comprehensive Software Distribution for the Life Sciences.* Nature Methods, 2018 doi:10.1038/s41592-018-0046-7.

Biopython:

**Cock PA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL** (2009) *Biopython: freely available Python tools for computational molecular biology and bioinformatics.* Bioinformatics, 25, 1422-1423




