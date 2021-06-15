## mTAGs: taxonomic profiling using degenerate <br> consensus reference sequences of ribosomal RNA genes

<p align="center">
<img src="img/mtags_logo.png" width="500"  />
</p>

Profile metagenomes by finding rRNA sequences and annotating them using a degenerate consensus references for higher sensitivity.


- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
  * [Profile](#profile)
  * [Merge profiles](#merge-profiles)
- [Output](#output)
  * [Profile file](#profile-file)
  * [Single insert annotation](#single-insert-annotation)
- [Database](#database)
- [References](#references)




## Overview


The tool is designed by Hans-Joachim Ruscheweyh and Guillem Salazar and distributed under the GPLv3 license. 

Questions/Comments? We're happy to help you via GitHub issues.

If you use mTAGs in a published work, please cite:

```
Guillem Salazar, Hans-Joachim Ruscheweyh, TODO, Shinichi Sunagawa. "mTAGs: accurate OTU 
level taxonomic profiling of metagenomes using full-length rRNA degenerate consensus 
references"
```



## Installation


mTAGs is written in Python3 and can be installed via `pip`. The tool has 3 dependencies that need to be installed in advance:




1. `vsearch>=2.14.0`
2. `Python>=3.7`
3. `hmmer>=3.3`

This can be done manually or using conda:

```
$ conda create -n mtags python=3.7 hmmer vsearch
$ source activate mtags
# or
$ conda activate mtags
```



mTAGs can then be installed via `pip`:

```
$ git clone https://github.com/SushiLab/mTAGs.git
$ cd mTAGs
$ pip install -r requirements.txt -e . 

# Test if installation was successful

$ mtags

usage: mtags <command> [<args>]

    [General]

        profile     Extract and taxonomically annotate rRNA reads in metagenomic samples
        merge       Merge profiles

    [Expert]

        extract     Extract rRNA reads in metagenomic samples
        annotate    Annotate and quantify rRNA reads

    [Installation]

        download    Download the mTAGs database - Once after download of the tool

```

The database needs to be downloaded in the last step of the installation. This needs to be done once and before the first metagenomic samples can be processed:

```
$ mtags download

2020-07-03 12:02:40,294 INFO: Starting mTAGs
2020-07-03 12:02:40,294 INFO: Start downloading the mTAGs database. ~600MB
2020-07-03 12:05:17,883 INFO: Finished downloading the mTAGs database.
2020-07-03 12:05:17,883 INFO: Finishing mTAGs
```



## Usage

The tool is split in to two steps. The first step uses HMM models to **extract** potential rRNA sequences from metagenomic data and aligns these sequences against a modified SILVA database and **annotate**s sequences taxonomically. The second step is a function that **merge**s profiles. (The steps **extract** and **annotate** are grouped into a single command **profile**)



### Profile

This step uses precomputed HMM models to extract rRNA sequences from a metagenomic sample. The rRNA sequences are then aligned against a clustered Silva database to annotate sequences and profile samples.



```bash

$ mtags profile
2020-08-31 17:58:18,259 INFO: Starting mTAGs
2020-08-31 17:58:18,260 INFO: Executing command profile
usage: mtags [-h] [-i1 I1 [I1 ...]] [-i2 I2 [I2 ...]] [-is IM [IM ...]] -o
             OUTPUT -s SAMPLE [-t THREADS] [-ma MA] [-mr MR]

Extract and taxonomically annotate rRNA reads in metagenomic samples

optional arguments:
  -h, --help       show this help message and exit
  -i1 I1 [I1 ...]  Input forward reads files. Can be fasta/fastq and gzipped.
                   The order of the files has to match the order in -i2. Read
                   pairs in files are identified by the identical name in -i2.
  -i2 I2 [I2 ...]  Input reverse reads files. Can be fasta/fastq and gzipped.
                   The order of the files has to match the order in -i1. Read
                   pairs in files are identified by the identical name in -i1.
  -is IM [IM ...]  Input merged/single-end reads files. Can be fasta/fastq and
                   gzipped.
  -o OUTPUT        Output folder.
  -s SAMPLE        Samplename
  -t THREADS       Number of threads.
  -ma MA           Maxaccepts, vsearch parameter. Larger numbers increase
                   sensitivity and runtime. default=1000
  -mr MR           Maxrejects, vsearch parameter. Larger numbers increase
                   sensitivity and runtime. default=100
2020-08-31 17:58:18,263 INFO: Finishing mTAGs with status:


$ mtags profile -i1 sample.1.fq.gz -i2 sample.2.fq.gz -is sample.s.fq.gz sample.m.fq.gz -o output -t 4 -s sample -ma 100 -mr 100
2020-08-31 18:00:09,940 INFO: Starting mTAGs
2020-08-31 18:00:09,940 INFO: Executing command profile
2020-08-31 18:00:09,942 INFO: Extracting FastA and revcomp FastA from sample.1.fq.gz
2020-08-31 18:00:20,207 INFO: Processed reads:	1000000
2020-08-31 18:00:30,471 INFO: Processed reads:	2000000
2020-08-31 18:00:40,969 INFO: Processed reads:	3000000
2020-08-31 18:00:51,208 INFO: Processed reads:	4000000
2020-08-31 18:00:53,419 INFO: Processed reads:	4223475
2020-08-31 18:00:53,419 INFO: Finished extracting. Found 4223475 sequences.
2020-08-31 18:00:53,419 INFO: Start detecting rRNA sequences in FastA files
2020-08-31 18:00:53,419 INFO: Start detecting rRNA sequences for molecule=ssu
2020-08-31 18:00:53,419 INFO: Executing:	hmmsearch --cpu 4 -o output/sample.1.fq.gz_fw.fasta_ssu.hmmer --domtblout output/sample.1.fq.gz_fw.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm output/sample.1.fq.gz_fw.fasta
2020-08-31 18:01:24,679 INFO: Finished hmmsearch
2020-08-31 18:01:24,775 INFO: Executing:	hmmsearch --cpu 4 -o output/sample.1.fq.gz_rev.fasta_ssu.hmmer --domtblout output/sample.1.fq.gz_rev.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm output/sample.1.fq.gz_rev.fasta
2020-08-31 18:01:57,804 INFO: Finished hmmsearch
2020-08-31 18:01:57,898 INFO: Finished detecting rRNA sequences for molecule=ssu
2020-08-31 18:01:57,898 INFO: Start detecting rRNA sequences for molecule=lsu
2020-08-31 18:01:57,898 INFO: Executing:	hmmsearch --cpu 4 -o output/sample.1.fq.gz_fw.fasta_lsu.hmmer --domtblout output/sample.1.fq.gz_fw.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm output/sample.1.fq.gz_fw.fasta
2020-08-31 18:03:48,873 INFO: Finished hmmsearch
2020-08-31 18:03:48,960 INFO: Executing:	hmmsearch --cpu 4 -o output/sample.1.fq.gz_rev.fasta_lsu.hmmer --domtblout output/sample.1.fq.gz_rev.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm output/sample.1.fq.gz_rev.fasta
2020-08-31 18:05:30,688 INFO: Finished hmmsearch
2020-08-31 18:05:30,808 INFO: Finished detecting rRNA sequences for molecule=lsu
2020-08-31 18:05:30,808 INFO: Found 56171 potential rRNA sequences.
2020-08-31 18:05:30,808 INFO: Finished detecting rRNA sequences from FastA files.
2020-08-31 18:05:30,808 INFO: Finding best molecule for each read
2020-08-31 18:05:30,888 INFO: Finished finding best molecule for each read
2020-08-31 18:05:30,888 INFO: Start extracting reads/writing output
2020-08-31 18:05:34,105 INFO: Processed reads:	1000000
2020-08-31 18:05:37,198 INFO: Processed reads:	2000000
2020-08-31 18:05:40,274 INFO: Processed reads:	3000000
2020-08-31 18:05:43,292 INFO: Processed reads:	4000000
2020-08-31 18:05:43,981 INFO: Processed reads:	4223475
2020-08-31 18:05:43,982 INFO: Finished extracting reads/writing output
2020-08-31 18:05:43,983 INFO: bac_lsu	37112
2020-08-31 18:05:43,983 INFO: bac_ssu	18877
2020-08-31 18:05:43,983 INFO: arc_ssu	97
2020-08-31 18:05:43,983 INFO: arc_lsu	83
2020-08-31 18:05:43,983 INFO: euk_lsu	2
2020-08-31 18:05:44,066 INFO: Extracting FastA and revcomp FastA from sample.2.fq.gz
2020-08-31 18:05:54,509 INFO: Processed reads:	1000000
2020-08-31 18:06:04,598 INFO: Processed reads:	2000000
2020-08-31 18:06:14,645 INFO: Processed reads:	3000000
2020-08-31 18:06:25,100 INFO: Processed reads:	4000000
2020-08-31 18:06:27,371 INFO: Processed reads:	4223475
2020-08-31 18:06:27,371 INFO: Finished extracting. Found 4223475 sequences.
2020-08-31 18:06:27,371 INFO: Start detecting rRNA sequences in FastA files
2020-08-31 18:06:27,371 INFO: Start detecting rRNA sequences for molecule=ssu
2020-08-31 18:06:27,371 INFO: Executing:	hmmsearch --cpu 4 -o output/sample.2.fq.gz_fw.fasta_ssu.hmmer --domtblout output/sample.2.fq.gz_fw.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm output/sample.2.fq.gz_fw.fasta
2020-08-31 18:06:57,558 INFO: Finished hmmsearch
2020-08-31 18:06:57,687 INFO: Executing:	hmmsearch --cpu 4 -o output/sample.2.fq.gz_rev.fasta_ssu.hmmer --domtblout output/sample.2.fq.gz_rev.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm output/sample.2.fq.gz_rev.fasta
2020-08-31 18:07:28,056 INFO: Finished hmmsearch
2020-08-31 18:07:28,153 INFO: Finished detecting rRNA sequences for molecule=ssu
2020-08-31 18:07:28,153 INFO: Start detecting rRNA sequences for molecule=lsu
2020-08-31 18:07:28,153 INFO: Executing:	hmmsearch --cpu 4 -o output/sample.2.fq.gz_fw.fasta_lsu.hmmer --domtblout output/sample.2.fq.gz_fw.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm output/sample.2.fq.gz_fw.fasta
2020-08-31 18:09:08,830 INFO: Finished hmmsearch
2020-08-31 18:09:08,938 INFO: Executing:	hmmsearch --cpu 4 -o output/sample.2.fq.gz_rev.fasta_lsu.hmmer --domtblout output/sample.2.fq.gz_rev.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm output/sample.2.fq.gz_rev.fasta
2020-08-31 18:10:52,109 INFO: Finished hmmsearch
2020-08-31 18:10:52,198 INFO: Finished detecting rRNA sequences for molecule=lsu
2020-08-31 18:10:52,199 INFO: Found 55897 potential rRNA sequences.
2020-08-31 18:10:52,199 INFO: Finished detecting rRNA sequences from FastA files.
2020-08-31 18:10:52,199 INFO: Finding best molecule for each read
2020-08-31 18:10:52,277 INFO: Finished finding best molecule for each read
2020-08-31 18:10:52,278 INFO: Start extracting reads/writing output
2020-08-31 18:10:55,212 INFO: Processed reads:	1000000
2020-08-31 18:10:58,132 INFO: Processed reads:	2000000
2020-08-31 18:11:01,049 INFO: Processed reads:	3000000
2020-08-31 18:11:03,975 INFO: Processed reads:	4000000
2020-08-31 18:11:04,629 INFO: Processed reads:	4223475
2020-08-31 18:11:04,636 INFO: Finished extracting reads/writing output
2020-08-31 18:11:04,636 INFO: bac_lsu	36956
2020-08-31 18:11:04,636 INFO: bac_ssu	18786
2020-08-31 18:11:04,636 INFO: arc_ssu	83
2020-08-31 18:11:04,636 INFO: arc_lsu	67
2020-08-31 18:11:04,636 INFO: euk_lsu	5
2020-08-31 18:11:04,752 INFO: Extracting FastA and revcomp FastA from sample.s.fq.gz
2020-08-31 18:11:15,075 INFO: Processed reads:	1000000
2020-08-31 18:11:15,153 INFO: Processed reads:	1007313
2020-08-31 18:11:15,153 INFO: Finished extracting. Found 1007313 sequences.
2020-08-31 18:11:15,153 INFO: Start detecting rRNA sequences in FastA files
2020-08-31 18:11:15,153 INFO: Start detecting rRNA sequences for molecule=ssu
2020-08-31 18:11:15,154 INFO: Executing:	hmmsearch --cpu 4 -o output/sample.s.fq.gz_fw.fasta_ssu.hmmer --domtblout output/sample.s.fq.gz_fw.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm output/sample.s.fq.gz_fw.fasta
2020-08-31 18:11:21,030 INFO: Finished hmmsearch
2020-08-31 18:11:21,049 INFO: Executing:	hmmsearch --cpu 4 -o output/sample.s.fq.gz_rev.fasta_ssu.hmmer --domtblout output/sample.s.fq.gz_rev.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm output/sample.s.fq.gz_rev.fasta
2020-08-31 18:11:27,699 INFO: Finished hmmsearch
2020-08-31 18:11:27,716 INFO: Finished detecting rRNA sequences for molecule=ssu
2020-08-31 18:11:27,716 INFO: Start detecting rRNA sequences for molecule=lsu
2020-08-31 18:11:27,716 INFO: Executing:	hmmsearch --cpu 4 -o output/sample.s.fq.gz_fw.fasta_lsu.hmmer --domtblout output/sample.s.fq.gz_fw.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm output/sample.s.fq.gz_fw.fasta
2020-08-31 18:11:45,314 INFO: Finished hmmsearch
2020-08-31 18:11:45,336 INFO: Executing:	hmmsearch --cpu 4 -o output/sample.s.fq.gz_rev.fasta_lsu.hmmer --domtblout output/sample.s.fq.gz_rev.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm output/sample.s.fq.gz_rev.fasta
2020-08-31 18:12:08,698 INFO: Finished hmmsearch
2020-08-31 18:12:08,720 INFO: Finished detecting rRNA sequences for molecule=lsu
2020-08-31 18:12:08,720 INFO: Found 12344 potential rRNA sequences.
2020-08-31 18:12:08,720 INFO: Finished detecting rRNA sequences from FastA files.
2020-08-31 18:12:08,720 INFO: Finding best molecule for each read
2020-08-31 18:12:08,738 INFO: Finished finding best molecule for each read
2020-08-31 18:12:08,738 INFO: Start extracting reads/writing output
2020-08-31 18:12:11,616 INFO: Processed reads:	1000000
2020-08-31 18:12:11,638 INFO: Processed reads:	1007313
2020-08-31 18:12:11,639 INFO: Finished extracting reads/writing output
2020-08-31 18:12:11,639 INFO: bac_lsu	8127
2020-08-31 18:12:11,639 INFO: bac_ssu	4153
2020-08-31 18:12:11,639 INFO: arc_ssu	33
2020-08-31 18:12:11,639 INFO: arc_lsu	29
2020-08-31 18:12:11,639 INFO: euk_lsu	2
2020-08-31 18:12:11,659 INFO: Extracting FastA and revcomp FastA from sample.m.fq.gz
2020-08-31 18:12:22,009 INFO: Processed reads:	1000000
2020-08-31 18:12:32,164 INFO: Processed reads:	2000000
2020-08-31 18:12:42,332 INFO: Processed reads:	3000000
2020-08-31 18:12:52,700 INFO: Processed reads:	4000000
2020-08-31 18:13:02,863 INFO: Processed reads:	5000000
2020-08-31 18:13:11,104 INFO: Processed reads:	5802126
2020-08-31 18:13:11,104 INFO: Finished extracting. Found 5802126 sequences.
2020-08-31 18:13:11,104 INFO: Start detecting rRNA sequences in FastA files
2020-08-31 18:13:11,104 INFO: Start detecting rRNA sequences for molecule=ssu
2020-08-31 18:13:11,104 INFO: Executing:	hmmsearch --cpu 4 -o output/sample.m.fq.gz_fw.fasta_ssu.hmmer --domtblout output/sample.m.fq.gz_fw.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm output/sample.m.fq.gz_fw.fasta
2020-08-31 18:13:58,458 INFO: Finished hmmsearch
2020-08-31 18:13:58,522 INFO: Executing:	hmmsearch --cpu 4 -o output/sample.m.fq.gz_rev.fasta_ssu.hmmer --domtblout output/sample.m.fq.gz_rev.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm output/sample.m.fq.gz_rev.fasta
2020-08-31 18:14:44,553 INFO: Finished hmmsearch
2020-08-31 18:14:44,679 INFO: Finished detecting rRNA sequences for molecule=ssu
2020-08-31 18:14:44,679 INFO: Start detecting rRNA sequences for molecule=lsu
2020-08-31 18:14:44,680 INFO: Executing:	hmmsearch --cpu 4 -o output/sample.m.fq.gz_fw.fasta_lsu.hmmer --domtblout output/sample.m.fq.gz_fw.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm output/sample.m.fq.gz_fw.fasta
2020-08-31 18:17:19,336 INFO: Finished hmmsearch
2020-08-31 18:17:19,456 INFO: Executing:	hmmsearch --cpu 4 -o output/sample.m.fq.gz_rev.fasta_lsu.hmmer --domtblout output/sample.m.fq.gz_rev.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm output/sample.m.fq.gz_rev.fasta
2020-08-31 18:19:50,845 INFO: Finished hmmsearch
2020-08-31 18:19:50,989 INFO: Finished detecting rRNA sequences for molecule=lsu
2020-08-31 18:19:50,989 INFO: Found 74109 potential rRNA sequences.
2020-08-31 18:19:50,989 INFO: Finished detecting rRNA sequences from FastA files.
2020-08-31 18:19:50,989 INFO: Finding best molecule for each read
2020-08-31 18:19:51,094 INFO: Finished finding best molecule for each read
2020-08-31 18:19:51,094 INFO: Start extracting reads/writing output
2020-08-31 18:19:54,149 INFO: Processed reads:	1000000
2020-08-31 18:19:57,103 INFO: Processed reads:	2000000
2020-08-31 18:20:00,127 INFO: Processed reads:	3000000
2020-08-31 18:20:03,163 INFO: Processed reads:	4000000
2020-08-31 18:20:06,130 INFO: Processed reads:	5000000
2020-08-31 18:20:08,495 INFO: Processed reads:	5802126
2020-08-31 18:20:08,497 INFO: Finished extracting reads/writing output
2020-08-31 18:20:08,497 INFO: bac_lsu	48331
2020-08-31 18:20:08,497 INFO: bac_ssu	25619
2020-08-31 18:20:08,497 INFO: arc_ssu	80
2020-08-31 18:20:08,497 INFO: arc_lsu	64
2020-08-31 18:20:08,498 INFO: euk_lsu	15
2020-08-31 18:20:08,666 INFO: Input files:
2020-08-31 18:20:08,666 INFO: 	Forward readfiles:
2020-08-31 18:20:08,666 INFO: 		output/sample.1.fq.gz_bac_ssu.fasta
2020-08-31 18:20:08,666 INFO: 		output/sample.1.fq.gz_arc_ssu.fasta
2020-08-31 18:20:08,666 INFO: 	Reverse readfiles:
2020-08-31 18:20:08,666 INFO: 		output/sample.2.fq.gz_bac_ssu.fasta
2020-08-31 18:20:08,666 INFO: 		output/sample.2.fq.gz_arc_ssu.fasta
2020-08-31 18:20:08,666 INFO: 	Single/Merged readfiles:
2020-08-31 18:20:08,666 INFO: 		output/sample.s.fq.gz_bac_ssu.fasta
2020-08-31 18:20:08,666 INFO: 		output/sample.s.fq.gz_arc_ssu.fasta
2020-08-31 18:20:08,666 INFO: 		output/sample.m.fq.gz_bac_ssu.fasta
2020-08-31 18:20:08,666 INFO: 		output/sample.m.fq.gz_arc_ssu.fasta
2020-08-31 18:20:08,666 INFO: Output pattern:
2020-08-31 18:20:08,667 INFO: 	output/sample
2020-08-31 18:20:08,667 INFO: Reading and testing FastA files
2020-08-31 18:20:08,818 INFO: Paired files stats:
2020-08-31 18:20:08,818 INFO: 	Found 20873 inserts in paired files of which 16970 are actually paired and 3903 are singletons
2020-08-31 18:20:08,818 INFO: Singleton/Merged files stats:
2020-08-31 18:20:08,818 INFO: 	Found 29885 inserts in singleton files.
2020-08-31 18:20:08,823 INFO: Total stats:
2020-08-31 18:20:08,823 INFO: 	Paired inserts:	16970
2020-08-31 18:20:08,823 INFO: 	Singleton inserts:	33788
2020-08-31 18:20:08,823 INFO: Finished reading and testing FastA files
2020-08-31 18:20:08,824 INFO: Writing paired sequences to FastA:	output/sample.pe.fasta
2020-08-31 18:20:08,842 INFO: Writing finished
2020-08-31 18:20:08,842 INFO: Writing singleton sequences to FastA:	output/sample.se.fasta
2020-08-31 18:20:08,862 INFO: Writing finished
2020-08-31 18:20:08,865 INFO: Aligning paired sequence file:	vsearch --usearch_global output/sample.pe.fasta --db mTAGs/db/SILVA-138_NR-97_complink_cons.vsearch.udb --id 0.97 --maxaccepts 100 --maxrejects 100 --strand both --userout output/sample.pe.m8 --userfields query+target+qcov+id+qilo+qihi+tilo+tihi+alnlen+ids+mism+gaps+qstrand --query_cov 0.5 --threads 4 --output_no_hits
vsearch v2.15.0_macos_x86_64, 16.0GB RAM, 4 cores
https://github.com/torognes/vsearch

Reading UDB file mTAGs/db/SILVA-138_NR-97_complink_cons.vsearch.udb 100%
Reorganizing data in memory 100%
Creating bitmaps 100%
Parsing abundances 100%
341380158 nt in 233519 seqs, min 900, max 3718, avg 1462
Searching 100%
Matching unique query sequences: 33814 of 33940 (99.63%)
2020-08-31 18:21:58,616 INFO: Finished alignment
2020-08-31 18:21:58,616 INFO: Aligning singleton sequence file:	vsearch --usearch_global output/sample.se.fasta --db mTAGs/db/SILVA-138_NR-97_complink_cons.vsearch.udb --id 0.97 --maxaccepts 100 --maxrejects 100 --strand both --userout output/sample.se.m8 --userfields query+target+qcov+id+qilo+qihi+tilo+tihi+alnlen+ids+mism+gaps+qstrand --query_cov 0.5 --threads 4 --output_no_hits
vsearch v2.15.0_macos_x86_64, 16.0GB RAM, 4 cores
https://github.com/torognes/vsearch

Reading UDB file mTAGs/db/SILVA-138_NR-97_complink_cons.vsearch.udb 100%
Reorganizing data in memory 100%
Creating bitmaps 100%
Parsing abundances 100%
341380158 nt in 233519 seqs, min 900, max 3718, avg 1462
Searching 100%
Matching unique query sequences: 33561 of 33788 (99.33%)
2020-08-31 18:24:12,407 INFO: Finished alignment
2020-08-31 18:24:12,407 INFO: Parsing alignment files
2020-08-31 18:24:45,860 INFO: Finished parsing alignment files
2020-08-31 18:24:45,860 INFO: Filtering alignments by score and length
2020-08-31 18:24:50,990 INFO: Finished filtering
2020-08-31 18:24:50,990 INFO: Performing LCA calculations
2020-08-31 18:24:55,099 INFO: Finished LCA calculations
2020-08-31 18:24:55,099 INFO: Start writing output files
2020-08-31 18:24:55,099 INFO: Start writing:	output/sample.bins
2020-08-31 18:24:55,726 INFO: Start writing:	output/sample.otu.tsv
2020-08-31 18:24:55,822 INFO: Start writing:	output/sample.genus.tsv
2020-08-31 18:24:55,826 INFO: Start writing:	output/sample.family.tsv
2020-08-31 18:24:55,827 INFO: Start writing:	output/sample.order.tsv
2020-08-31 18:24:55,828 INFO: Start writing:	output/sample.class.tsv
2020-08-31 18:24:55,829 INFO: Start writing:	output/sample.phylum.tsv
2020-08-31 18:24:55,829 INFO: Start writing:	output/sample.domain.tsv
2020-08-31 18:24:55,829 INFO: Start writing:	output/sample.root.tsv
2020-08-31 18:24:55,830 INFO: Finished writing output files
2020-08-31 18:24:57,815 INFO: Finished annotate command
2020-08-31 18:24:57,816 INFO: Finishing mTAGs with status:	0

```

### Merge profiles

The previous step produces profiles for different samples. It is useful for analysis purposes to have smultiple profiles in one file. The merge function combines multiple profiles into a tab-separated file with each column representating a column:

```bash
$  mtags merge --help
2020-07-07 16:34:59,447 INFO: Starting mTAGs
2020-07-07 16:34:59,447 INFO: Executing command merge


Merge rRNA profiles generated by the mTAGs tool

optional arguments:
  -h, --help          show this help message and exit
  -i BINS [BINS ...]  The *bins files that should be merged into one table
  -o O                Output pattern. This pattern will be prepended to all
                      output files
                      
$ mtags merge -i *bins -o merged_profile

 mtags merge -i mtags2/*bins -o out
2020-07-07 16:36:06,642 INFO: Starting mTAGs
2020-07-07 16:36:06,642 INFO: Executing command merge
2020-07-07 16:36:06,676 INFO: Finished reading sample1.mTAG.bins. Found 35405 inserts.
2020-07-07 16:36:06,691 INFO: Finished reading sample2.mTAG.bins. Found 16265 inserts.
2020-07-07 16:36:06,718 INFO: Finished reading sample3.mTAG.bins. Found 25680 inserts.
2020-07-07 16:36:06,749 INFO: Finished reading sample4.mTAG.bins. Found 34954 inserts.

```



## Output

The tool produces two types of output, the taxonomic profiles and the annotation of every insert.

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

The output of a file looks the following:

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

The otu level is comprised of a complete taxonomic path plus the name of the otu (e.g. `otu__silva_138_complink_cons_otu_215825`) and the number of inserts that map to this otu. The two last rows represent all sequences that could not be aligned or that could be aligned but could not be assigned at this taxonomic level. 


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

The `.bins` file with the taxonomic annotation of each insert looks the following:

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

The default database of mTAGs is a modified version of the Silva release 138. We cluster the database at genus level using a complete linkage algorithm on pairwise distances and then build one consensus sequence per cluster. We also generated a database from Silva release 128 using the same strategy for reproducibility reasons. However, we encourage users to perform their analysis with most recent database version.

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




