# mTAGs: accurate OTU level taxonomic profiling of metagenomes using full-length rRNA degenerate consensus references

Profiles metagenomes by finding rRNA sequences and annotating them using a degenerate consensus references for higher sensitivity.



## Overview


The tool is designed by Hans-Joachim Ruscheweyh and Guillem Salazar and distributed under the GPLv3 license. 

Questions/Comments? Write a github issue.

If you use mTAGs in a published work, please cite:

```
Guillem Salazar, Hans-Joachim Ruscheweyh, TODO, Shinichi Sunagawa. "mTAGs: accurate OTU 
level taxonomic profiling of metagenomes using full-length rRNA degenerate consensus 
references"
```



## Installation:


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



The actual can then be installed via `pip`:

```
$ git clone https://github.com/SushiLab/mTAGs.git
$ cd mTAGs
$ pip install -r requirements.txt -e . 

# Test if installation was successful

$ mtags

usage: mtags <command> [<args>]
    Command options
        download    Download the mTAGs database
        find        Searches for potential rRNA reads in metagenomic samples. (First step)
        annotate    Annotate and quantify rRNA reads. (Second Step)
        merge       Merge profiles

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

The tool is split in to three steps. The first step uses HMM models to extract potential rRNA sequences from metagenomic data. The second step aligns these sequences against a modified SILVA database and annotates sequences taxonomically. The third step is a function that merges profiles.


### 1. Find rRNA sequences in metagenomic data

This step uses precomputed HMM models to extract rRNA sequences from a metagenomic sample. 

The output are 6 fasta files, one for each molecule (ssu/lsu) and domain (archaea, bacteria, eukaryotes).

```bash

$ mtags find --help
2020-07-03 13:49:00,737 INFO: Starting mTAGs
2020-07-03 13:49:00,737 INFO: Executing command find
usage: mtags [-h] -i INPUT -o OUTPUT [-t {1,2,4,8}]

Takes a sequence file and extracts sequences that likely come from SSU/LSU
genes.

optional arguments:
  -h, --help    show this help message and exit
  -i INPUT      Input fq/fa file. Can be gzipped.
  -o OUTPUT     Output folder.
  -t {1,2,4,8}  Number of threads for hmmsearch.
```

```bash

# Input files:
# sample.1.fq.gz
# sample.2.fq.gz

# Output files:

# output_folder/sample.1.fq.gz_arc_lsu.fasta
# output_folder/sample.1.fq.gz_bac_lsu.fasta
# output_folder/sample.1.fq.gz_euk_lsu.fasta
# output_folder/sample.1.fq.gz_arc_ssu.fasta
# output_folder/sample.1.fq.gz_bac_ssu.fasta
# output_folder/sample.1.fq.gz_euk_ssu.fasta

# output_folder/sample.2.fq.gz_arc_lsu.fasta
# output_folder/sample.2.fq.gz_bac_lsu.fasta
# output_folder/sample.2.fq.gz_euk_lsu.fasta
# output_folder/sample.2.fq.gz_arc_ssu.fasta
# output_folder/sample.2.fq.gz_bac_ssu.fasta
# output_folder/sample.2.fq.gz_euk_ssu.fasta

$ mtags find -i sample.1.fq.gz -o output_folder -t 4
2020-07-03 13:52:28,100 INFO: Starting mTAGs
2020-07-03 13:52:28,100 INFO: Executing command find
2020-07-03 13:52:28,101 INFO: Extracting FastA and revcomp FastA from sample.1.fq.gz
2020-07-03 13:52:29,459 INFO: Processed reads:	100000
2020-07-03 13:52:30,908 INFO: Processed reads:	200000
2020-07-03 13:52:32,266 INFO: Processed reads:	300000
2020-07-03 13:52:33,567 INFO: Processed reads:	400000
2020-07-03 13:52:34,179 INFO: Finished extracting. Found 445320 sequences.
2020-07-03 13:52:34,179 INFO: Start detecting rRNA sequences in FastA files
2020-07-03 13:52:34,179 INFO: Start detecting rRNA sequences for molecule=ssu
2020-07-03 13:52:34,179 INFO: Executing:	hmmsearch --cpu 4 -o output_folder/sample.1.fq.gz_fw.fasta_ssu.hmmer --domtblout output_folder/sample.1.fq.gz_fw.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm output_folder/sample.1.fq.gz_fw.fasta
2020-07-03 13:52:41,437 INFO: Finished hmmsearch
2020-07-03 13:52:41,444 INFO: Executing:	hmmsearch --cpu 4 -o output_folder/sample.1.fq.gz_rev.fasta_ssu.hmmer --domtblout output_folder/sample.1.fq.gz_rev.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm output_folder/sample.1.fq.gz_rev.fasta
2020-07-03 13:52:48,655 INFO: Finished hmmsearch
2020-07-03 13:52:48,661 INFO: Finished detecting rRNA sequences for molecule=ssu
2020-07-03 13:52:48,661 INFO: Start detecting rRNA sequences for molecule=lsu
2020-07-03 13:52:48,661 INFO: Executing:	hmmsearch --cpu 4 -o output_folder/sample.1.fq.gz_fw.fasta_lsu.hmmer --domtblout output_folder/sample.1.fq.gz_fw.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm output_folder/sample.1.fq.gz_fw.fasta
2020-07-03 13:53:12,704 INFO: Finished hmmsearch
2020-07-03 13:53:12,718 INFO: Executing:	hmmsearch --cpu 4 -o output_folder/sample.1.fq.gz_rev.fasta_lsu.hmmer --domtblout output_folder/sample.1.fq.gz_rev.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm output_folder/sample.1.fq.gz_rev.fasta
2020-07-03 13:53:38,879 INFO: Finished hmmsearch
2020-07-03 13:53:38,895 INFO: Finished detecting rRNA sequences for molecule=lsu
2020-07-03 13:53:38,895 INFO: Found 4851 potential rRNA sequences.
2020-07-03 13:53:38,895 INFO: Finished detecting rRNA sequences from FastA files.
2020-07-03 13:53:38,895 INFO: Finding best molecule for each read
2020-07-03 13:53:38,904 INFO: Finished finding best molecule for each read
2020-07-03 13:53:38,904 INFO: Start extracting reads/writing output
2020-07-03 13:53:39,275 INFO: Processed reads:	100000
2020-07-03 13:53:39,644 INFO: Processed reads:	200000
2020-07-03 13:53:40,011 INFO: Processed reads:	300000
2020-07-03 13:53:40,382 INFO: Processed reads:	400000
2020-07-03 13:53:40,549 INFO: Finished extracting reads/writing output
2020-07-03 13:53:40,549 INFO: euk_lsu	2640
2020-07-03 13:53:40,549 INFO: euk_ssu	1332
2020-07-03 13:53:40,549 INFO: bac_lsu	617
2020-07-03 13:53:40,549 INFO: bac_ssu	248
2020-07-03 13:53:40,549 INFO: arc_lsu	12
2020-07-03 13:53:40,549 INFO: arc_ssu	2
2020-07-03 13:53:40,578 INFO: Finishing mTAGs

$ mtags find -i sample.2.fq.gz -o output_folder -t 4
2020-07-03 13:55:22,562 INFO: Starting mTAGs
2020-07-03 13:55:22,562 INFO: Executing command find
2020-07-03 13:55:22,562 INFO: Extracting FastA and revcomp FastA from sample.2.fq.gz
2020-07-03 13:55:23,908 INFO: Processed reads:	100000
2020-07-03 13:55:25,166 INFO: Processed reads:	200000
2020-07-03 13:55:26,431 INFO: Processed reads:	300000
2020-07-03 13:55:27,799 INFO: Processed reads:	400000
2020-07-03 13:55:28,403 INFO: Finished extracting. Found 445320 sequences.
2020-07-03 13:55:28,403 INFO: Start detecting rRNA sequences in FastA files
2020-07-03 13:55:28,403 INFO: Start detecting rRNA sequences for molecule=ssu
2020-07-03 13:55:28,403 INFO: Executing:	hmmsearch --cpu 4 -o output_folder/sample.2.fq.gz_fw.fasta_ssu.hmmer --domtblout output_folder/sample.2.fq.gz_fw.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm output_folder/sample.2.fq.gz_fw.fasta
2020-07-03 13:55:35,746 INFO: Finished hmmsearch
2020-07-03 13:55:35,753 INFO: Executing:	hmmsearch --cpu 4 -o output_folder/sample.2.fq.gz_rev.fasta_ssu.hmmer --domtblout output_folder/sample.2.fq.gz_rev.fasta_ssu.dom -E 0.01 mTAGs/data/ssu.hmm output_folder/sample.2.fq.gz_rev.fasta
2020-07-03 13:55:41,823 INFO: Finished hmmsearch
2020-07-03 13:55:41,834 INFO: Finished detecting rRNA sequences for molecule=ssu
2020-07-03 13:55:41,834 INFO: Start detecting rRNA sequences for molecule=lsu
2020-07-03 13:55:41,834 INFO: Executing:	hmmsearch --cpu 4 -o output_folder/sample.2.fq.gz_fw.fasta_lsu.hmmer --domtblout output_folder/sample.2.fq.gz_fw.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm output_folder/sample.2.fq.gz_fw.fasta
2020-07-03 13:56:03,676 INFO: Finished hmmsearch
2020-07-03 13:56:03,689 INFO: Executing:	hmmsearch --cpu 4 -o output_folder/sample.2.fq.gz_rev.fasta_lsu.hmmer --domtblout output_folder/sample.2.fq.gz_rev.fasta_lsu.dom -E 0.01 mTAGs/data/lsu.hmm output_folder/sample.2.fq.gz_rev.fasta
2020-07-03 13:56:25,024 INFO: Finished hmmsearch
2020-07-03 13:56:25,039 INFO: Finished detecting rRNA sequences for molecule=lsu
2020-07-03 13:56:25,039 INFO: Found 4739 potential rRNA sequences.
2020-07-03 13:56:25,039 INFO: Finished detecting rRNA sequences from FastA files.
2020-07-03 13:56:25,040 INFO: Finding best molecule for each read
2020-07-03 13:56:25,047 INFO: Finished finding best molecule for each read
2020-07-03 13:56:25,047 INFO: Start extracting reads/writing output
2020-07-03 13:56:25,424 INFO: Processed reads:	100000
2020-07-03 13:56:25,797 INFO: Processed reads:	200000
2020-07-03 13:56:26,176 INFO: Processed reads:	300000
2020-07-03 13:56:26,547 INFO: Processed reads:	400000
2020-07-03 13:56:26,720 INFO: Finished extracting reads/writing output
2020-07-03 13:56:26,720 INFO: euk_lsu	2567
2020-07-03 13:56:26,720 INFO: euk_ssu	1315
2020-07-03 13:56:26,720 INFO: bac_lsu	600
2020-07-03 13:56:26,720 INFO: bac_ssu	238
2020-07-03 13:56:26,720 INFO: arc_lsu	17
2020-07-03 13:56:26,720 INFO: arc_ssu	2
2020-07-03 13:56:26,741 INFO: Finishing mTAGs

```



### 2. Annotate SSU sequences using the Silva database


In this step you use the previously generated sequence files and align them against a modified Silva database and extract the most meaningful alignments using a percent identy/coverage cutoff and applying an LCA.


This examples shows usage for paired-end data. The program also works for single-end or merged inserts.


```bash
$ mtags annotate --help
2020-07-03 14:10:16,022 INFO: Starting mTAGs
2020-07-03 14:10:16,022 INFO: Executing command annotate
usage: mtags [-h] [-i1 I1 [I1 ...]] [-i2 I2 [I2 ...]] [-is IM [IM ...]] -o O
             [-t T] [-ma MA] [-mr MR] [-b]

Assigns LCA based taxonomy to rRNA ssu reads.

optional arguments:
  -h, --help       show this help message and exit
  -i1 I1 [I1 ...]  Input R1 fasta file(s). Sequences in paired files will be
                   matched by name. Unmatched sequences will be counted as
                   singletons.
  -i2 I2 [I2 ...]  Input R2 fasta file(s). Sequences in paired files will be
                   matched by name. Unmatched sequences will be counted as
                   singletons.
  -is IM [IM ...]  Input singletons/merged fasta file(s).
  -o O             Output pattern. This pattern will be prepended to all
                   output files
  -t T             Number of threads for vsearch alignment
  -ma MA           Maxaccepts, vsearch parameter. Larger numbers increase
                   sensitivity and runtime. default=1000
  -mr MR           Maxrejects, vsearch parameter. Larger numbers increase
                   sensitivity and runtime. default=100
  -b               Enable binning - write taxonomic assignment for each
                   read/insert
```

```
# Input files

# output_folder/sample.1.fq.gz_arc_ssu.fasta
# output_folder/sample.1.fq.gz_bac_ssu.fasta
# output_folder/sample.1.fq.gz_euk_ssu.fasta

# output_folder/sample.2.fq.gz_arc_ssu.fasta
# output_folder/sample.2.fq.gz_bac_ssu.fasta
# output_folder/sample.2.fq.gz_euk_ssu.fasta

$ mtags annotate -i1 output_folder/sample.1.fq.gz_*ssu.fasta -i2 output_folder/sample.2.fq.gz_*ssu.fasta -o sample -ma 100 -mr 100 -t 4 -b
2020-07-03 14:11:01,288 INFO: Starting mTAGs
2020-07-03 14:11:01,288 INFO: Executing command annotate
2020-07-03 14:11:01,289 INFO: Reading and testing FastA files
2020-07-03 14:11:01,299 INFO: Paired files stats:
2020-07-03 14:11:01,299 INFO: 	Found 1859 inserts in paired files of which 1278 are actually paired and 581 are singletons
2020-07-03 14:11:01,299 INFO: Singleton/Merged files stats:
2020-07-03 14:11:01,299 INFO: 	Found 0 inserts in singleton files.
2020-07-03 14:11:01,300 INFO: Total stats:
2020-07-03 14:11:01,300 INFO: 	Paired inserts:	1278
2020-07-03 14:11:01,300 INFO: 	Singleton inserts:	581
2020-07-03 14:11:01,300 INFO: Finished reading and testing FastA files
2020-07-03 14:11:01,300 INFO: Writing paired sequences to FastA:	sample.pe.fasta
2020-07-03 14:11:01,302 INFO: Writing finished
2020-07-03 14:11:01,303 INFO: Writing singleton sequences to FastA:	sample.se.fasta
2020-07-03 14:11:01,303 INFO: Writing finished
2020-07-03 14:11:01,304 INFO: Aligning paired sequence file:	vsearch --usearch_global sample.pe.fasta --db mTAGs/db/SILVA-138_NR-97_complink_cons.vsearch.udb --id 0.97 --maxaccepts 100 --maxrejects 100 --strand both --userout sample.pe.m8 --userfields query+target+qcov+id+qilo+qihi+tilo+tihi+alnlen+ids+mism+gaps+qstrand --query_cov 0.5 --threads 4 --output_no_hits
vsearch v2.15.0_macos_x86_64, 16.0GB RAM, 4 cores
https://github.com/torognes/vsearch

Reading UDB file mTAGs/db/SILVA-138_NR-97_complink_cons.vsearch.udb 100%
Reorganizing data in memory 100%
Creating bitmaps 100%
Parsing abundances 100%
341380158 nt in 233519 seqs, min 900, max 3718, avg 1462
Searching 100%
Matching unique query sequences: 1994 of 2556 (78.01%)
2020-07-03 14:11:36,989 INFO: Finished alignment
2020-07-03 14:11:36,990 INFO: Aligning singleton sequence file:	vsearch --usearch_global sample.se.fasta --db mTAGs/db/SILVA-138_NR-97_complink_cons.vsearch.udb --id 0.97 --maxaccepts 100 --maxrejects 100 --strand both --userout sample.se.m8 --userfields query+target+qcov+id+qilo+qihi+tilo+tihi+alnlen+ids+mism+gaps+qstrand --query_cov 0.5 --threads 4 --output_no_hits
vsearch v2.15.0_macos_x86_64, 16.0GB RAM, 4 cores
https://github.com/torognes/vsearch

Reading UDB file mTAGs/db/SILVA-138_NR-97_complink_cons.vsearch.udb 100%
Reorganizing data in memory 100%
Creating bitmaps 100%
Parsing abundances 100%
341380158 nt in 233519 seqs, min 900, max 3718, avg 1462
Searching 100%
Matching unique query sequences: 349 of 581 (60.07%)
2020-07-03 14:11:46,846 INFO: Finished alignment
2020-07-03 14:11:46,846 INFO: Parsing alignment files
2020-07-03 14:11:48,206 INFO: Finished parsing alignment files
2020-07-03 14:11:48,206 INFO: Filtering alignments by score and length
2020-07-03 14:11:48,375 INFO: Finished filtering
2020-07-03 14:11:48,375 INFO: Performing LCA calculations
2020-07-03 14:11:50,629 INFO: Finished LCA calculations
2020-07-03 14:11:50,630 INFO: Start writing output files
2020-07-03 14:11:50,630 INFO: Start writing:	sample.bins
2020-07-03 14:11:50,986 INFO: Start writing:	sample.otu.tsv
2020-07-03 14:11:51,087 INFO: Start writing:	sample.genus.tsv
2020-07-03 14:11:51,091 INFO: Start writing:	sample.family.tsv
2020-07-03 14:11:51,092 INFO: Start writing:	sample.order.tsv
2020-07-03 14:11:51,093 INFO: Start writing:	sample.class.tsv
2020-07-03 14:11:51,093 INFO: Start writing:	sample.phylum.tsv
2020-07-03 14:11:51,094 INFO: Start writing:	sample.domain.tsv
2020-07-03 14:11:51,094 INFO: Start writing:	sample.root.tsv
2020-07-03 14:11:51,094 INFO: Finished writing output files
2020-07-03 14:11:51,178 INFO: Finished annotate command
2020-07-03 14:11:51,178 INFO: Finishing mTAGs
```

### 2. Merge profiles

The previous step produces profiles for different samples. For comparision it is useful to have multiple profiles in one file. For that we provide the merge function:

```bash
$  mtags merge --help
2020-07-07 16:34:59,447 INFO: Starting mTAGs
2020-07-07 16:34:59,447 INFO: Executing command merge
usage: mtags [-h] -i BINS [BINS ...] -o O

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
2020-07-07 16:36:13,705 INFO: Finishing mTAGs with status:	0

```

## The Output

The tool produces two types of output, the taxonomic profiles and the annotation of every insert.

### The profile

The tool produces by default profiles for 8 taxonomic ranks:

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

On the otu level there is the complete taxonomic path plus the name of the otu (e.g. `otu__silva_138_complink_cons_otu_215825`) and the number of inserts that map to this otu. The two last rows represent all sequences that could not be aligned or that could be aligned but could not be assigned at this taxonomic level. 


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




