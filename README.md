# RNAprep
Process for taking a single sample of Illumina RNA expression data, in fastq format, and processing it into FPKM values that are comparable to publicly available datasets.

## Reference Files
The **hg19** (GRCh37) build of the USCS human reference genome can be found [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/), which should be downloaded in 2bit format. The utility program for Linux distributions, twoBitToFa, can be used to extract .fa file(s) from this 2bit binary. A precompiled version of this tool can be found [here](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/).

The **comprehensive gene annotation list** can be found [here](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19). These are GENCODE annotations specific to the hg19 build of the human genome reference chromosomes and were downloaded in gtf format. 

## Tools
..* [STAR](https://github.com/alexdobin/STAR/releases/tag/STAR_2.4.2a) aligner v2.4.2
..* RSEM was used to quantify expression and can be cloned from this [git repo](https://github.com/deweylab/RSEM)

NOTE: If you plan to append the sample of RNAseq data to an existing matrix of samples, [GEMprep](https://github.com/SystemsGenetics/GEMprep) offers an easy method of performing quantile normalization.
