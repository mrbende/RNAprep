# RNAprep
Process for taking a single sample of Illumina RNA expression data, in fastq format, and processing it into FPKM values that are comparable to publicly available datgasets.

# Reference Files
The **hg19** (GRCh37) build of the USCS human reference genome can be found [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/), which should be downloaded in 2bit format. The utility program for Linux distributions, twoBitToFa, can be used to extract .fa file(s) from this 2bit binary. A precompiled version of this tool can be found [here](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/).

The **comprehensive gene annotation list** can be downloaded via a wget to this download [link](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19). These are GENCODE annotations specific to the hg19 build of the human genome reference chromosomes and were downloaded in gtf format. 
