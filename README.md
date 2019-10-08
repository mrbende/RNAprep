# RNAprep
Process for taking a single sample of Illumina RNA expression data, in fastq format, and processing it into FPKM values that are comparable to publicly available datasets.

## Reference Files
The **hg19** (GRCh37) build of the USCS human reference genome can be found [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/), which should be downloaded in 2bit format. The utility program for Linux distributions, twoBitToFa, can be used to extract .fa file(s) from this 2bit binary. A precompiled version of this tool can be found [here](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/).

The **comprehensive gene annotation list** can be found [here](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19). These are GENCODE annotations specific to the hg19 build of the human genome reference chromosomes and were downloaded in gtf format. 

## Tools
* [STAR](https://github.com/alexdobin/STAR/releases/tag/STAR_2.4.2a) aligner v2.4.2
* [RSEM](https://github.com/deweylab/RSEM) quantification
* [GEMprep](https://github.com/SystemsGenetics/GEMprep) repository tools

NOTE: If you plan to append the sample of RNAseq data to an existing matrix of samples, [GEMprep](https://github.com/SystemsGenetics/GEMprep) offers an easy method of performing quantile normalization.

## Process
This process was developed using Clemson University's Palmetto Cluster, which utilizes the Portable Batch Scheduling system (PBS) to manage job submission. Most commands were wrapped into independent scripts that specified resource allocation and the exact command line parameters. The code is copied below, but considering the large resource requirements of these processes it is recommended that they be submitted as batch jobs if possible. It is also important that the STAR executable as well as the RSEM commands are recognizable within the system path, or otherwise are explicitly directed.

**1. Generate the Genome Index**
```
STAR --runThreadN 24 --runMode genomeGenerate \
--genomeDir /path/to/desired/output/directory \
--genomeFastaFiles /path/to/genome/fasta/hg19.fa \
--sjdbGTFfile /path/to/annotations/gencode.v19.annotation.gtf \
--sjdbOverhang 99
```
The `--runThreadN` flag should be set to the number of avaiable cores on the node. The genome index that this process creates will be stored in a new directory, designated by the `--genomeDir` flag. It will henseforth be referred to as `/path/to/genome`.

**2. Align the Expression Data to the Reference Genome**
```
STAR --runThreadN 24 --runMode alignReads \
--outSAMtype BAM Unsorted SortedByCoordinate \
--genomeDir /path/to/genome \
--outFileNamePrefix /path/to/output/DESIRED_FILE_PREFIX \
--readFilesIn /path/to/lane1/read1,/path/to/lane2/read1 /path/to/lane1/read2,/path/to/lane2/read2
--outFilterType BySJout \
--outSAMattributes NH HI AS NM MD \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--quantMode TranscriptomeSAM \
--alignSJoverhangMin 1 \
--alignSJDBoverhangMin 8 
```
By specifying the option  `--quantMode TranscriptomeSAM`, STAR will output a file `Aligned.toTranscriptome.bam`. This is what we will use for RSEM. For more information regarding these paramters, refer to the [STAR manual](http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf).

If your fastq files are not pair-end reads, you will only have one read to input. Not all illumina results will appear in different 'lane' files, as this largely depends on the size of the inputs and numebr of genes sequenced. Alter this input to represent your data.

**3. RSEM Prepare Genome Reference**
```
/path/to/RSEM/rsem-prepare-reference -p 24 --star \
--gtf /path/to/annotations/gencode.v19.annotation.gtf \
/path/to/genome/fasta/hg19.fa \
/path/to/desired/output/human_ref/hg19
```
The -p flag replaces the `--runThreadN` flag before, and still represnts the number of threads available. This command takes the same inputs as when genrating the genome index with STAR, however greates a unique reference directory for use with RSEM. Ensure that this output directory does not overwrite the directory generated with STAR, as it will be used in the next step... 

**4. RSEM calculate expression**

Before calculating expression, it is important to verify that the input files are valid because we used an alternate aligner (RSEM defaults to Bowtie aligner, this process uses STAR). RSEM requires that the two mates of any paired-end alignements be adjacent. To check this, run the following:
```
/path/to/RSEM/rsem-sam-validator DESIRED_FILE_PREFIXAligned.toTranscriptome.out.bam
```
This could take some time, depending on the size of the bam file. If the command returns `The input file is valid!` then you are good to go! Now for the ultimate step of finally calculating expression...

```
/path/to/RSEM/rsem-calculate-expression --alignments --paired-end -p 24 \
/path/to/STARaligned/DESIRED_FILE_PREFIXAligned.toTranscriptome.out.bam \
/path/to/output/human_ref/hg19 \
/path/to/desired/FPKM/outputs/SAMPLE_NAME
```
If the initial expression data was not paired-end, remove the `--paired-end` flag. The end result will be a file `SAMPLE_NAME.genes.results`. This will contain ensembl gene IDs, FPKM values, along with other extraneous information. 

To isolate the gene IDs and FPKM values, the simplest way is to run the following in a bash environment:
```
cat SAMPLE_NAME.genes.results | awk '{print $1,$7}' > SAMPLE_NAME.fpkm.txt
```

This is the ultimate result, a vector of FPKM expression values for the sample, matched with Ensembl gene IDs. For instructions on appending this vector to a pre-existing Gene Expression Matrix (GEM), continue on!

**5. Downloading and Processing Tissue Data**

This [study](https://www.nature.com/articles/sdata201861) published the initial method for combining GTEx and TCGA data by a common processing method.

![alt text](https://github.com/mrbende/RNAprep/blob/master/sdata201861-f1.jpg "Processing Method")

Here, we download the tissue datasets of interest from the public database. These files present normalized gene expression levels, in FPKM format, after being corrected for batch variations between GTEx and TCGA. These files, divided into independent files by tissue type, can be viewed in [figshare](https://figshare.com/articles/Data_record_3/5330593). Particular datasets of interest can be downloaded independently, or all can be downloaded together.

The rest of this documentation will use the kidney data for example. The following files were then downloaded from the fileshare and extracted from their gunziped format:
```
kidney-rsem-fpkm-gtex.txt
kirp-rsem-fpkm-tcga-t.txt
kirp-rsem-fpkm-tcga.txt
kich-rsem-fpkm-tcga-t.txt
kich-rsem-fpkm-tcga.txt
kirc-rsem-fpkm-tcga-t.txt
kirc-rsem-fpkm-tcga.txt
```
KIRP, KICH, and KIRC are three subtypes of renal cell carcinoma. The tcga data is divided into the tumor samples, denoted as `*-tcga-t.txt`, and healthy tissue samples, denotes as `*-tcga.txt`. We also downloaded the gtex samples, which can now be direcly evaluated with the tcga normal samples. 

The recommended way to use the scripts in this repository is with an Anaconda environment. This is the Anaconda environemnt used by GEMprep later, but we will need some of these packages before then. To create an Anaconda environment:
```bash
module add anaconda3/5.1.0

conda create -n myenv python=3.6 matplotlib mpi4py numpy pandas r scikit-learn seaborn

source deactivate
```
All of the downloaded Gene Expression Matrices (GEMs) should be added to the same directory, containing **only** these files. Then, to activate the conda environment for running the python script...
```
conda activate myenv
```
The `gemmaker.py` script in this repo is used to merge all of the GEMs located in the directory while aligning the genes. It is important to make note of the beginnign of this python script.
```
# This is the absolute path to the directory containing your files
dir_path = '/scratch2/mrbende/KIDNEY/'
# This is the name of the output labels file for tsne plotting
labels_file = 'kidney-labels.txt'
# This is the name of the output GEM
GEM_file = 'kidney_FPKM.txt'
```
These variables should be changed to reflect your environment. The labels file being generated will be relevant in downstream applications so just make note of it here.
