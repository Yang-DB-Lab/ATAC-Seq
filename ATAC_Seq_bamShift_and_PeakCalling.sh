#!/bin/bash


##***************** Description of this script *****************##
## This script is adapted from 
## https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/
## It is used to analyze ATAC-Seq data from fastq files
##
##**************************************************************##


## Frist get into the folder contains cleaned bam file of ATAC-Seq ##
sample_folder='/data/Exercise/ATAC_Seq/SRR891269'
cd $sample_folder


## ***************** How to use this script ****************** ##
## This script is intended to analyze one sample(includes paired
## reads)
## To work with multiple samples, use for loop of bash to run
## this script on different samples.
##
## The script should be used in this way:
## bash ATAC_Seq_Pipeline.sh <sample_name>
## One cleaned bam file(without mitochondrial reads, PCR duplicates removed, etc)
## get from "ATAC_Seq_fastq_to_bam.sh" named $sample_name.cleaned.high.quality.bam
## will be used for following analysis.
## ********************************************************** ##


# sample_name='SRR891269'
sample_name=$1
## Then check conda env list. 
## alignmentSieve tool from deepTools will be used to shift ATAC-Seq bam file.

## conda env list
### conda environments:
###
### base                  *  /home/guang/anaconda3
### NGS_Py2.7                /home/guang/anaconda3/envs/NGS_Py2.7


## To use alignmentSeive too. First activate the conda environment containing
## deepTools
conda activate NGS_Py2.7

## ************* Shift ATAC-Seq bam file *************** ##
# index the cleaned bam file of ATAC-Seq first
samtools index -@ 6 $sample_name.cleaned.high.quality.bam

# use --ATACshift to do bam file read shift
alignmentSieve --numberOfProcessors 6 --ATACshift --bam $sample_name.cleaned.high.quality.bam -o $sample_name.shifted.tmp.bam

# the bam file needs to be sorted again
samtools sort -@ 6 -O bam -o $sample_name.shifted.bam $sample_name.shifted.tmp.bam
samtools index -@ 6 $sample_name.shifted.bam
rm $sample_name.shifted.tmp.bam
## **************************************************** ##

## ************* Peak calling using MACS2 ******************* ##

# -t: The IP data file (this is the only REQUIRED parameter for MACS)!! The REAL experimental data.
# -c: The control or mock data file!! Data from control group.
# -f: format of input file; Default is “AUTO” which will allow MACS to decide the format automatically.
# -g: mappable genome size which is defined as the genome size which can be sequenced; some precompiled values provided.


# -f BAMPE, use paired-end information
# --keep-dup all, keep all duplicate reads.
# --outdir: MACS2 will save all output files into speficied folder for this option
# -n: The prefix string for output files
# -B/--bdg: store the fragment pileup, control lambda, -log10pvalue and -log10qvalue scores in bedGraph files

macs2 callpeak -f BAMPE -g hs --keep-dup all --cutoff-analysis -n $sample_name \
  -t $sample_name.shifted.bam --outdir macs2/$sample_name 2> macs2.of.$sample_name.log



## ************ Visualization of bam file **************** ##

## ***** Bonus: Change chromosome notation (re-head) of bam file  **** ##

### Check chromosome names used in the bam file
### samtools idxstats $sample_name.shifted.bam | head -10
### 1	248956422	1854622	0
### 10	133797422	1045792	0
### 11	135086622	1073272	0
### 12	133275309	1041088	0
### 13	114364328	683574	0
### 14	107043718	581912	0
### 15	101991189	475290	0
### 16	90338345	587082	0
### 17	83257441	503336	0
### 18	80373285	557980	0
############### IMPORTANT NOTE: ############################
## Chromosome names in ENSEMBL genome are noted as: 1, 2, 3 ...
## while chorosome name used by UCSC are: Chr1, Chr2, Chr3 ...
## if one file created using ENSEMBL notation is checked by 
## tool using UCSC notation(like IGV).
## You will get nothing!!! Because IGV can never find any chromosome
## names(Chr1, Chr2 ...)

## To use UCSC realted tools, the chromosome names in bam file 
## generated using ENSEML reference genome must be transformed
## into UCSC format.

## The head of bam file looks like this:
## @HD	VN:1.6	SO:coordinate
## @SQ	SN:1	LN:248956422
## @SQ	SN:10	LN:133797422
## @SQ	SN:11	LN:135086622
## @SQ	SN:12	LN:133275309
## Need to change "SN:<num>" to "SN:chr<num>" to match UCSC notation

## Explanation of the piped script below:
## Get the old header(chromosome notation) of bam file
## Then using sed, change the bam header:
## from SN:[0-9XY] to SN:chr[0-9XY]
## from SN:MT to SN:chrM.
## Finally, use reheader function of samtools to re-assign the changed
## header to original bam file to produce a UCSC-headed new bam file.
samtools view -H $sample_name.shifted.bam \
| sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' \
| samtools reheader - $sample_name.shifted.bam > $sample_name.shifted.UCSC.header.bam
## index the reheaded bam file
samtools index -@ 6 $sample_name.shifted.UCSC.header.bam

## ***** Now, this bam chromosome names change from ENSEMBL to UCSC ***** ##
## NOTE: IGV will automatically deal with different notation.
##       If nothing appeared after open a bam file(which you know has many reads),
##       kepp zooming until you see short reads in the file.
## ********************************************************************* ##


# create bigWig files for visualizing using bamCoverage in deepTools





## ************* ATAC-Seq data quality check ************** ##





