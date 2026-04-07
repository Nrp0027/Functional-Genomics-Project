#!/bin/bash

######### FunGen Course Instructions ############
## Purpose: The purpose of this script is to 
## 	learn to make scratch directory
## 	learn to define variables
## 	download data from NCBI SRA using the SRAtoolkit and the SRA run IDs: https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
## 	use FASTQC to evaluate the quality of the data: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## Download from SRA: Input Data: NA
## 			Output: Downloaded read files, R1 and R2 files for each sample if paired-end data (FASTQ)
## FASTQC 	InPut: Downloade SRA files .fastq
##		Output: is a folder for each file that contains a .html file to visualize the quality, and .txt files of quality statistics.
##			The last line of this script will make a tarball of the output directory to bring back to your computer
## 	After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
## 	then run the script by using "run_script [script name]"
## 	suggested paramenters are below to submit this script.
## 		queue: class
##		core: 1
##		time limit (HH:MM:SS): 04:00:00 
##		Memory: 4gb
##		run on asax
###############################################


########## Load Modules
source /apps/profiles/modules_asax.sh.dyn
module load sra
module load fastqc/0.10.1

##########  Define variables and make directories
## Replace the numbers in the brackets with Your specific information
  ## make variable for your ASC ID so the directories are automatically made in YOUR directory
  ## These are represented in the code by [#] replace these according to the examples provided
MyID=[1]          ## Example: MyID=aubrmg001

  ## Make variable that represent YOUR working directory(WD) in scratch, your Raw data directory (DD) and the pre or postcleaned status (CS).
DD=/scracth/${MyID}/FinalProject/RawData   			## Example: DD=/scratch/${MyID}/PracticeRNAseq/RawData
WD=/scratch/${MyID}/FinalProject				## Example: WD=/scratch/${MyID}/PracticeRNAseq
RDQ=RawDataQuality
 
##  make the directories in SCRATCH for holding the raw data 
## -p tells it to make any upper level directories that are not there. Notice how this will also make the WD.
mkdir -p ${DD}
## move to the Data Directory
cd ${DD}

##########  Download data files from NCBI: SRA using the Run IDs
  ### from SRA use the SRA tool kit - see NCBI website https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
	## this downloads the SRA file and converts to fastq
	## -F 	Defline contains only original sequence name.
	## -I 	Append read id after spot id as 'accession.spot.readid' on defline.
	## splits the files into R1 and R2 (forward reads, reverse reads)

## These samples are from Bioproject PRJNA437447. An experiment on Daphnia pulex, 5 samples on ad lib feed, 5 samples on caloric restriction diet
## https://www.ncbi.nlm.nih.gov/bioproject?LinkName=sra_bioproject&from_uid=5206312
## For class only do the 6 that you are assigned, delete the other 4 from this list

## We must run vdb-config to force creation of the default config file, otherwise we will get an error. This is a 'hack'. 

vdb-config --interactive > /dev/null 2>&1 <<EOF
q
EOF

fastq-dump -F --split-files SRR22459445
fastq-dump -F --split-files SRR22459088


##### Extra ####
## If you are downloading data from a sequencing company instead of NCBI, using wget for example, then calculate the md5sum values of all the files in the folder (./*), and read into a text file.
## then you can compare the values in this file with the ones provided by the company.
#md5sum ./* > md5sum.txt

##### Extra ####
## If you data comes with multiple R1 and R2 files per individual. You can contatenate them together using "cat" before running FASTQC
## see examples below for one file. You will probably want to use a loop to process through all the files.
#cat SRR22459445*_R1_*.fastq.gz > SRR6819014_All_R1.fastq.gz
#cat SRR22459445*_R2_*.fastq.gz > SRR6819014_All_R2.fastq.gz


############## FASTQC to assess quality of the sequence data
## FastQC: run on each of the data files that have 'All' to check the quality of the data
## The output from this analysis is a folder of results and a zipped file of results and a .html file for each sample
## Note, if the parent directory already exists using -p will result in an error.

mkdir ${WD}/${RDQ}
fastqc *.fastq --outdir=${WD}/${RDQ}

#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
cd ${WD}/${RDQ}
tar cvzf ${RDQ}.tar.gz  ${WD}/${RDQ}/*
## when finished use scp or rsync to bring the tarballed .gz results file to your computer and open the .html file to evaluate the quality of your raw data.
