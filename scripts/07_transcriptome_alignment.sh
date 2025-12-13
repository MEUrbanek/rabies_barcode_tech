#Bash script for processing transcriptome FASTQs using Pipseeker
#Requires pipseeker and accompanying reference genomes
    #For Shin & Urbanek, used Pipseeker version 
    #with human reference genome: pipseeker-gex-reference-GRCh38-2022.04

#Input files include:
#   1) paired FASTQs containing transcriptome libraries

#Output files include:
#   1) filtered gene count matrices at 5 different sensitivities

#Note: this is compatible with Pip-seq V4 chemistry--other chemistries may need a different version of Pipseeker

#Arguments
PIPSEEKER=$1 #complete path to pipseeker software
FASTQ=$2 #complete path to paired FASTQ files and either the FASTQ pair's prefix or . for all FASTQs in the directory
REFGENOME=$3 #complete path to pipseeker compatible reference genome

#Example use
#bash transcriptome_alignment.sh ~/pipseeker/pipseeker  ~/data/fastqs/c1/txn/. ~/pipseeker/pipseeker-gex-reference-GRCh38-2022.04

#Set ulimit to 10000 to keep STAR from halting
echo Increasing open file limit to 10,000
ulimit -n 10000

#Execute pipseeker
$PIPSEEKER full --fastq $FASTQ --star-index-path $REFGENOME --chemistry v4 --output $FASTQ --retain-barcoded-fastqs --skip-version-check --skip-preflight-check
