#Bash script for processing randomer FASTQs
#Requires bbduk.sh, samtools, and UMI_tools
#Input files include:
#   1) paired FASTQs containing randomer barcode sequence and UMI information
#Output files include
#   1) randomer_flat.txt, which contains the exact randomer and UMI sequence for each read (which can be used to collapse similar barcodes)
#   2) randomer_completecounts.tsv, which contains the collapsed UMI counts for each randomer sequence (without any collapsing of similar randomer sequences)
#Last amended by Maddie Urbanek on 12/12/2025

#Note: the UMI for each transcript was encoded on the reverse transcription primer and replaced the i5 index. As a result, the sequencing run was split into four FASTQ files depending on the first letter of the UMI so it could be sequenced in tandem with a paired end library.

#Arguments
FASTQS=$1 #complete path to directory with either single barcode FASTQ or FASTQ pair

#Be sure that bbmap and UMI-tools are available through path:
#export PATH=/Users/maddieurbanek/Desktop/Rabies/localanalysis/myenv:$PATH

#Example use:
#bash 02_randomer_processing.sh /Users/maddieurbanek/Desktop/revision_data/resubmission/data/fastqs/diversity_libs/randomer/ 

#Switch into data directory containing randomer FASTQs
cd $FASTQS

#Concatenate separate randomer FASTQs together into one pair of FASTQs
cat *R1_001.fastq.gz > barcode_R1.fq.gz
cat *R2_001.fastq.gz > barcode_R2.fq.gz

#gunzip read 2 FASTQ file to manually check for UMI and randomer sites
gunzip barcode_R2.fq.gz

#Randomer starts 29 base pairs in and is 30 base pairs long, UMI starts 92 base pairs in and is 15bp long
#Eliminate low quality reads with bbduk.sh
bbduk.sh in=barcode_R2.fq out=temp.fq maq=20 int=f

#To reduce contaminaton, bbduk.sh's literal function can be used to pull all sequences that have approximately the same sequence upstream of the randomer barcodes that matches our plasmid backbone:
bbduk.sh in=temp.fq literal=ACGGCAATTAGGTA k=6 restrictright=20 hdist=2 out=garbage.fq outm1=randomer_R2.fq
#And remove the non-matching sequences
rm garbage.fq
rm temp.fq

#randomer_R2.fq can be checked with a program like FASTQC to ensure base pairs at start of the sequence match!

#Trim randomer_R2.fq file down to the randomer sequence
bbduk.sh in=randomer_R2.fq out=temp1.fq ftl=29 ftr=58
#Reformat randomer for input into flat.txt file 
/Users/maddieurbanek/Desktop/Rabies/localanalysis/bbmap/reformat.sh in=temp1.fq out=temp1.sam int=f
samtools sort temp1.sam -n -o temp5.sam
awk 'BEGIN {OFS="\t"}; {print $1,$11}' temp5.sam  > randomer.txt
rm temp1*
rm temp5*

#Pull UMI from read 2 file
bbduk.sh in=randomer_R2.fq out=temp1.fq ftl=92 ftr=106
#Reformat UMI for input into flat.txt file
/Users/maddieurbanek/Desktop/Rabies/localanalysis/bbmap/reformat.sh in=temp1.fq out=temp1.sam int=f
samtools sort temp1.sam -n -o temp5.sam
awk 'BEGIN {OFS="\t"}; {print $1,$11}' temp5.sam  > cbc_Insert_UMI.txt
rm temp1*
rm temp5*

#Remove any missing sequences from .txt files
awk '!/@/' randomer.txt > temp && mv temp randomer.txt
awk '!/@/' cbc_Insert_UMI.txt > temp && mv temp cbc_Insert_UMI.txt

#Join randomer and UMI into one text file
join -1 1 -2 1 -o 1.1,1.2,2.2 randomer.txt cbc_Insert_UMI.txt > part1.txt

#Print UMI_tools compatible cell barcode sequence into flat.txt file
awk 'BEGIN {OFS="\t"}; {print $1"_"$3"_AAAAAAAAAAAA",$2}' part1.txt > full1.txt

#Remove line missing components after joining
awk '!/\*/' full1.txt > full.txt

#Sort by read ID
sort -k2  full.txt > randomer_flat.txt

#SHOULDN'T NEED ANY OF THE FOLLOWING FOR COLLAPSING
#Input flat file into UMI_tools to count UMIs per unique barcode sequence with a hamming distance on UMIs set to 1
#umi_tools count_tab --per-cell --edit-distance-threshold=1 -I randomer_flat.txt -S counts.tsv -L counts.log

#Drop file suffix
#tail -n +2 counts.tsv > completecounts.tsv

#Add new header describing column outputs
#awk '!/cell/' completecounts.tsv > tmpfile && mv tmpfile completecounts.tsv
#{ printf 'CBC\tbarcode\tUMI_Count\n'; cat completecounts.tsv; } > counts.tsv

#Send file to completecounts.tsv
#mv counts.tsv randomer_completecounts.tsv

#Remove intermediate files
rm randomer.txt
rm cbc_Insert_UMI.txt
rm full*
rm part*

rm randomer_R2.fq
rm barcode_R1.fq.gz
rm barcode_R2.fq

#rm completecounts.tsv
