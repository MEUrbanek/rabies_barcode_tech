#Bash script for processing barcode libraries from transcriptome-paired experiments

#Requires bowtie2, bbmap, and samtools to be loaded into path like so:
#export PATH=$PATH:/path/to/bbmap

#Requires a python virtual environment with UMI_tools installed

#Input files include:
#   1) A pair of FASTQs corresponding to Read 1 and Read 2 of a sequenced RVdG barcode library from a transcriptome-paired experiment
#   2) bowtie2-compatible reference genome with list of 500 possible sequences for each bit

#Output files include
#   1) flat.txt, which contains the bit sequence and UMI for each read
#   2) completecounts.tsv, which contains the collapsed UMI counts for each barcode per cell
#   3) helperindex.tsv, which contains collapsed UMI counts for each helper transcript per cell

#Last amended by Maddie Urbanek on 1/5/25

#NOTE: this requires each unique FASTQ pair to be in its own directory--otherwise, it will overwrite any other libraries present
#FASTQs should have also been run through CellRanger to derive a cbc_whitelist.txt file

#Arguments
FASTQS=$1 #complete path to directory with FASTQ pair
CBCWHITELIST=$2 #path to cbc_whitelist.txt output from CellRanger
BARCODESTART=$3 #base pair number in sequence where the first barcode starts (for 10X, this should be set to 28)
BOWTIEREFS=$4 #complete path to bowtie2 reference genomes made from bit lists
BBMAP=$5 #path to bbmap scripts
HELPERSEQUENCE=$6 #Sequence to parse for identifying helpers, for lenti used AGCCATCTGTTGTTT

#Example use:
#bash 11_10X_experiment_barcode_alignment.sh /Users/maddieurbanek/Desktop/test /Users/maddieurbanek/Desktop/test/barcodes 28 /Users/maddieurbanek/Desktop/Rabies/localanalysis/refgenomes/bits /Users/maddieurbanek/Desktop/Rabies/localanalysis/bbmap

#Switch into data directory containing randomer FASTQs
cd $FASTQS


echo Making CBC reference genome
#Make dataset-specific CBC reference genome to collapse any sequencing errors or mutations in cell barcodes into
cp $CBCWHITELIST/barcode_whitelist.txt $FASTQS/barcodes.txt   
awk '1;1' barcodes.txt >cbc.txt
awk 'NR % 2 == 1 {sub(/^/,">")} {print}' cbc.txt >cbc.fa  
bowtie2-build cbc.fa cbc 
rm barcodes.txt  
rm cbc.txt  

#Generate rabies barcode and helper count matrices
#Filter read 2 to remove low-quality barcodes
cp *R2* temp.fq.gz
echo Filtering low-quality barcodes
bbduk.sh in=temp.fq.gz out=temp.fq maq=20 int=f

#Pull bit 1 and align to bit reference genome
echo Aligning first bit

BITEND=$(($BARCODESTART + 19))
bbduk.sh in=temp.fq out=temp1.fq ftl=$BARCODESTART ftr=$BITEND
bowtie2 -x $BOWTIEREFS/bit1 temp1.fq -N 1 -S bit1.sam
samtools sort bit1.sam -n -o bit1sorted.sam
awk 'BEGIN {OFS="\t"}; {print $1,$3}' bit1sorted.sam > bit1_hit.txt
rm temp1.fq
rm bit1.sam
rm bit1sorted.sam

BITSTART=$(($BARCODESTART + 21))
BITEND=$(($BARCODESTART + 40))
#Repeat with bit 2 and bit 3
echo Aligning second bit
bbduk.sh in=temp.fq out=temp1.fq ftl=$BITSTART ftr=$BITEND
bowtie2 -x $BOWTIEREFS/bit2 temp1.fq -N 1 -S bit2.sam
samtools sort bit2.sam -n -o bit2sorted.sam
awk 'BEGIN {OFS="\t"}; {print $1,$3}' bit2sorted.sam > bit2_hit.txt
rm temp1.fq
rm bit2.sam
rm bit2sorted.sam

BITSTART=$(($BARCODESTART + 41))
BITEND=$(($BARCODESTART + 60))
echo Aligning third bit
bbduk.sh in=temp.fq out=temp1.fq ftl=$BITSTART ftr=$BITEND
bowtie2 -x $BOWTIEREFS/bit3 temp1.fq -N 1 -S bit3.sam
samtools sort bit3.sam -n -o bit3sorted.sam
awk 'BEGIN {OFS="\t"}; {print $1,$3}' bit3sorted.sam > bit3_hit.txt
rm temp1.fq
rm bit3.sam
rm bit3sorted.sam

#Align cell barcodes
echo Aligning CBCs
cp *R1* temp.fq.gz
bbduk.sh in=temp.fq.gz out=cbc.fq ftr=15 qin=33
bowtie2 -x ../barcodes/cbc cbc.fq -N 1 -S cbc.sam
samtools sort cbc.sam -n -o cbcsorted.sam
awk 'BEGIN {OFS="\t"}; {print $1,$10}' cbcsorted.sam > cbc.txt

#Pull UMIs
echo Pulling UMIs
bbduk.sh in=temp.fq.gz out=temp4.fq ftr=28 qin=33
$BBMAP/reformat.sh in=temp4.fq out=temp4.sam int=f
samtools sort temp4.sam -n -o temp5.sam
awk 'BEGIN {OFS="\t"}; {print $1,$11}' temp5.sam  > cbc_Insert_UMI.txt

#Remove any lines with missing bits or errors
echo Dropping lines with missing bits
awk '!/@/' bit1_hit.txt > temp && mv temp bit1_hit.txt
awk '!/@/' bit2_hit.txt > temp && mv temp bit2_hit.txt
awk '!/@/' bit3_hit.txt > temp && mv temp bit3_hit.txt
awk '!/@/' cbc.txt > temp && mv temp cbc.txt
awk '!/@/' cbc_Insert_UMI.txt > temp && mv temp cbc_Insert_UMI.txt
    
#Join bit sequences, cell barcode, and UMI
echo Making flat file for input into UMI-tools
join -1 1 -2 1 -o 1.1,1.2,2.2 bit1_hit.txt bit2_hit.txt > part1.txt
join -1 1 -2 1 -o 1.1,1.2,1.3,2.2 part1.txt bit3_hit.txt > part2.txt
join -1 1 -2 1 -o 1.1,1.2,1.3,1.4,2.2 part2.txt cbc.txt > part3.txt
join -1 1 -2 1 -o 1.1,1.2,1.3,1.4,1.5,2.2 part3.txt cbc_Insert_UMI.txt > full.txt

#Organize file for input into UMI_tools
awk 'BEGIN {OFS="\t"}; {print $1"_"$6"_"$5,$2"-"$3"-"$4}' full.txt > full1.txt

#Remove any sequences with missing barcode bits
awk '!/\*/' full1.txt > full.txt
sort -k2  full.txt > flat.txt

#Run UMI_tools on input file, calculating UMIs (collapsed by 1 mismatch) per barcode per cell
echo Running UMI-tools
umi_tools count_tab --per-cell --edit-distance-threshold=1 -I flat.txt -S counts.tsv -L counts.log
tail -n +2 counts.tsv > completecounts.tsv
awk '!/cell/' completecounts.tsv > tmpfile && mv tmpfile completecounts.tsv
{ printf 'CBC\tbarcode\tUMI_Count\n'; cat completecounts.tsv; } > counts.tsv

#Output completecounts.tsv file
mv counts.tsv completecounts.tsv
echo Done generating barcode count matrix!

#Process helper counts
echo Pulling helper transcripts
cp *R2* temp.fq.gz
bbduk.sh in=temp.fq.gz out=temp.fq maq=20 int=f
bbduk.sh in=temp.fq literal=$HELPERSEQUENCE k=12 skipr1 restrictright=125 hdist=1 out=helper.fq outm1=helper_R1.fq outm2=helper_R2.fq int=t
$BBMAP/reformat.sh in=helper_R2.fq out=helper_R2.sam int=f
samtools sort helper_R2.sam -n -o helper_R2sorted.sam
awk 'BEGIN {OFS="\t"}; {print $1,$11}' helper_R2sorted.sam  > helper.txt

awk '!/@/' helper.txt > temp && mv temp helperhit.txt  
join -1 1 -2 1 -o 1.1,1.2,2.2 cbc.txt helperhit.txt > helperindex.txt
join -1 1 -2 1 -o 1.1,1.2,1.3,2.2 helperindex.txt cbc_Insert_UMI.txt > helperindexfull.txt
awk 'BEGIN {OFS="\t"}; {print $1"_"$4"_"$2,$3}' helperindexfull.txt > helperindex.txt
awk 'BEGIN {OFS="\t"}; {print $1,"helper"}' helperindex.txt > helperUMI.txt
umi_tools count_tab --per-cell --edit-distance-threshold=1 -I helperUMI.txt -S helperindex.tsv -L counts.log
awk '!/cell/' helperindex.tsv > tmpfile && mv tmpfile helperindex.tsv
tail -n +2 helperindex.tsv > helpercounts.tsv
{ printf 'CBC\thelper\tUMI_Count\n'; cat helpercounts.tsv; } > helperindex.tsv
echo Done generating helper count matrix!

#Clean up
rm *temp*
rm *full*
rm *.sam
rm *hit*
rm part*
rm cbc*
rm *UMI*
rm *.fq
mv helperindex.tsv xhelperindex.tsv
rm helper*
mv xhelperindex.tsv helperindex.tsv
