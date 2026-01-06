#Bash script for processing barcode complexity libraries

#Requires bowtie2, bbmap, and samtools to be loaded into path like so:
#export PATH=$PATH:/path/to/bbmap

#Requires a python virtual environment with UMI_tools installed

#Input files include:
#   1) either single end or paired FASTQs containing bit-based barcode diversity libraries
#   2) bowtie2-compatible reference genome with list of 500 possible sequences for each bit

#Output files include
#   1) library_flat.txt, which contains the bit sequence and UMI for each read
#   2) library_completecounts.tsv, which contains the collapsed UMI counts for each barcode

#Last amended by Maddie Urbanek on 1/5/26

#NOTE: this each unique library to be in its own directory--otherwise, it will overwrite any other libraries present

#Arguments
FASTQS=$1 #complete path to FASTQs
PAIRED=$2 #whether the FASTQs are single or paired, n for single and y for paired
BARCODESTART=$3 #base pair number in sequence where the first barcode starts (for Pip-seq4, this should be set to 28)
BOWTIEREFS=$4 #complete path to bowtie2 reference genomes made from bit lists
BBMAP=$5 #path to bbmap scripts

#Example use:
#For single-end
#bash diversity_barcode_alignment.sh /Users/maddieurbanek/Desktop/test/single n 36 /Users/maddieurbanek/Desktop/Rabies/localanalysis/refgenomes/bits /Users/maddieurbanek/Desktop/Rabies/localanalysis/bbmap

#For paired-end
#bash diversity_barcode_alignment.sh /Users/maddieurbanek/Desktop/test/paired y 19 /Users/maddieurbanek/Desktop/Rabies/localanalysis/refgenomes/bits /Users/maddieurbanek/Desktop/Rabies/localanalysis/bbmap


cd $FASTQS

if [ $PAIRED == n ]; then
  echo Single-ended FASTQ
  echo Filtering out low-quality reads
  cp *.fastq.gz barcodes_R1.fq.gz
  bbduk.sh in=barcodes_R1.fq.gz out=temp.fq maq=20 int=f

  echo Aligning first bit
  BITEND=$(($BARCODESTART + 19))
  bbduk.sh in=temp.fq out=temp1.fq ftl=$BARCODESTART ftr=$BITEND
  bowtie2 -x $BOWTIEREFS/bit1 temp1.fq -N 1 -S bit1.sam
  samtools sort bit1.sam -n -o bit1sorted.sam
  awk 'BEGIN {OFS="\t"}; {print $1,$3}' bit1sorted.sam > bit1_hit.txt
  rm temp1.fq
  rm bit1.sam
  rm bit1sorted.sam

  echo Aligning second bit
  BITSTART=$(($BARCODESTART + 21))
  BITEND=$(($BARCODESTART + 40))
  bbduk.sh in=temp.fq out=temp1.fq ftl=$BITSTART ftr=$BITEND
  bowtie2 -x $BOWTIEREFS/bit2 temp1.fq -N 1 -S bit2.sam
  samtools sort bit2.sam -n -o bit2sorted.sam
  awk 'BEGIN {OFS="\t"}; {print $1,$3}' bit2sorted.sam > bit2_hit.txt
  rm temp1.fq
  rm bit2.sam
  rm bit2sorted.sam
  
  echo Aligning second bit
  BITSTART=$(($BARCODESTART + 41))
  BITEND=$(($BARCODESTART + 60))
  bbduk.sh in=temp.fq out=temp1.fq ftl=$BITSTART ftr=$BITEND
  bowtie2 -x $BOWTIEREFS/bit3 temp1.fq -N 1 -S bit3.sam
  samtools sort bit3.sam -n -o bit3sorted.sam
  awk 'BEGIN {OFS="\t"}; {print $1,$3}' bit3sorted.sam > bit3_hit.txt
  rm temp1.fq
  rm bit3.sam
  rm bit3sorted.sam

  echo Pulling UMIs
  bbduk.sh in=barcodes_R1.fq.gz out=temp4.fq ftr=15 qin=33
  $BBMAP/reformat.sh in=temp4.fq out=temp4.sam int=f
  samtools sort temp4.sam -n -o temp5.sam
  awk 'BEGIN {OFS="\t"}; {print $1,$11}' temp5.sam  > cbc_Insert_UMI.txt
fi

if [ $PAIRED == y ]; then
  echo Paired FASTQs
  echo Filtering out low-quality reads
  cp *R2* barcodes_R2.fq.gz
  bbduk.sh in=barcodes_R2.fq.gz out=temp.fq maq=20 int=f

  echo Aligning first bit
  BITEND=$(($BARCODESTART + 19))
  bbduk.sh in=temp.fq out=temp1.fq ftl=$BARCODESTART ftr=$BITEND
  bowtie2 -x $BOWTIEREFS/bit3r temp1.fq -N 1 -S bit1.sam
  samtools sort bit1.sam -n -o bit1sorted.sam
  awk 'BEGIN {OFS="\t"}; {print $1,$3}' bit1sorted.sam > bit1_hit.txt
  rm temp1.fq
  rm bit1.sam
  rm bit1sorted.sam

  echo Aligning second bit
  BITSTART=$(($BARCODESTART + 21))
  BITEND=$(($BARCODESTART + 40))
  bbduk.sh in=temp.fq out=temp1.fq ftl=$BITSTART ftr=$BITEND
  bowtie2 -x $BOWTIEREFS/bit2r temp1.fq -N 1 -S bit2.sam
  samtools sort bit2.sam -n -o bit2sorted.sam
  awk 'BEGIN {OFS="\t"}; {print $1,$3}' bit2sorted.sam > bit2_hit.txt
  rm temp1.fq
  rm bit2.sam
  rm bit2sorted.sam
  
  echo Aligning second bit
  BITSTART=$(($BARCODESTART + 41))
  BITEND=$(($BARCODESTART + 60))
  bbduk.sh in=temp.fq out=temp1.fq ftl=$BITSTART ftr=$BITEND
  bowtie2 -x $BOWTIEREFS/bit1r temp1.fq -N 1 -S bit3.sam
  samtools sort bit3.sam -n -o bit3sorted.sam
  awk 'BEGIN {OFS="\t"}; {print $1,$3}' bit3sorted.sam > bit3_hit.txt
  rm temp1.fq
  rm bit3.sam
  rm bit3sorted.sam

  echo Pulling UMIs
  cp *R1* barcodes_R1.fq.gz
  bbduk.sh in=barcodes_R1.fq.gz out=temp4.fq ftr=15 qin=33
  $BBMAP/reformat.sh in=temp4.fq out=temp4.sam int=f
  samtools sort temp4.sam -n -o temp5.sam
  awk 'BEGIN {OFS="\t"}; {print $1,$11}' temp5.sam  > cbc_Insert_UMI.txt
fi

echo Dropping non-read lines
awk '!/@/' bit1_hit.txt > temp && mv temp bit1_hit.txt
awk '!/@/' bit2_hit.txt > temp && mv temp bit2_hit.txt
awk '!/@/' bit3_hit.txt > temp && mv temp bit3_hit.txt
awk '!/@/' cbc_Insert_UMI.txt > temp && mv temp cbc_Insert_UMI.txt
    
echo Making flat file for input into UMI-tools    
join -1 1 -2 1 -o 1.1,1.2,2.2 bit1_hit.txt bit2_hit.txt > part1.txt
join -1 1 -2 1 -o 1.1,1.2,1.3,2.2 part1.txt bit3_hit.txt > part2.txt
join -1 1 -2 1 -o 1.1,1.2,1.3,1.4,2.2 part2.txt cbc_Insert_UMI.txt > part3.txt
awk 'BEGIN {OFS="\t"}; {print $1"_"$5"_AAAAAAAAAAAA",$2"-"$3"-"$4}' part3.txt > full1.txt
awk '!/\*/' full1.txt > full.txt
sort -k2  full.txt > flat.txt


echo Calculating UMIs per barcode
umi_tools count_tab --per-cell --edit-distance-threshold=1 -I flat.txt -S counts.tsv -L counts.log
tail -n +2 counts.tsv > completecounts.tsv
awk '!/cell/' completecounts.tsv > tmpfile && mv tmpfile completecounts.tsv
{ printf 'CBC\tbarcode\tUMI_Count\n'; cat completecounts.tsv; } > counts.tsv
mv counts.tsv completecounts.tsv

#Clean up
rm *temp*
rm *full*
rm *hit*
rm part*
rm *UMI*
rm *.fq*
