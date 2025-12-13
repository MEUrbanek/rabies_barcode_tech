#Bash script for processing barcode diversity libraries

#Requires bowtie2, bbmap, and samtools to be loaded into path like so:
#export PATH=$PATH:/path/to/bbmap

#Requires a python virtual environment with UMI_tools installed

#Input files include:
#   1) either single end or paired FASTQs containing bit-based barcode diversity libraries
#   2) bowtie2-compatible reference genome with list of 500 possible sequences for each bit

#Output files include
#   1) library_flat.txt, which contains the bit sequence and UMI for each read
#   2) library_completecounts.tsv, which contains the collapsed UMI counts for each barcode

#Last amended by Maddie Urbanek on !!!!!

#NOTE: this each unique library to be in its own directory--otherwise, it will overwrite any other libraries present

#Arguments
FASTQS=$1 #complete path to directory with either single barcode FASTQ or FASTQ pair
BOWTIEREFS=$2 #complete path to bowtie2 reference genomes made from bit lists
BARCODESTART=$3 #base pair number in sequence where the first barcode starts
OUTPUT=$4 #prefix to add to completecounts.tsv file to differentiate libraries

#Example use:
#bash diversity_barcode_alignment.sh 

#Switch into data directory containing randomer FASTQs
cd $FASTQS

--genomeDir $BOWTIEREFS
--readFilesIn "${i}_R1_001.fastq"

#FOR PAIRED END READS
#Filter low quality R2 with BBDUK
bbduk.sh in=*_R2* out=temp.fq maq=20 int=f

#Trim to bit 1, align, and pull sequence
bbduk.sh in=temp.fq out=temp1.fq ftl=$BARCODESTART ftr=$BARCODESTART+19
#For single end reads, the rabies genome is flipped from the index list, so we use the reverse bowtie reference indices in flipped order
bowtie2 -x $BOWTIEREFS/bit3r temp1.fq -N 1 -S bit1.sam
samtools sort bit1.sam -n -o bit1sorted.sam
awk 'BEGIN {OFS="\t"}; {print $1,$3}' bit1sorted.sam > bit1_hit.txt
rm temp1.fq                                   
rm bit1.sam                                                                                                 
rm bit1sorted.sam                          

bbduk.sh in=temp.fq out=temp1.fq ftl=$BARCODESTART+20 ftr=$BARCODESTART+39
bowtie2 -x $BOWTIEREFS/bit2r temp1.fq -N 1 -S bit2.sam
samtools sort bit2.sam -n -o bit2sorted.sam
awk 'BEGIN {OFS="\t"}; {print $1,$3}' bit2sorted.sam > bit2_hit.txt
rm temp1.fq                                              
rm bit2.sam                                      
rm bit2sorted.sam                        

bbduk.sh in=temp.fq out=temp1.fq ftl=$BARCODESTART+40 ftr=$BARCODESTART+59
bowtie2 -x $BOWTIEREFS/bit1r temp1.fq -N 1 -S bit3.sam
samtools sort bit3.sam -n -o bit3sorted.sam                
awk 'BEGIN {OFS="\t"}; {print $1,$3}' bit3sorted.sam > bit3_hit.txt                                
rm temp1.fq                            
rm bit3.sam                                                           
rm bit3sorted.sam

#Pull UMI from read 1, which should be the first 28 base pairs of that read
bbduk.sh in=*_R1* out=temp4.fq ftr=28 qin=33
/Users/maddieurbanek/Desktop/Rabies/localanalysis/bbmap/reformat.sh in=temp4.fq out=temp4.sam int=f
samtools sort temp4.sam -n -o temp5.sam
awk 'BEGIN {OFS="\t"}; {print $1,$11}' temp5.sam  > cbc_Insert_UMI.txt
    
#Drop any non-mapping lines
awk '!/@/' bit1_hit.txt > temp && mv temp bit1_hit.txt
awk '!/@/' bit2_hit.txt > temp && mv temp bit2_hit.txt
awk '!/@/' bit3_hit.txt > temp && mv temp bit3_hit.txt
awk '!/@/' cbc_Insert_UMI.txt > temp && mv temp cbc_Insert_UMI.txt

#Join text files for input into UMI_tools    
join -1 1 -2 1 -o 1.1,1.2,2.2 bit1_hit.txt bit2_hit.txt > part1.txt
join -1 1 -2 1 -o 1.1,1.2,1.3,2.2 part1.txt bit3_hit.txt > part2.txt
join -1 1 -2 1 -o 1.1,1.2,1.3,1.4,2.2 part2.txt cbc_Insert_UMI.txt > part3.txt

#Print fake cell barcode so UMI_tools treats all barcodes as coming from one entry
awk 'BEGIN {OFS="\t"}; {print $1"_"$5"_AAAAAAAAAAAA",$2"-"$3"-"$4}' part3.txt > full1.txt

#Drop any lines with missing bits or UMIs
awk '!/\*/' full1.txt > full.txt

#Sort UMI_tools input file 
sort -k2  full.txt > flat.txt

#Run UMI_tools for whole library, collapsing UMIs with 1 or fewer UMI mismatches
umi_tools count_tab --per-cell --edit-distance-threshold=1 -I flat.txt -S counts.tsv -L counts.log

#Drop non-barcode lines at bottom of output file
tail -n +2 counts.tsv > completecounts.tsv

#Fix headings
awk '!/cell/' completecounts.tsv > tmpfile && mv tmpfile completecounts.tsv
{ printf 'CBC\tbarcode\tUMI_Count\n'; cat completecounts.tsv; } > counts.tsv

#Rename final file to "output"_completecounts.tsv
mv counts.tsv "$BOWTIEREFS_completecounts.tsv"
