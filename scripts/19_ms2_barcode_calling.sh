#Bash script for processing unbarcoded rabies MS2 comparisons

#Requires bowtie2, bbmap, and samtools to be loaded into path like so:
#export PATH=$PATH:/path/to/bbmap

#Requires a python virtual environment with UMI_tools installed

#Input files include:
#   1) Barcoded rabies FASTQs for libraries from MS2 experiments

#Output files include
#   1) library_completecounts.tsv, which contains the collapsed UMI counts for each barcode per cell

#Last amended by Maddie Urbanek on 12/12/25

#NOTE: this each unique library to be in its own directory--otherwise, it will overwrite any other libraries present

#Arguments
FASTQS=$1 #complete path to directory with either single barcode FASTQ or FASTQ pair
PIPSEEKER=$2 #path to pipseeker

#Be sure that bbmap and UMI-tools are available through path:
#export PATH=/Users/maddieurbanek/Desktop/Rabies/localanalysis/myenv:$PATH
#export PATH=$PATH:/Users/maddieurbanek/Desktop/Rabies/localanalysis/bbmap/

#Example use:
#bash ms2_barcode_calling.sh /Users/maddieurbanek/Desktop/revision_data/resubmission/data/fastqs/slice_experiments/ms2_comparison/ms2/bc/ms2_wpre/ /Users/maddieurbanek/pipseeker/pipseeker

cd $FASTQS

#Reformat raw FASTQs for 10X-compatible CBC

$PIPSEEKER barcode --fastq . --chemistry v4 --output . --skip-version-check

cd barcoded_fastqs
#gunzip *.gz
cat barcoded_*_R1.fastq.gz > barcode_R1.fq.gz
cat barcoded_*_R2.fastq.gz > barcode_R2.fq.gz
gunzip *.gz
rm *fastq*

##Make dataset-specific CBC reference genome
cd ../barcodes
awk '1;1' barcode_whitelist.txt > cbc.txt
awk 'NR % 2 == 1 {sub(/^/,">")} {print}' cbc.txt >cbc.fa
bowtie2-build cbc.fa cbc
cd ../barcoded_fastqs

grep -B 1 CGGTAGCTTTTCAGTCGAG barcode_R2.fq > perfect_hit.txt 

bbduk.sh in=barcode_R1.fq out=cbc.fq ftr=15 qin=33
bowtie2 -x ../barcodes/cbc cbc.fq -N 1 -S cbc.sam
samtools sort cbc.sam -n -o cbcsorted.sam
awk 'BEGIN {OFS="\t"}; {print $1,$10}' cbcsorted.sam > cbc.txt

bbduk.sh in=barcode_R1.fq out=temp4.fq ftr=28 qin=33
/Users/maddieurbanek/Desktop/Rabies/localanalysis/bbmap/reformat.sh in=temp4.fq out=temp4.sam int=f
samtools sort temp4.sam -n -o temp5.sam
awk 'BEGIN {OFS="\t"}; {print $1,$11}' temp5.sam  > cbc_Insert_UMI.txt

awk '!/@/' cbc.txt > temp && mv temp cbc.txt
awk '!/@/' cbc_Insert_UMI.txt > temp && mv temp cbc_Insert_UMI.txt
    
join -1 1 -2 1 -o 1.1,1.2,2.2 cbc.txt cbc_Insert_UMI.txt > part1.txt
awk 'BEGIN {OFS="\t"}; {print $1"_"$3"_"$2,"barcode"}' part1.txt > full1.txt
awk '!/\*/' full1.txt > full.txt
sort -k2  full.txt > flat.txt

umi_tools count_tab --per-cell --edit-distance-threshold=1 -I flat.txt -S counts.tsv -L counts.log
tail -n +2 counts.tsv > completecounts.tsv
awk '!/cell/' completecounts.tsv > tmpfile && mv tmpfile completecounts.tsv
{ printf 'CBC\tbarcode\tUMI_Count\n'; cat completecounts.tsv; } > counts.tsv
mv counts.tsv completecounts.tsv

rm *hit*
rm *sam
rm *full*
rm *cbc*
rm counts.tsv
