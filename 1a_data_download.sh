!/bin/sh -x
set -e
set -u
set -o pipefail

# ----------------------------------------------------------------- #
#	Download of FASTQ-files from the European Nucleotide Archive.	#
#	Download of reference genome from UCSC							#
# ----------------------------------------------------------------  #

# PART 1: DOWNLOAD OF FASTQ-FILES

# The sequencing datasets of the publication from Lugowska et al. (2018) are available in the European Nucleotide Archive repository under the accession number PRJEB20222:
# https://www.ebi.ac.uk/ena/data/view/PRJEB20222
# I wanted to start my own analysis with the FASTQ files of this project. At first, I wanted to create a short script to automatically download all FASTQ-files from the site.
# In order to create the script, I first downloaded an overview table. For this, I opened the "select columns" option which is just above the main table on the website.
# I selected only the "Fastq md5" and "FASTQ files (FTP)" columns. Using the "TEXT" button situated just above, I downloaded the overview table.

# Using this table, I created a download script called "1b_fastq_download.sh":
touch 1b_fastq_download.sh # Create empty file to file with script
grep -v "^fastq" PRJEB20222.txt | awk '$1="wget"'  > 1b_fastq_download.sh # Exclude first row (which is the header) and add wget to each line. Write all rows to empty file.
chmod u+x 1b_fastq_download.sh # Make script executable

# Running the download script (This was done on 2018-04-25):
./1b_fastq_download.sh # Download of the FASTQ-files (this will take a while as the total size is approximately 2 GB)
mkdir data/fastq # Create directory for the FASTQ-files
mv ERR* data/fastq/ # transfer all downloaded FASTQ-files to the new folder

# Check if the downloaded FASTQ-files are not corrupted:
grep -v "^fastq" PRJEB20222.txt | sed -e 's!ftp.sra.ebi.ac.uk/vol1/fastq/ERR193/.*/.*/!!' > md5sum.txt # Create separate file with all md5sums
mv md5sum.txt data/fastq/md5sum.txt # Move md5sum.txt file to folder with FASTQ-files
cd data/fastq # Move to folder with FASTQ-files
md5sum -c md5sum.txt # Check file integrity
# ERR1934147.fastq.gz: OK
# ERR1934148.fastq.gz: OK
# ERR1934149.fastq.gz: OK
# ERR1934150.fastq.gz: OK
# ERR1934151.fastq.gz: OK
# ERR1934152.fastq.gz: OK
# --> all were OK
rm md5sum.txt
cd ../..


# PART 2: DOWNLOAD OF HUMANE REFERENCE GENOME

# Download entire human genome (hg38):
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz # downloaded 2018-05-22 | This will download approximately 984 MB of data
mkdir data/genome # Create folder for the genome file
mv hg38.fa.gz data/genome/hg38.fa.gz # move file

# Check if the downloaded genomic file is not corrupted:
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt # Download md5sum information file
cat md5sum.txt | grep hg38.fa.gz > md5sum_gz.txt # Extract md5sum value for hg38.fa.gz
rm md5sum.txt
mv md5sum_gz.txt data/genome/md5sum_gz.txt
cd data/genome
md5sum -c md5sum_gz.txt # Check file integrity
# hg38.fa.gz: OK
rm md5sum_gz.txt
cd ../..


# PART 3: DOWNLOAD OF EXAMPLE BAM FILE AND ITS INDEX (also available on the PRJEB20222 project site of the European Nucleotide Archive)
cd data
# Download bam:
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA885/ERA885352/bam/IonXpress_001_R_2015_05_06_19_09_56_user_SN2-55-Ion_AmpliSeq_Cancer_Hotspot_Panel_v2_Auto_user_SN2-55-Ion_AmpliSeq_Cancer_Hotspot_Panel_v2_196.bam
# Download index:
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA885/ERA885352/bam/IonXpress_001_R_2015_05_06_19_09_56_user_SN2-55-Ion_AmpliSeq_Cancer_Hotspot_Panel_v2_Auto_user_SN2-55-Ion_AmpliSeq_Cancer_Hotspot_Panel_v2_196.bam.bai
