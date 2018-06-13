!/bin/sh -x
set -e
set -u
set -o pipefail

# ----------------------------------------------------------------- #
#		Map FASTQ-files using the Burrows-Wheeler Aligner			#
# ----------------------------------------------------------------- #


# PART 1: EXTRACT TARGET SEQUENCES FROM THE HUMAN GENOME
# The "Ion AmpliSeq Cancer Hotspot Panel v2" comes with a Browser Extensible Data file (called CHP2.20131001.hotspots.bed) which indicates the DNA positions to which the FASTQ-files shall be matched. I don't have this BED file and I couldn't find it on the internet. So I will create a BED-file from the example BAM file that I downloaded from the European Nucleotide Archive.

# Determine against which build the read sequences of the BAM file had been mapped:
samtools view -H data/IonXpress_001_R_2015_05_06_19_09_56_user_SN2-55-Ion_AmpliSeq_Cancer_Hotspot_Panel_v2_Auto_user_SN2-55-Ion_AmpliSeq_Cancer_Hotspot_Panel_v2_196.bam | grep "tmap.*mapall"
# --> @PG	ID:tmap	CL:mapall -n 12 -f /results/referenceLibrary/tmap-f3/hg19/hg19.fasta -r basecaller_results/IonXpress_001_rawlib.basecaller.bam -v -Y -u --prefix-exclude 5 -o 2 stage1 map4	VN:4.4.8 (e29ce66) (201502092025)
# --> The "tmap mapall" aligner used hg19 as the reference genome

# Extract the BED-file (estimated on starting position for each run + 500 bp):
samtools view data/IonXpress_001_R_2015_05_06_19_09_56_user_SN2-55-Ion_AmpliSeq_Cancer_Hotspot_Panel_v2_Auto_user_SN2-55-Ion_AmpliSeq_Cancer_Hotspot_Panel_v2_196.bam | awk '{print $3 "\t" $4 "\t" $4+500}' > data/token1.bed
# Keep only unique rows (exclude all rows starting with "*"), remove duplicates:
grep -v "^\*" data/token1.bed | sort -u > data/token2.bed
# sort according to chromosome and then position:
sort -k1,1V -k2,2n data/token2.bed > data/token3.bed
# merge overlapping entries in bed file:
mergeBed -i data/token3.bed > data/panel_hg19.bed
rm data/token1.bed data/token2.bed data/token3.bed

# Convert BED-file to hg38 genome build:
./liftOver data/panel_hg19.bed data/hg19ToHg38.over.chain data/panel_hg38.bed unMapped
# sort resulting hg38 bed file:
sort -k1,1V -k2,2n data/panel_hg38.bed > data/panel_hg38-2.bed
# extract genome coordinates from genome file index:
awk '{print $1"\t"$2}'head data/genome/hg38.fa.fai > data/genome/genome_coordinates.txt
# sort genome coordinates file according to the same rules as the bed file:
sort -k1,1V -k2,2n data/genome/genome_coordinates.txt > data/genome/genome_coordinates-2.txt
# create bed file with complement sections for bed file:
bedtools complement -i data/panel_hg38-2.bed -g data/genome/genome_coordinates-2.txt > data/complement.bed
# create genomic sequence where only the panel regions are visible:
bedtools maskfasta -fi data/genome/hg38.fa -bed data/complement.bed -fo data/genome/hg38_masked.fa

# Extract sequences defined by panel_hg38 from genome:
# bedtools getfasta -fi data/genome/hg38.fa -bed data/panel_hg38.bed -fo data/genome/hg38_panel.fa


# PART2: TMAP-MAPPING
# annotate genome with TMAP:
# ./TMAP/tmap index -vf data/genome/hg38_panel.fa
./TMAP/tmap index -vf data/genome/hg38_masked.fa
# unzip all trimmed fastq-files (as tmap mapall only accects decompressed fastq files):
ls data/fastq_trimmed/* | xargs -n1 bgzip -d
# create folder for results:
mkdir data/bam
# The for-loop will map the files, convert them from SAM to BAM and sort them:
ls data/fastq_trimmed/* > fastq_names.txt # create txt-file with fastq filenames
FASTQ_FILES=($(cat fastq_names.txt)) # load names into bash array
echo ${FASTQ_FILES[@]} # just for verification
for fastq_file in ${FASTQ_FILES[@]}
do
	RESULTS_FILE="$(basename $fastq_file .trimmed.fastq).bam" # prepare results files
	./TMAP/tmap mapall \
	-f data/genome/hg38_masked.fa \
	-r $fastq_file \
	-v -Y -u -o 0 stage1 map4 |\
	samtools view -b |\
	samtools sort -o data/bam/$RESULTS_FILE
done
rm fastq_names.txt
# zip all fastq files again to save disk space:
ls data/fastq_trimmed/* | xargs -n1 bgzip


# PART 3: QUALITY FILTERING OF BAM FILES
# Check which decimal number corresponds to the bitwise flag "unmapped":
samtools flags unmap
# 0x4	4	UNMAP
# Create for-loop that will filter out (1) all unmapped reads and (2) all reads with a mapping quality below 10
mkdir data/bam_filtered
ls data/bam/* > bam_names.txt
BAM_FILES=($(cat bam_names.txt)) # load names into bash array
echo ${BAM_FILES[@]} # just for verification
for bam_file in ${BAM_FILES[@]}
do
	results_file="$(basename $bam_file .bam).filtered.bam"
	samtools view -b -F 4 $bam_file | \
	samtools view -bq 10 > data/bam_filtered/$results_file
done
rm bam_names.txt
