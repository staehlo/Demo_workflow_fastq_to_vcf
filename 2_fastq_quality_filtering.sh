!/bin/sh -x
set -e
set -u
set -o pipefail

# ----------------------------------------------------------------- #
#		Quality filtering of FASTQ-files with Trimmomatic			#
# ----------------------------------------------------------------  #

 # Create folder for trimmed FASTQ-files:
mkdir data/fastq_trimmed

# Trim all files with trimmomatic:
ls data/fastq/* > fastq_names.txt
FASTQ_FILES=($(cat fastq_names.txt)) # load names into bash array
echo ${FASTQ_FILES[@]} # just for verification
for fastq_file in ${FASTQ_FILES[@]}
do
	results_file="$(basename $fastq_file .fastq.gz).trimmed.fastq.gz"
	java -jar trimmomatic-0.38.jar \
		SE \
		-phred33 \
		$fastq_file \
		data/fastq_trimmed/$results_file \
		LEADING:5 TRAILING:5 \
		SLIDINGWINDOW:4:15 \
		MINLEN:40
done
# The following options were chosen:
# SE = single end
# phred33 = phred33 quality coding (which is used by the Ion Torrent PGM software)
# LEADING:5 TRAILING:5 = cut leading and trailing with a PHRED-score < 5
# SLIDINGWINDOW:windowSize=4:requiredQuality=15
# MINLEN:40 = remove all sequences whose length < 40 after the previous operations
rm fastq_names.txt

