!/bin/sh -x
set -e
set -u
set -o pipefail

# ----------------------------------------------------------------- #
#				Variant Calling with bcftools call					#
# ----------------------------------------------------------------- #


# PART 1: FORMAT CONVERSION FROM BAM TO VCF
# I switch off the Base Alignment Quality (BAQ) algorithm (which is a standard part of samtools mpileup), as I already quality filtered my BAM-files.
mkdir data/vcf
ls data/bam_filtered/* > bam_names.txt
BAM_FILES=($(cat bam_names.txt)) # load names into bash array
echo ${BAM_FILES[@]} # just for verification
for bam_file in ${BAM_FILES[@]}
do
	results_file="$(basename $bam_file .filtered.bam).mpileup.vcf.gz"
	samtools mpileup -v -t AD,DP --no-BAQ --fasta-ref data/genome/hg38_masked.fa $bam_file \
	> data/vcf/$results_file
done
rm bam_names.txt
# -v = output is a VCF-file
# -t AD = output includes Allelic Depth (number of reads that support the reported allele)
# -t DP = output includes filtered Depth (=number of filtered reads that support the reported allele)


# PART 2: CHANGE SAMPLE NAME IN TABLE HEADER FROM "NOSM" TO MORE MEANINGFUL SAMPLE NAMES
# create directory
mkdir data/vcf_1_sample_names
# retrieve sample names
ls data/vcf/* > vcf_names.txt
VCF_FILES=($(cat vcf_names.txt)) # load names into bash array
echo ${VCF_FILES[@]} # just for verification
for vcf_file in ${VCF_FILES[@]}
do 
	results_file="$(basename $vcf_file .mpileup.vcf.gz).renamed.mpileup.vcf.gz"
	basename $vcf_file | sed 's/.mpileup.vcf.gz//' > sample_name.txt
	bcftools reheader -s sample_name.txt $vcf_file -o data/vcf_1_sample_names/$results_file
done
rm vcf_names.txt sample_name.txt


# PART 3: VARIANT CALLING WITH bcftools call
mkdir data/vcf_called
ls data/vcf_1_sample_names/* > vcf_names.txt
VCF_FILES=($(cat vcf_names.txt)) # load names into bash array
echo ${VCF_FILES[@]} # just for verification
for vcf_file in ${VCF_FILES[@]}
do
	results_file="$(basename $vcf_file .renamed.mpileup.vcf.gz).called.vcf"
	bcftools call -m -g 10 -f GQ $vcf_file > data/vcf_called/$results_file
done
rm vcf_names.txt
# -m = multiallelic calling model
# -f GQ = include a quality score in case of variant genotypes


# PART 4: FILTER CALLED VARIANTS
mkdir data/vcf_filtered
ls data/vcf_called/* > vcf_names.txt
VCF_FILES=($(cat vcf_names.txt)) # load names into bash array
echo ${VCF_FILES[@]} # just for verification
for vcf_file in ${VCF_FILES[@]}
do
	results_file="$(basename $vcf_file .called.vcf).filtered.vcf"
	bcftools view -i 'MIN(FMT/DP)>10' $vcf_file > data/vcf_filtered/$results_file
done
rm vcf_names.txt
# -i 'MIN(FMT/DP)>10' = filter out all reads for which filtered depth < 11
