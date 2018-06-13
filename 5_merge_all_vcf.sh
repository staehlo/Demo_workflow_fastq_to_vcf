!/bin/sh -x
set -e
set -u
set -o pipefail

# ----------------------------------------------------------------- #
#							Merge all VCFs							#
# ----------------------------------------------------------------- #

# PART 1: MERGE ALL VCF-files
# Compress all called VCF-files with BGZIP:
ls data/vcf_filtered/*.vcf | xargs -n1 bgzip
# Index all compressed VCF-files with bcftools index:
ls data/vcf_filtered/*.vcf.gz | xargs -n1 bcftools index
# join all VCF files into on big file
bcftools merge data/vcf_filtered/*.vcf.gz > data/merged.vcf
# Count number of positions which could be determined after filtering:
grep -v "^#" data/merged.vcf | wc -l
# -> 3038

# PART 2: CREATE OVERVIEW TABLE WITH MOST COMMON SNPs:
# extract the following data from the combined table:
# Chrom - Pos - Number of undefined genotypes - Number of homozygote wildtype:
grep -v "^#" data/merged.vcf | awk -F'\t' 'BEGIN{print "#CHROM", "POS", "UNDEFINED", "homozygote_wt"}{print $1, $2, gsub("\\.\\/\\.:\\.","\t" NR), gsub("0\\/0","\t" NR)}' > aaa.txt
# Add another columns with the number of defined alternative genotypes:
grep -v "^#" aaa.txt | awk 'BEGIN{print "#CHROM", "POS", "UNDEFINED", "homozygote_wt", "homo/hetero_mt"}{print $1, $2, $3, $4, 89-$3-$4}' > bbb.txt
# reformat table:
tr ' ' \\t < bbb.txt > ccc.txt
# sort table:
grep "^#" ccc.txt > data/variant_frequency.txt
grep -v "^#" ccc.txt | sort -k5,5rV >> variant_frequency.txt
rm aaa.txt bbb.txt ccc.txt
