#!/usr/bin/bash

# Download and process hg38 genome assembly gap file (lists regions of the genome that are unassembled or unknown)
if [ ! -f ../data/raw/gap.txt ]; then
    echo "Gap file not found. Downloading..."
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gap.txt.gz -P ../data/raw/
    gunzip ../data/raw/gap.txt.gz
else
    echo "Gap file already exists. Skipping download."
fi







# Extract the relevant columns (chromosome, start, end) and convert to BED format
awk 'BEGIN {OFS="\t"} {print $2, $3, $4}' ../data/raw/gap.txt > ../data/processed/hg38_gaps.bed
bedtools sort -i ../data/processed/hg38_gaps.bed > ../data/processed/sorted_hg38_gaps.bed

echo "Sorting and merging the gap regions..."

# cat ../data/processed/sorted_hg38_gaps.bed ../data/raw/authors_bed_hg38/sorted_peaks.bed  > ../data/processed/combined_exclusion.bed
awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3}' ../data/processed/sorted_hg38_gaps.bed ../data/raw/authors_bed_hg38/sorted_peaks.bed > ../data/processed/combined_exclusion.bed

echo "Filtering combined exclusion regions to standard chromosomes..."
cut -f1 ../data/raw/hg38.chrom.sizes > ../data/processed/allowed_chroms.txt
awk 'NR==FNR { chroms[$1]=1; next } chroms[$1]' ../data/processed/allowed_chroms.txt ../data/processed/combined_exclusion.bed > ../data/processed/filtered_combined_exclusion.bed




echo "Sorting combined exclusion regions..."

bedtools sort -i ../data/processed/filtered_combined_exclusion.bed -g ../data/raw/hg38.chrom.sizes > ../data/processed/sorted_combined_exclusion.bed

echo "Merging combined exclusion regions... "
# Merge any overlapping intervals in the combined exclusion file to create a single set of regions to avoid:

bedtools merge -i ../data/processed/sorted_combined_exclusion.bed > ../data/processed/merged_excluded_regions.bed

# Compute safe regions
bedtools complement -i ../data/processed/merged_excluded_regions.bed -g ../data/raw/hg38.chrom.sizes > ../data/processed/genome_nonpeak_regions.bed

grep -E '^(chr[0-9]+\b|chr[XYM]\b)' ../data/processed/genome_nonpeak_regions.bed > ../data/processed/genome_nonpeak_regions_.bed 
rm ../data/processed/genome_nonpeak_regions.bed 
mv ../data/processed/genome_nonpeak_regions_.bed ../data/processed/genome_nonpeak_regions.bed


# bedtools shuffle -i ../data/processed/peaks_for_analysis.bed -g ../data/raw/hg38.chrom.sizes -incl ../data/processed/genome_nonpeak_regions.bed > ../data/processed/unbound_regions.bed
bedtools shuffle -i ../data/processed/peaks_for_analysis_all.bed -g ../data/raw/hg38.chrom.sizes -incl ../data/processed/genome_nonpeak_regions.bed > ../data/processed/unbound_seq.bed


# bedtools getfasta -fi ../data/raw/hg38.fa -bed ../data/processed/unbound_regions.bed -fo ../data/processed/unbound_regions.fa
bedtools getfasta -fi ../data/raw/hg38.fa -bed ../data/processed/unbound_seq.bed -fo ../data/processed/unbound_seq.fa

