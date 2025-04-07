#!/usr/bin/bash
# Download the hg38 reference genome
if [ ! -f ../data/raw/hg38.fa ]; then
    echo "hg38 reference genome not found. Downloading..."
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P ../data/raw/
    gunzip ../data/raw/hg38.fa.gz
else
    echo "hg38 reference genome already exists. Skipping download."
fi



# Sort by chromosome (col 1) and start position (col 2, numeric)
# sort -k1,1 -k2,2n ../data/raw/authors_bed_hg38/GSM1346529_peaks_hg38_chrfiltered.bed > ../data/raw/authors_bed_hg38/sorted_peaks.bed

bedtools sort -i ../data/raw/authors_bed_hg38/GSM1346529_peaks_hg38_chrfiltered.bed -g ../data/raw/hg38.chrom.sizes > ../data/raw/authors_bed_hg38/sorted_peaks.bed

# Use all peaks
# cp sorted_peaks.bed peaks_for_analysis.bed

# Select top 10,000 peaks by score 
# sort -k5,5nr ../data/raw/authors_bed_hg38/sorted_peaks.bed | head -n 10000 > ../data/processed/peaks_for_analysis.bed
sort -k5,5nr ../data/raw/authors_bed_hg38/sorted_peaks.bed > ../data/processed/peaks_for_analysis_all.bed


# bedtools getfasta -fi ../data/raw/hg38.fa -bed ../data/processed/peaks_for_analysis.bed -fo ../data/processed/bound_top10000.fa
bedtools getfasta -fi ../data/raw/hg38.fa -bed ../data/processed/peaks_for_analysis_all.bed -fo ../data/processed/bound_seq.fa
