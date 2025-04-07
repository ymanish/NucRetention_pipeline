#!/usr/bin/bash

# Check if the chain file already exists
if [ ! -f ../data/raw/hg19ToHg38.over.chain.gz ]; then
    echo "Chain file not found. Downloading..."
    # Download chain file for liftover
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -P ../data/raw/
else
    echo "Chain file already exists. Skipping download."
fi


# Liftover peaks from hg19 → hg38
mkdir -p ../data/raw/authors_bed_hg38

liftOver ../data/raw/authors_bed/GSM1346529_Human_146bp_peaks.bed ../data/raw/hg19ToHg38.over.chain.gz ../data/raw/authors_bed_hg38/GSM1346529_peaks_hg38.bed ../data/raw/authors_bed_hg38/GSM1346529_unmapped.bed

grep -E '^(chr[0-9]+\b|chr[XYM]\b)' ../data/raw/authors_bed_hg38/GSM1346529_peaks_hg38.bed > ../data/raw/authors_bed_hg38/GSM1346529_peaks_hg38_chrfiltered.bed
	