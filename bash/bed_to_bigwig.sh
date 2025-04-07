#!/usr/bin/bash

# BED → BedGraph
awk '{print $1 "\t" $2 "\t" $3 "\t" $5}' ../data/raw/authors_bed_hg38/GSM1346529_peaks_hg38_chrfiltered.bed > scores.bedgraph

# Sort the BEDGraph file
sort -k1,1 -k2,2n  scores.bedgraph > scores.sorted.bedgraph

bedtools merge -i scores.sorted.bedgraph -c 4 -o mean > scores.sorted.merged.bedgraph


#Download chromosome sizes for hg38
if [ ! -f ../data/raw/hg38.chrom.sizes ]; then
    echo "Chromosome sizes file not found. Downloading..."
    # Download chromosome sizes for hg38
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes -P ../data/raw/
else
    echo "Chromosome sizes file already exists. Skipping download."
fi

# wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
# chmod a+x ./bedToBigBed

# BedGraph → BigWig
bedGraphToBigWig scores.sorted.merged.bedgraph ../data/raw/hg38.chrom.sizes ../data/raw/authors_bed_hg38/GSM1346529_peaks_hg38_chrfiltered.bw

rm scores.bedgraph scores.sorted.bedgraph scores.sorted.merged.bedgraph



