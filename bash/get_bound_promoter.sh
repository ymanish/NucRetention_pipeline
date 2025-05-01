#!/usr/bin/bash

Annotation_dir="../data/annotations"
raw_dir="../data/raw"
processed_dir="../data/processed"

# Download the hg38 reference genome
if [ ! -f ${Annotation_dir}/hg38.fa ]; then
    echo "hg38 reference genome not found. Downloading..."
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P ${Annotation_dir}
    gunzip ${Annotation_dir}/hg38.fa.gz
else
    echo "hg38 reference genome already exists. Skipping download."
fi

# Download Gencode GTF file
if [ ! -f ${Annotation_dir}/gencode.v47.annotation.gtf ]; then
    echo "Downloading Gencode v47 annotation GTF..."
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz -P ${Annotation_dir}
    gunzip ${Annotation_dir}/gencode.v47.annotation.gtf.gz
else
    echo "Gencode v43 annotation GTF already exists. Skipping download."
fi

# Sort peaks by chromosome and start position
bedtools sort -i ${raw_dir}/authors_bed_hg38/GSM1346529_peaks_hg38_chrfiltered.bed -g ${Annotation_dir}/hg38.chrom.sizes > ${raw_dir}/authors_bed_hg38/sorted_peaks.bed

# Sort all peaks by score (column 5, descending)
sort -k5,5nr ${raw_dir}/authors_bed_hg38/sorted_peaks.bed > ${processed_dir}/peaks_for_analysis_all.bed


awk '$3 == "transcript" && $0 ~ /transcript_type "protein_coding"/' ${Annotation_dir}/gencode.v47.annotation.gtf > ${processed_dir}/protein_coding_transcripts.gtf


# Create promoter regions BED file 
echo "Creating promoter regions BED file..."
python ../src/utils/create_promoter_bed.py





# Intersect peaks with promoter regions
echo "Identifying peaks in promoter regions..."
bedtools intersect -a ${processed_dir}/promoter_regions.bed -b ${processed_dir}/peaks_for_analysis_all.bed -u > ${processed_dir}/promoter_peaks.bed

# Extract sequences for promoter peaks
echo "Extracting sequences for promoter peaks..."
bedtools getfasta -fi  ${Annotation_dir}/hg38.fa -bed ${processed_dir}/promoter_peaks.bed -fo ${processed_dir}/promoter_bound_seq.fa

echo "Done! Promoter peak sequences saved to ${processed_dir}/promoter_bound_seq.fa"

# Only report those entries in A that have _no overlaps_ with B.

bedtools intersect -a ${processed_dir}/promoter_regions.bed -b ${processed_dir}/peaks_for_analysis_all.bed -v > ${processed_dir}/unbound_promoters.bed
bedtools getfasta -fi  ${Annotation_dir}/hg38.fa -bed ${processed_dir}/unbound_promoters.bed -fo ${processed_dir}/unbound_promoters_seq.fa

echo "Done! Unbound Promoter peak sequences saved to ${processed_dir}/unbound_promoters_seq.fa"
