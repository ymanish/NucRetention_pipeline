#!/usr/bin/bash

mkdir -p ../data/raw/authors_bed

# Download the BED file archive (already called peaks)
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE47nnn/GSE47843/suppl/GSE47843_RAW.tar

tar -xvf GSE47843_RAW.tar -C ../data/raw/authors_bed

rm GSE47843_RAW.tar