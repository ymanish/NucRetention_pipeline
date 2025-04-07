#!/usr/bin/bash

# Download SRA data for SRP025289 project
mkdir -p ../data/raw

# Define the SRR accessions for the project
SRR_ACCESSIONS=("SRR896001" "SRR1187944") # SRR896001 is human sperm and SRR1187944 four anonymised fertile donors

# Loop through each accession
for SRR in "${SRR_ACCESSIONS[@]}"; do
  # Download the SRA file with prefetch
#   prefetch "$SRR" -O ../data/raw

  # Extract FASTQ files using fasterq-dump
  fasterq-dump ../data/raw/"$SRR" -O ../data/raw/"$SRR" 

  # Optional: Remove the .sra file to save space
  rm ../data/raw/"$SRR"/"$SRR".sra
done