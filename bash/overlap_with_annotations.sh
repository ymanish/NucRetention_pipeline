# Split input BED into clusters
awk 'BEGIN {OFS="\t"} $3 == 0 {print $6,$7,$8}' ../notebooks/energies_with_clusters.txt > ../notebooks/cluster0.bed
awk 'BEGIN {OFS="\t"} $3 == 1 {print $6,$7,$8}' ../notebooks/energies_with_clusters.txt > ../notebooks/cluster1.bed

# Calculate overlaps
annotations=("cpg_islands" "repeats" "genes")

for ann in "${annotations[@]}"; do
    for cluster in 0 1; do
        bedtools intersect -a ../notebooks/cluster${cluster}.bed -b ../data/annotations/${ann}.bed -u > ../notebooks/cluster${cluster}_${ann}.bed
        count=$(wc -l < ../notebooks/cluster${cluster}_${ann}.bed)
        total=$(wc -l < ../notebooks/cluster${cluster}.bed)
        echo "cluster${cluster},${ann},${count},${total}" >> ../notebooks/annotation_counts.csv
    done
done

# # Create BED files for visualization
# bedtools merge -i cluster0.bed > cluster0_merged.bed
# bedtools merge -i cluster1.bed > cluster1_merged.bed
