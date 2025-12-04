#!/bin/bash
set -e

bam_dir="../SDR"
window_bed="windows.bed"
matrix_dir="coverage_matrix"

min_zero_samples=55
merged_output="merged_candidate_regions.bed"
genome_fa="male_hap1.fa" 
fasta_output="merged_candidate_regions.fa"

mkdir -p "$matrix_dir"

for bam in "$bam_dir"/*.bam; do
    sample=$(basename "$bam" .bam)
    out_cov="$matrix_dir/${sample}.cov"
    echo " -> $sample"
    bedtools coverage -a "$window_bed" -b "$bam" -mean \
      | awk -v OFS='\t' '{zero=($4==0)?1:0; print $1, $2, $3, zero}' > "$out_cov"
done

paste "$matrix_dir"/*.cov | \
awk -v OFS='\t' '{
    chr=$1; start=$2; end=$3;
    zero_count=0;
    for(i=4; i<=NF; i+=4) {
        zero_count += $(i)
    }
    print chr, start, end, zero_count
}' > window_zero_counts.txt

awk -v th="$min_zero_samples" '$4 >= th { print $1, $2, $3 }' window_zero_counts.txt > candidate_windows.bed

bedtools merge -i candidate_windows.bed > "$merged_output"

bedtools getfasta -fi "$genome_fa" -bed "$merged_output" -fo "$fasta_output"

