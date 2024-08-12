#!/bin/bash

# STAR v2.7.9a
# Cell Ranger Reference, 2020-A
# Mouse reference, mm10 (GENCODE vM23/Ensembl 98)
# Input: single-cell RNA-seq data generated using 10x Genomics Chromium Next GEM Single-Cell 3’ Reagent Kits for v2 chemistry

# Generate genome indices
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir /path/to/GENOME/mm10_cellranger_2020-A/star/ \
--genomeFastaFiles /path/to/GENOME/mm10_cellranger_2020-A/fasta/genome.fa \
--sjdbGTFfile /path/to/GENOME/mm10_cellranger_2020-A/genes/genes.gtf \
--sjdbOverhang 99 \
--genomeSAsparseD 3

# Align reads
STAR --runThreadN 6 \
--genomeDir /path/to/GENOME/mm10_cellranger_2020-A/star/ \
--readFilesManifest /path/to/mmCortex_scRNA-seq/file_manifests/SAMPLE_readFiles.tsv \
--readFilesCommand zcat \
--outFileNamePrefix /path/to/mmCortex_scRNA-seq/alignments_AllBestScore/SAMPLE/ \
--outSAMtype BAM Unsorted \
--outSAMattributes NH HI AS nM CR CY UR UY GX GN \
--outSAMprimaryFlag AllBestScore \
--outSAMmultNmax 10 \
--outBAMcompression 10 \
--soloType CB_UMI_Simple \
--soloCBwhitelist /path/to/737K-august-2016.txt \
--soloBarcodeReadLength 1 \
--soloCBmatchWLtype 1MM_multi \
--soloMultiMappers Unique \
--soloUMIdedup 1MM_CR \
--soloUMIfiltering - \
--soloCellFilter CellRanger2.2 3000 0.99 10
