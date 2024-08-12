#!/bin/bash

# STAR v2.7.9a
# Generate count matrices from filtered BAM files
STAR --runThreadN 6 \
--genomeDir /path/to/GENOME/mm10_cellranger_2020-A/star/ \
--readFilesType SAM SE \
--readFilesIn Aligned_unique.bam,Aligned_multi_primary_n1.bam \
--readFilesPrefix /path/to/mmCortex_scRNA-seq/alignments_AllBestScore/SAMPLE/  \
--readFilesCommand samtools view -F 0x100 \
--outFileNamePrefix /filter_multimap_counts/SAMPLE/ \
--outSAMtype None \
--soloType CB_UMI_Simple \
--soloCBwhitelist /path/to/737K-august-2016.txt \
--soloBarcodeReadLength 1 \
--soloCBmatchWLtype 1MM_multi \
--soloInputSAMattrBarcodeSeq CR UR \
--soloInputSAMattrBarcodeQual CY UY \
--soloMultiMappers Unique \
--soloUMIdedup 1MM_CR \
--soloUMIfiltering - \
--soloCellFilter CellRanger2.2 3000 0.99 10
