# Split BAM into uniquely mapped and multi-mapped reads
samtools view -@ 5 -q 255 -o Aligned_unique.bam -U Aligned_multi.bam Aligned.out.bam

# For multi-mapped reads, count the number of primary alignments (-F 0x100)
samtools view -@ 5 -F 0x100 Aligned_multi.bam | cut -f1,12,14 | uniq -c > count_primary_alignments.txt

# Get reads names for single top hit alignments
awk '$1 == 1 { print $2 }' count_primary_alignments.txt > keep_reads_n1.txt

# Get reads that match filter criteria
samtools view -@ 5 -F 0x100 -N keep_reads_n1.txt -o Aligned_multi_primary_n1.bam Aligned_multi.bam