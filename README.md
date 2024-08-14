# Custom Code for Single-Cell and Single-Nucleus Analyses in Yim et al., 2024
Scripts used for analyzing embryonic single-cell RNA-sequencing and juvenile single-nucleus RNA-sequencing datasets generated from wild type and *Chd8* heterozygous cortical samples, alongside scripts for analyzing published human fetal cortex single-cell RNA-sequencing data, as reported in:

Kristina M Yim*, Marybeth Baumgartner*, Martina Krenzer*, María F. Rosales Larios, Guillermina Hill-Terán, Timothy Nottoli, Rebecca A Muhle, and James P Noonan (2024). *bioRxiv* date. [doi](link).
# Read Alignment and Filtering
For each single-cell or single-nucleus sample, fastq reads were aligned to the mm10 reference genome using [STAR](https://doi.org/10.1093/bioinformatics/bts635/) v2.7.9a. Reads mapping to more than 10 loci were then discarded, as were multi-mapping reads with more than one best scoring alignment. Aligned reads passing these criteria were used to generate single-cell and single-nucleus count matrices using STAR. 

The `read_alignment` directory contains all scripts for single-cell and single-nucleus read alignment, read filtering, and generating counts. Input fastq files are available through the Gene Expression Omnibus (GEO) under accession numbers GSE273271, for the single-cell dataset, and GSE273765, for the single-nucleus dataset.
## Single-Cell Data
- `align_reads_scRNA-seq.sh`
- `filter_multimap_reads.sh`
- `generate_counts_scRNA-seq.sh`
## Single-Nucleus Data
- `align_reads_snRNA-seq.sh`
- `filter_multimap_reads.sh`
- `generate_counts_snRNA-seq.sh`
# Single-Cell Analysis of Embryonic Mouse Cortex
The `R/mmCortex_scRNA-seq` directory contains all R scripts for preliminary data filtering, downsampling, normalization, integration, cluster analysis, pseudotime analysis, and metagene analysis for the scRNA-seq dataset. [Seurat](https://doi.org/10.1016/j.cell.2021.04.048) v4.0.4 was used to perform preliminary analysis, normalization, integration, and cluster analysis. Input counts matrices are available through GEO under accession number GSE273271. 
## Preliminary Analysis
For each sample in the scRNA-seq dataset, low quality cells and potential doublets were filtered by thresholding on the total feature counts per cell, the number of genes detected per cell, and the percentage of counts originating from mitochondrial RNA. Filtered count matrices were downsampled based on median feature counts by time point, for use in downstream metagene and differential expression analysis scripts, using [DropletUtils](https://doi.org/10.1038/s41467-018-05083-x) v1.14.2.
- `00_create_seurat_objects.R`
- `01_prelim_filter.R`
- `02_downsample_counts.R`
## Normalize and Integrate Data
Using filtered count matrices, percentage mitochondrial content and the difference between G2/M and S phase cell cycle scores were regressed out, and feature counts were normalized across each batch. Normalized data across all samples were integrated using a reference dataset, corresponding to the largest batch with equal sex representation. These integrated data were used for UMAP embedding and preliminary cluster identification.
- `03_SCTransform.R`
- `04_integrate_SCT_cca_RefF.R`
## Cluster Analysis
Based on preliminary clusters identified in the integrated dataset, cluster marker genes were computed. Marker genes with known cell type-specific expression were used to assign a cell type to each cluster (cluster labels). For each cell type, we compared the number of cells between wild type and *Chd8* heterozygous samples at each time point.
- `05a_cluster_markers_res0.4.R`
- `05b_cluster_labels.R`
- `05c_cluster_ttest.R`
## Pseudotime Analysis
The integrated dataset was used for PHATE embedding, using [phateR](https://doi.org/10.1038/s41587-019-0336-3) v1.0.7; trajectory identification; and computation of pseudotime along the trajectory of interest. Pseudotime was divided into bins to identify gene expression trajectories. For each time point and genotype, downsampled count matrices were used to calculate the mean expression of each gene in each pseudotime bin. Genes rarely detected in any of these pseudotime bins (i.e., expressed in few cells per bin) were excluded from downstream analysis.
- `06_compute_phate_pseudotime.R`
- `07_bin_expression.R`
## Metagene Analysis
Gene expression trajectories were centered, scaled, and transformed with symbolic aggregate approximation ([SAX](https://doi.org/10.1007/s10618-007-0064-z)) using [jMotif](https://github.com/jMotif/jmotif-R) v1.1.1. SAX-transformed gene expression trajectories in the wild type dataset were used for k-means clustering, and similar expression trajectories were aggregated into metagene centers by hierarchical clustering. Each SAX-transformed gene expression trajectory from the wild type and *Chd8* heterozygous datasets was then assigned to a metagene based on maximum Pearson correlation with the metagene centers.
- `08_compute_sax.R`
- `09_compute_metagenes.R`


# Single-Cell Analysis of Human Cortical Data
The `R/hgCortex` directory contains all R scripts for data normalization, pseudotime analysis, and metagene analysis of human fetal cortex scRNA-seq data published in [Polioudakis et al., 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6831089/). The raw count matrix is available through the Geschwind laboratory’s [Cortical Development Expression (CoDEX) viewer](http://solo.bmap.ucla.edu/shiny/webapp/) webtool. 

We log-normalized the count matrix, identified variable features, and regressed out the percentage mitochondrial content, the difference between G2/M and S phase cell cycle scores, library, donor, and total counts per cell before proceeding to UMAP and PHATE embedding. We utilized the published clusters and cell type assignments in our analysis. We selected cell types corresponding to those in the mouse scRNA-seq dataset for PHATE embedding, trajectory identification, pseudotime analysis, and metagene analysis. These scripts mirror those used for pseudotime and metagene analyses of the mouse scRNA-seq dataset. We calculated Pearson correlation between human and mouse metagene centers, then labeled human metagenes based on the maximal correlation with mouse metagene centers.  
- `01_Polioudakis_compute_pseudotime.R`
- `02_Polioudakis_bin_expression.R`
- `03_Polioudakis_sax_metagenes.R`
# Single-Nucleus Analysis of Juvenile Mouse Cortex
The `R/mmCortex_snRNA-seq` directory contains all R scripts for preliminary data filtering, normalization, integration, preliminary clustering, doublet removal, and cluster analysis for the snRNA-seq dataset. Input counts matrices are available through GEO under accession number GSE273765. 
## Preliminary Analysis
For each sample in the snRNA-seq dataset, low quality cells were filtered by thresholding on the percentage of counts originating from mitochondrial RNA. 
- `00_create_seurat_objects.R`
- `01_prelim_filter.R`
## Normalize and Integrate Data
Using filtered count matrices, percentage mitochondrial content was regressed out, and feature counts were normalized across each batch. Normalized data across all samples were integrated, and these integrated data were used for UMAP embedding and preliminary cluster identification.
- `02_SCTransform.R`
- `03_integrate_SCT_cca.R`
## Preliminary Cluster Analysis
Based on clusters identified in the integrated dataset, cluster marker genes were computed. Genes with known cell type-specific expression were used to assign a preliminary cell type to each cluster (cluster labels). 
- `04a_cluster_markers_res0.5.R`
- `04b_cluster_labels.R`
## Doublet Identification
Preliminary clusters identified in the integrated dataset were used for by-sample, cluster-informed doublet identification using [scDblFinder](https://doi.org/10.12688/f1000research.73600.2) v.1.18.0. Doublets were removed from the integrated dataset, and the resulting singlet-only integrated dataset was used for downstream analyses. 
- `05_doublet_removal.R`
## Cluster Filtering and Analysis
Using the singlet-only integrated dataset, cluster marker genes were re-computed. Marker genes with known cell type-specific expression were used to verify the preliminary cell type assigned to each cluster. For cluster filtering, we assessed the ambiguity of marker genes identified in each cluster, the total number of nuclei in each cluster, and the proportion of nuclei in each cluster that originated from each batch. For each cluster that passed these filters, we compared the number of nuclei between wild type and *Chd8* heterozygous samples.
- `06a_cluster_markers_singlets_res0.5.R`
- `06b_cluster_ttest.R`
# Differential Expression Analysis of Mouse Cortical Datasets
The `R/monocle` directory contains all R scripts for data pre-processing and differential expression analysis of the embryonic mouse scRNA-seq and juvenile mouse snRNA-seq datasets, using [Monocle 3](https://doi.org/10.1038/s41586-019-0969-x).
## Data Pre-Processing
For the embryonic scRNA-seq dataset, downsampled counts were used to generate a Monocle 3-compatible cell_data_set (cds) object for each time point, each of which was split further by trajectory or cell type. Cell types were annotated based on pseudotime bins. For the juvenile snRNA-seq dataset, the singlet-only data were used to generate a cds object, which was split by clusters/cell types that passed all upstream filters. For each cds object, genes rarely detected in the wild type and *Chd8* heterozygous cells/nuclei were excluded from differential expression testing.
- `create_cds.R`
- `subset_cds_comparisons.R`
- `create_cds_P25.R`
## Differential Expression Testing
For each time point, Monocle 3 was run on each filtered cds object from upstream scripts using an additive model. Output average log2(fold change) for each gene was calculated by Seurat’s FoldChange function in the pre-processing scripts.
- `monocle_E12.R`
- `monocle_E14.R`
- `monocle_E16.R`
- `monocle_E17.R`
- `monocle_P25.R`
