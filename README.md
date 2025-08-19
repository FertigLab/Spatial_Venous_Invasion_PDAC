# Spatial profiling of human pancreatic ductal adenocarcinoma reveals molecular alterations associated with venous invasion 

## Alexander T.F. Bell, Peter Chianchiano, Katsuya Hirose, Daniel Salas-Escabillas,  Doreen Gisder, Jacob T. Mitchell, Ludmila Danilova, Alexander I. Damanakis, Joseph Tandurella, Jeanette Johnson, Qingfeng Zhu, Robert A. Anders, Jing Zhu, James R. Eshleman, Luciane T. Kagohara, Elana J. Fertig, Laura D. Wood

### Abstract
In pancreatic ductal adenocarcinoma (PDAC), venous invasion (VI) is a critical step in metastasis and is associated with poor survival. However, little is known about the molecular features of VI. To investigate, we performed spatial transcriptomic analysis of 95 human PDAC tissue samples from eight treatment-naïve patients. Our analysis revealed that compared to PDAC in stroma, VI demonstrated upregulation of genes associated with epithelial differentiation, classical subtype, and benign exocrine function. Conversely, VI demonstrated downregulation of genes associated with mesenchymal differentiation, basal-like subtype, and disease aggression. Additionally, we uncovered characteristics of VI morphology that correlated with these molecular features. VI-intraepithelial-neoplasia-like (VI-IN-Like) had preserved venous architecture and possessed a classical, epithelial molecular phenotype. VI-Destructive had destroyed venous architecture and possessed a more basal-like, mesenchymal phenotype than VI-IN-Like. We contextualized our findings using public RNA-seq data and observed that metastatic PDAC had greater similarity to PDAC in stroma than to VI, whereas circulating tumor cells showed no preferential association. We confirmed our findings by spatial proteomic analysis of VI in an independent cohort of 19 treatment-naïve patients with PDAC. Overall, our work provides a reference atlas of spatial transcriptomics and proteomics of VI and reveals paradoxical unexpected increases in molecular features associated with better patient outcomes.

### Data Acquisition

- GeoMx Digital Spatial Profiler Patient Cohort
  - The following files are on Zenodo:
    - GeoMx_raw_counts.csv which contains the raw unprocessed count matrix 
    - ProbeQC_merged_batches.csv which contains the count matrix which underwent filtration of low quality probes in the Nanostring cloud based data             analysis portal (this is the input data matrix to 01_Preprocessing)
    - fully_batch_corrected_vsd.csv which contains the fully processed matrix with low quality probes removed, low count genes removed, one mislabeled           sample removed, normalized through variance stabilizing transformation and inter-patient and inter-batch differences regressed out using ComBat.
    - initial_metadata.csv which contains incomplete cohort metadata and is used as an input file for 01_Preprocessing. 
    - metadata_with_VI_subtypes_and_NGS.csv which contains the complete cohort metadata (Tissue type, VI subtypes, SMAD4/p53/KRAS/p16 status).
    - batch_merged_cogapsP4.rds contains CoGAPS outputs and is used as an input for several scripts.
    - annotations.csv is an input used for running script 03_VI_subtypes.R. This file contains VI subtype annotations coded as A (VI-Destructive), B             (VI-IN-Like), C (VI-Conventional), X (Not VI).
    - NGS_table.csv contains base and amino acid changes for all the called DNA mutations in the clinical cancer panel

- PDAC atlas data as a monocle3 object
  - The single-cell RNA-seq data set is comprized of six aggregated data sets. Aggregation and annotation of the data is described in [Guinn et al](https://pubmed.ncbi.nlm.nih.gov/38587552/). This analysis utilized cells contributed 
    from [Steele et al](https://pubmed.ncbi.nlm.nih.gov/34296197/) and [Peng et al](https://pubmed.ncbi.nlm.nih.gov/31273297/) to the atlas. Scripts used to generate these data are avaialble on [GitHub](https://github.com/FertigLab/PDAC_Atlas).

- Bulk RNA-Seq of metastatic and primary PDAC
  - Bulk-RNA-Seq samples from primary and metastatic PDAC from the cohort described in [Connor et al.](https://pubmed.ncbi.nlm.nih.gov/30686769/), were retrieved from the European Genome Phenome Archive (EGAS00001002543) and from International Consortium Cancer Genome data portal (PACA-CA). Additional metadata was retrieved from Connor et al. supplementary materials. 

- Circulating pancreatic tumor cells
  - Bulk RNA-Seq samples from circulating pancreatic tumor cells from the cohort described in [Franses et al.](https://pubmed.ncbi.nlm.nih.gov/32620742/) were retrieved from NCBI-GEO under accession number GSE144561.
  - Both GSE144561_DESeq2_NormalizedCountsForAllSamples.txt.gz and GSE144561_rawCountsAllsamples.txt.gz were used as inputs for script 06_CTC_Projection.R

- Protein-RNA correlation 
  - Inputs for 09_CPTAC_RNA_Correlations.R were retrieved from Proteomic Data Commons ID: PDC000270 and downloaded from [LinkedOmics](https://www.linkedomics.org/data_download/CPTAC-PDAC/). Files needed to run code: mRNA_RSEM_UQ_log2_Tumor.tsv, proteomics_gene_level_MD_abundance_tumor.tsv, molecular_pheno_tumor.csv. 

- PhenoCycler multiplex protein analysis patient cohort
  - All files are available on Zenodo:
    - .qptiff output files from PhenoCycler that are used as inputs to HALO software
    - geojson_annotations.zip contains manual annotations for individual ROIs, to accompany each .qptiff file
    - The HALO output used for downstream analysis (as the input to script 10_PhenoCycler.R) is available as file raw_HALO_outputs.csv
      - NOTE: The .qptiff files contain markers (CD31, EpCAM, Annexin8, Urokinase) and tissue types (PNI, ND, duodenum, etc.) which are not included in           the downstream analysis, and are not present in file raw_HALO_outputs.csv.
    - rshiny_phenocycler_annotations.zip which contains the manually annotated VI foci. These files are used as inputs to script 10_PhenoCycler.R
    - processed_phenocyler_data.csv which contains PanCK positive cells with VI/PDAC foci annotated, cell-level Moffitt classifications, foci-level             moffitt classifications. 

## Scripts
### All scripts are located at SAKJDHAKSDHKAJSHDKAJSHDLAKSJHDAKSJDHALSKJDHALSKJHDALSKJDHALSKJDHASLKDJHALSKJDHALSKJHD

Execute scripts in the order that they are numbered, as outputs will be used as inputs for later scripts. The three scripts with the prefix "00" contain custom functions that will be imported at the beginning of each script. These scripts do not produce any objects and do not have to be run independently. 

For compatible software and package versions, please see 00_package_versions.html

### 00_Custom_Functions.R

Contains custom functions necessary for all analysis of GeoMx data

### 00_Custom_Functions_CTC.R

Contains custom functions necessary only for analysis of the projection of GeoMx CoGAPS patterns onto circulating tumor cell data (analysis performed in script 06_CTC_Projection.R)

### 00_Custom_Functions_atlas.R

Contains custom functions necessary for analysis involving the projection of GeoMx CoGAPS patterns onto scRNA-Seq PDAC atlas (analysis performed in script 07_Atlas_projection.R)

### 01_Preprocessing.R

Imports the GeoMx data as raw counts and performs normalization and batch correction and produces files necessary for many downstream analyses. This script produces Figure S1A and Figure 2B.

### 02_VI_PDAC.R

Performs the majority of the GeoMx analysis, including differential expression between combined PDAC and VI, VI and PDAC, CoGAPS analysis, and pathway analysis. This script produces Figures 2C-2D, 3A-3F, 4A-4C, 5A-5E, S5A-S5C, S3A, S2A.

### 03_VI_subtypes.R

Performs analysis of VI subtypes (VI-IN-Like, VI-Destructive). This script produces Figure 6B-6E, S6A-S6E, S7A-S7B, and S8A-S8C.

### 04_left_over_icgc_samples.R

Produces a file needed to organize the patient cohort analyzed in 05_metastatic_projection.R. 

### 05_metastatic_projection.R

Analyzes a projection of the GeoMx CoGAPS patterns onto the bulk RNA-Seq PDAC cohort described in Connor et al. 2019. Produces Figure 7A-7C, S10A-S10D.

### 06_CTC_Projection.R

Analyzes a projection of the GeoMx CoGAPS patterns onto the bulk RNA-Seq circulating tumor cell cohort described in Franses et al. 2020. Produces Figure S11A-S11C, 7D-7E. 

### 07_Atlas_projection.R

Analyzes a projection of the GeoMx CoGAPS patterns onto the scRNA-Seq PDAC atlas composed of data generated by Peng et al. 2019, Steele et al. 2020, and assembled by Guinn et al. 2024. Produces Figure S4A-S4C.

### 08_Clinical_NGS.R

Analyzes clinical next generation sequencing for each patient in the GeoMx cohort. Analyzes staining for p53 and SMAD4. Produces Figure S9A-S9E.

### 09_CPTAC_RNA_Protein_Correlations.R

Derives RNA-protein correlation in PDAC for the VI marker genes to help select candidate markers for protein validation. 

### 10_PhenoCycler.R

Analyzes multiplex protein staining (PhenoCycler) in independent PDAC cohort. Produces Figure 8B-8C and S12A-S12D. 