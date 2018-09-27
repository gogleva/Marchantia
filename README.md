# Marchantia

Scripts for Marchantia II paper, mainly transcriptomics

**Scripts:**

**``DEG_analysis_Mpolymorpha_Ppalmivora.R``**

- Analysis of differentially expressed genes in M.polymorpha-P.palmivora (ARI-tdTomato) timecourse (1,2,3 and 4 dpi with correspondent mock samples).
- Volcano plots for pair-wise comparisons;
- Upset plots for up- and down-regulated genes;
- Heatmap, showing overall expression patterns of DEGs, identified in all pair-wise comparisons.

  This script requires 2 input data files:
  - ``data/FeartureCounts_STAR_Mpolymorpha.tsv`` - raw counts
  - ``data/sample_table_Mpolymorpha.csv`` - description of samples and conditions in M.polymorpha-P.palmivora (ARI-tdTomato) timecourse


**``DEG_analysis_Nbenthamiana_Ppalmivora.R``**
- Analysis of differentially expressed genes in N.benthamiana-P.palmivora (ARI-tdTomato) timecourse (14,24,48 and 72 hpi with correspondent mock samples).
- Volcano plots for pair-wise comparisons;
- Upset plots for up- and down-regulated genes;
- Heatmap, showing overall expression patterns of DEGs, identified in all pair-wise comparisons.

  This script requires 2 input data files:
  - ``data/FeartureCounts_STAR_Nbenthamiana.csv`` - raw counts
  - ``data/sample_table_Nbenthamiana.csv`` - description of samples and conditions in M.benthamiana-P.palmivora (ARI-tdTomato) timecourse

**``ortholog_expression_analysis_Mpolymorpha_Nbenthamiana.R``**
- DEG analysis, filter genes by adjusted p-value, keep all LFC
- parse OrthoDinfer outputs
- visualise expression patterns of single-copy orthologs during comparable stages of infection

  This script requires the following input files:
  - ``data/FeartureCounts_STAR_Mpolymorpha.tsv`` - raw counts
  - ``data/sample_table_Mpolymorpha.csv`` - description of samples and conditions in M.polymorpha-P.palmivora (ARI-tdTomato) timecourse
  - ``data/FeartureCounts_STAR_Nbenthamiana.csv`` - raw counts
  - ``data/sample_table_Nbenthamiana.csv`` - description of samples and conditions in M.benthamiana-P.palmivora (ARI-tdTomato) timecourse
  - ``Orthogroups.csv`` - raw output from OrthoFinder, orthogroups between M.polymorpha and N.benthamiana
  - ``SingleCopyOrthogroups.txt`` - list of orthogroups with single-copy genes in M.polymorpha and N.benthamiana
  - ``Mpolymorpha_tidy_annotation.csv`` - summary tidy functional annotation for M.polymorpha

  
