# Analysis of novel cfRNAs in the plasma of esophageal cancer and Barrett's ersophagus
## Alex D. Hill

### Abstract

Liquid biopsies detect disease noninvasively by profiling cell-free nucleic acids that are secreted into the circulation. However, existing methods exhibit low sensitivity for detecting early stages of diseases such as cancer. Here we show that long-read nanopore sequencing of long cell-free RNAs in plasma from healthy individuals, precancerous Barrett’s esophagus patients with high-grade dysplasia, or patients with esophageal adenocarcinoma reveals a diverse cell-free RNA transcriptome that can be leveraged for detecting and treating disease. We discovered 286,936 novel, intergenic cell-free RNAs, which we used to build a custom transcriptome reference for quantification, feature selection, and machine learning to accurately classify both precancer and cancer. Moreover, we found potential therapeutic targets, including metabolic, signaling, and immune checkpoint pathways, that are highly upregulated in both precancer and cancer patients. Our findings highlight the utility of our RNA liquid biopsy platform technology for discovering and targeting early stages of disease with molecular precision.

### Analysis

The analysis, panels, and figures generated in this project are all available and annotated in `notebooks/analysis.rmd` notebook. The RMarkdown file has a few necessary code chunks that read in the data and setup the environment, and all subsequent panels generate a single (or a few identical) panel(s). As long as the required libraries exist the notebooks should be able to be run front to back relatively easily.

#### TODO

 - [X] Add CRATER (as creater) submodule which contains a series of helper function in R for plotting and theming
 - [X] Add conda yaml for replicable environment

### Citation

For now the manuscript is published as a preprint on [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.07.02.662774v1).

Hill, A.D., *et al.* RNA liquid biopsy via nanopore sequencing for novel biomarker discovery and cancer early detection. *bioRxiv* (2025)
