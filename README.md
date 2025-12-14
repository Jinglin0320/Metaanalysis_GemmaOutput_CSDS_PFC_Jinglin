# Metaanalysis_GemmaOutput_CSDS_PFC_Jinglin
This repository contains the **analysis code supporting the manuscript**:https://www.biorxiv.org/content/10.1101/2025.10.23.683091v1

This project performs a **meta-analysis of chronic stress–associated gene expression changes in the mouse prefrontal cortex (PFC)** across publicly available transcriptomic datasets curated in **Gemma**.  

The analysis focuses primarily on two widely used rodent chronic stress paradigms:
- **Chronic Social Defeat Stress (CSDS)**
- **Chronic Unpredictable Mild Stress (CUMS)**


- **`DEextraction-CSDS+CUMS-PFC.R`**  
  Extracts differential expression results from selected Gemma datasets.

- **`Meta-Analysis of CSDS and CUMS.R`**  
  Performs the random-effects meta-analysis across CSDS and CUMS studies.


### Gene set enrichment analyses

- **`BrainGMT_Stress_PFC.R`**  
  Runs gene set enrichment analysis (fgsea) using curated gene sets.

- **`BrainGMT_Stress_PFC_SingleCell.R`**  
  Performs enrichment analysis using single-cell–derived brain cell-type gene sets.
