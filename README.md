# COVID-19 Multi-Level Transcriptomic Analysis: Differential Expression, Network Biology, and Machine Learning

## Project Overview

This project presents a comprehensive transcriptomic analysis of COVID-19 patient gene expression data using RNA-seq workflows and systems biology approaches.

The study extends beyond conventional differential expression analysis by integrating:

* Differential Gene Expression Analysis (DESeq2)
* Functional Enrichment Analysis (GO & KEGG)
* Weighted Gene Co-expression Network Analysis (WGCNA)
* Protein-Protein Interaction (PPI) Network Analysis
* Disease Severity Stratification
* Cross-Disease Comparative Analysis
* Machine Learning-Based Classification
* Multi-Method Validation (DESeq2, edgeR, limma-voom)

## Objectives

1. Identify genes dysregulated in COVID-19 patients.
2. Discover biological pathways associated with infection.
3. Detect co-expression modules linked to disease severity.
4. Identify hub genes through network analysis.
5. Compare COVID-19 signatures with other respiratory infections.
6. Evaluate the robustness of findings using multiple statistical frameworks.
7. Develop predictive models capable of classifying disease status.

## Workflow

Data Acquisition
↓
Quality Control
↓
Normalization
↓
Differential Expression Analysis
↓
GO & KEGG Enrichment
↓
WGCNA Module Detection
↓
PPI Network Construction
↓
Hub Gene Identification
↓
Severity-Level Comparison
↓
Cross-Disease Analysis
↓
Machine Learning Classification
↓
Biological Interpretation

## Key Analyses

### Differential Expression Analysis

* DESeq2
* edgeR
* limma-voom

### Functional Enrichment

* Gene Ontology (GO)
* KEGG Pathway Analysis

### WGCNA

* Identification of co-expression modules
* Correlation of modules with disease severity
* Detection of biologically relevant gene clusters

### Protein-Protein Interaction Analysis

* STRING Database
* Cytoscape Network Visualization
* Hub Gene Identification

### Disease Severity Analysis

* Mild vs Severe
* Severe vs ICU
* Identification of severity-associated molecular signatures

### Cross-Disease Comparison

* COVID-19 vs Influenza
* COVID-19 vs RSV
* Identification of disease-specific biomarkers

### Machine Learning

* Random Forest
* Support Vector Machine (SVM)

Evaluation Metrics:

* Accuracy
* Precision
* Recall
* F1 Score
* ROC-AUC

## Major Findings

* Immune and inflammatory pathways were strongly dysregulated in COVID-19.
* Multiple co-expression modules showed significant correlation with disease severity.
* Hub genes identified through PPI analysis may represent potential biomarkers.
* Comparative analysis revealed signatures unique to COVID-19 relative to other respiratory infections.
* Machine learning models successfully distinguished patients from controls using transcriptomic signatures.
* Cross-validation with DESeq2, edgeR, and limma-voom demonstrated strong concordance of results.

## Technologies Used

* R
* DESeq2
* edgeR
* limma
* WGCNA
* clusterProfiler
* STRINGdb
* Cytoscape
* Random Forest
* e1071 (SVM)
* ggplot2

## Repository Structure

data/
scripts/
results/
figures/
notebooks/
README.md

## Future Directions

* Multi-omics integration
* Single-cell RNA-seq analysis
* Biomarker validation on external cohorts
* Drug repurposing analysis using hub genes
