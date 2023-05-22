# Comprehensive Single-Cell Proteomics Analysis of Cancer Cell-Line Model Systems

This MSc thesis project is comprehensively investigates the proteomic profile of Melanoma, Hela, Pancreatic Ductal Adenocarcinoma (PDAC), and Monocytes (U937) cancer cell-lines.

## Table of Contents

* Project Overview
* Analysis Worklow
   *  Create Venn diagram
   *  Gene Set Enrichment Analysis (GSEA)
   *  KEGG pathway analysis
   *  Network analysis
   *  Cell Clustering: 
      *  Principal Component Analysis (PCA)
      *  Uniform Manifold Approximation and Projection (UMAP)
* License 
* Acknowledgments

## Project Overview

In this project, we gained a comprehensive understanding of Four cancer cell lines based on their proteome profiles at the level of single-cell. By integrating various analytics approaches, we identified Gene Ontology (GO) terms, the pathways associated with each cancer type, hub proteins, and clustering subpopulations.

## Analysis Workflow

The analysis workflow consists of the following pipeline:

## Create Venn diagram

A venn diagram is created to visualize the overlap of the identified proteins in four cancer cell-lines using [SCoPE](https://scope2.slavovlab.net) approaches. This diagram provides the intersected proteins that will be used for further analysis across different cancer types. You can check the code and data files [here](https://github.com/IbrahimElzahaby/MSc_Thesis/tree/58deca008814c9b84af54945a028f0797cf28bcb/Venn_Diagram).

## Gene Set Enrichment Analysis (GSEA)

GSEA is performed to identify the enriched gene sets and biological pathways associated with each cancer cell line. The enrichment analysis helped uncover the biological processes and pathways that plays a vital role in each cancer type. More details can be found in [Enrichment_analysis](https://github.com/IbrahimElzahaby/MSc_Thesis/tree/58deca008814c9b84af54945a028f0797cf28bcb/Enrichment_analysis) file.

## KEGG pathway analysis

[KEGG pathway analysis](https://github.com/IbrahimElzahaby/MSc_Thesis/blob/main/Enrichment_analysis/KEGG_pathway_analysis.R) is conducted to identify the key signaling pathways that are dysregulated in the different cancer cell lines. By analyzing the pathway-level protein expression, we can gain insights into the underlying biological mechanisms driving each cancer type.

## Network analysis

[Protein-Protein network analysis](https://github.com/IbrahimElzahaby/MSc_Thesis/blob/main/Network_analysis/Cor_Net.R) is performed to identify protein-protein correlations and construct functional protein networks for each cancer cell line. This analysis helps uncover the interconnectedness of proteins and potential hubs that play critical roles in Melanoma and Hela cancer types.

## Cell Clustering

In this [file](https://github.com/IbrahimElzahaby/MSc_Thesis/tree/58deca008814c9b84af54945a028f0797cf28bcb/Dimentionality_reduction), Cell clustering techniques are applied to the proteomic data to identify distinct cell populations within each cancer cell line. Two commonly used methods are utilized:

## Principal Component Analysis (PCA)

[Here](https://github.com/IbrahimElzahaby/MSc_Thesis/blob/main/Dimentionality_reduction/PCA_ALLCELLTYPES.R) we applied PCA to reduce the dimensionality of the data and visualize the variation between different cell populations within each cancer type. This technique helps identify clusters and patterns in the proteomic data.

## Uniform Manifold Approximation and Projection (UMAP)

UMAP is a dimensionality reduction technique that preserves the local structure of the data. It is used to visualize and identify distinct cell populations within each cancer cell line based on their proteomic profiles. you can find the code in this [R file](https://github.com/IbrahimElzahaby/MSc_Thesis/blob/main/Dimentionality_reduction/UMAP_ALLCELLTYPES.R)

## License

This project is released under the MIT License.

## Acknowledgments

I would like to acknowledge the [Slavov laboratory](https://slavovlab.net) where the data was obtained for this project. This project was conducted as part of the master's thesis at [Systems and Personalized Medicine Division - Laboratory of Systems and Synthetic Biology](https://www.wur.nl/en/Research-Results/Chair-groups/Agrotechnology-and-Food-Sciences/Laboratory-of-Systems-and-Synthetic-Biology.htm) at Wageningen University and Research. I extend my gratitude to my supervisors, [Dr. Cristina Furlan](https://www.wur.nl/en/persons/cristina-dr.-c-cristina-furlan.htm) and [Dr. Edoardo Saccenti](https://www.wur.nl/en/persons/edoardo-dr.-e-edoardo-saccenti.htm) for their guidance and support throughout the process of completing this thesis project.





























