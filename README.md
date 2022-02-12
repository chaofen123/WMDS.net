# WMDS.net

## Introduction

We developed a weighted MDS network model (WMDS.net) to find the drive nodes of differential gene co-expression networks. WMDS.net integrates the degree of nodes in the network and the significance of gene co-expression difference between two physiological states into the measurement of node controllability of the transcriptional network. To confirm its validity, we applied WMDS.net to the discovery of cancer driver genes in RNA-seq datasets from The Cancer Genome Atlas. WMDS.net is robust and powerful among various cancer datasets and outperformed the other top-tier tools with a better balance between precision and recall. 

## Gene_Net

Gene_Net is the hub network of intersection of gene interaction network constructed in Dawnrank[1] and the genes in TCGA.

## Tools
 
In this study, we implement the algorithm using the Matlab program language and solve the binary integer-programming problem using function “intlinprog” which is available in the Optimization ToolBox of MatLab version R2020b. The inputs to the main program are the expression matrices of the two conditions (which are TCGA_BLCAtumor and TCGA_BLCAnormal in the tested data), and the output is the gene list of key nodes in the final network. The rows of the expression matrices represent gene expression, and the columns represent samples. The row and column names of the expression matrices are gene symbols and sample ids, respectively.

## Reference

[1]Hou JP, Ma J: DawnRank: discovering personalized driver genes in cancer. Genome Medicine 2014, 6.
