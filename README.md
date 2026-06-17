# gsGNA: A Universal Pipeline for Mining Plant Stress‑Resistance Regulators

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)

## Overview

gsGNA (General Stress‑resistance Gene Regulatory Network) is a universal computational pipeline for identifying plant stress‑resistance regulators using multi‑algorithm gene regulatory network (GRN) inference.

**Key features:**
- Integrates **six** GRN inference algorithms: GENIE3, KBoost, GRNBoost2, 3DCEMA, DeepRIG, IGEGRNS
- Supports multiple abiotic stress conditions (drought, alkaline, cold)
- Multi‑layer validation using ChIP‑seq, TFBS, and DEGs
- Network topology analysis (small‑world, scale‑free)
- Functional module detection and GO/KEGG enrichment
- Reverse screening for novel regulatory TFs

## Repository Structure
gsGNA/
├── code/ # All analysis scripts
│ ├── GRNBoost2网络评估.R
│ ├── GRNBoost筛选显著边.R
│ ├── IGE网络评估.R
│ ├── IGRGRN_Z分值筛选.R
│ ├── Kboost关键TF堆叠图.R
│ ├── RIG网络验证.R
│ ├── RNA-seq比较图.R
│ ├── TFBS和CHIP-seq共同验证.R
│ └── kboost关键TF富集分析.R
├── data/ # Input datasets
│ ├── GSE104928.xlsx # Drought stress
│ ├── GSE266657.xlsx # Alkaline stress
│ ├── GSE121303_Processed_data.xlsx # Cold stress
│ ├── Osj_TF_list.txt # Rice TF list
│ ├── chip_data.txt # ChIP‑seq validation data
│ └── BioMart_data.txt # Annotation data
└── results/ # Output results (to be added)
