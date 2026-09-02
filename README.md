# gsGNA: A General Pipeline for Mining Plant Stress-Resistance Regulators

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.8+-green.svg)](https://www.python.org/)
[![DOI](https://img.shields.io/badge/DOI-pending-blue.svg)]()

## Overview

gsGNA (General Stress-resistance Gene Regulatory Network) is a computational pipeline for identifying plant stress-resistance regulators via multi-algorithm ensemble gene regulatory network (GRN) inference. The pipeline integrates six base GRN inference algorithms and six ensemble fusion strategies, selects the optimal strategy via multi-layer validation, and constructs high-confidence GRNs for downstream regulator mining.

The pipeline was developed and benchmarked on **rice (Oryza sativa)** leaf RNA-seq data under three abiotic stresses (drought, alkaline, cold), with external validation on three biotic stress datasets.

**Key features:**
- Integrates **six** GRN inference algorithms: GENIE3, KBoost, GRNBoost2, 3DCEMA, DeepRIG, IGEGRNS
- Compares **six** ensemble fusion strategies: AdaBoost, Bagging, Stacking, XGBoost, hard voting, quantile rank fusion
- Multi-layer validation using ChIP-seq (AUROC/AUPR), TFBS overlap, and TF overexpression/interference DEGs
- Network topology analysis (scale-free, small-world)
- Functional module detection (8 clustering algorithms) and GO/KEGG enrichment
- Key TF identification via degree centrality with permutation-based significance testing
- Reverse screening for novel regulatory TFs based on core pathway target genes

## Pipeline Workflow
                          <img width="693" height="853" alt="image" src="https://github.com/user-attachments/assets/007cd571-a00d-41a6-af21-8508f4fe8b5f" />

