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
