# gsGNA: An Ensemble Pipeline for Mining Rice Abiotic Stress-Associated Regulators

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R 4.0+](https://img.shields.io/badge/R-4.0%2B-blue.svg)](https://www.r-project.org/)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-green.svg)](https://www.python.org/)

## Overview

**gsGNA** (**G**eneral **S**tress-resistance **G**ene Regulatory **N**etwork **A**nalysis) is an ensemble-based computational pipeline for inferring stress-responsive gene regulatory networks (GRNs) and identifying candidate regulatory factors in rice.

The pipeline integrates six base GRN inference algorithms and six ensemble fusion strategies, systematically evaluates different algorithm combinations, and selects suitable ensemble configurations through multi-layer validation. The resulting GRNs support downstream network characterization, functional module analysis, and pathway-guided candidate regulator screening.

The pipeline was developed and evaluated using **rice** (*Oryza sativa*) leaf RNA-seq data under three abiotic stress conditions: drought, alkaline, and cold stress. Three additional biotic-stress datasets were used to explore the application of the ensemble framework in different biological contexts.

## Key Features

- Integrates **six GRN inference algorithms**: GENIE3, KBoost, GRNBoost2, 3DCEMA, DeepRIG, and IGEGRNs.
- Compares **six ensemble fusion strategies**: AdaBoost, Bagging, Stacking, XGBoost, Hard Voting, and Quantile Rank Fusion (QRF).
- **Systematic algorithm evaluation**: Evaluates all 57 combinations of two to six base algorithms.
- **Multi-layer validation**: TF-wise nested cross-validation, independent TF-disjoint testing, and evaluation using ChIP-seq, TFBS, and TF perturbation-associated DEG evidence.
- **Network topology analysis**: Degree distribution, network connectivity, clustering coefficient, and small-world properties.
- **Hub TF identification**: Based on the top 10% of positive TF outdegree, with candidate-constrained randomization analysis.
- **Functional module detection**: Comparison of eight clustering algorithms, module reproducibility assessment, and GO/KEGG enrichment.
- **Cross-stress network analysis**: Directed-edge conservation, module preservation, TF–target rewiring, and stress-specific regulatory edge analysis.
- **Pathway-guided reverse screening**: Identification and prioritization of candidate upstream TFs using pathway-associated target genes and ensemble GRNs.

## Pipeline Workflow

<p align="center">
  <img src="figures/gsGNA_workflow.png" alt="gsGNA pipeline workflow" width="85%">
</p>

The gsGNA pipeline consists of four main stages:

1. **Input Data**: Gene expression matrices and species-specific TF lists from PlantTFDB.
2. **Base GRN Construction**: Six independent inference algorithms generate base GRNs.
3. **GRN Fusion and Evaluation**: Six ensemble strategies and 57 algorithm combinations are compared using TF-wise validation and independent testing, together with ChIP-seq, TFBS, and DEG evidence. XGBoost is selected for downstream ensemble GRN construction.
4. **Key TF Identification and Network Analysis**: Hub TF identification, network topology analysis, functional module detection, cross-stress network comparison, and pathway-guided candidate regulator screening.

## Dependencies

### Base GRN Inference Algorithms

The six base GRN inference algorithms are **not bundled** in this repository. Please install them separately, run them on your expression matrix with the TF list, and place the output edge lists (columns: `TF`, `Target`, `Weight`) into the corresponding method folders before running the ensemble fusion step.

| Algorithm | Language | Installation / Source | Reference |
|---|---|---|---|
| GENIE3 | R | `BiocManager::install("GENIE3")` ([link](https://bioconductor.org/packages/GENIE3)) | Huynh-Thu et al. 2010 |
| KBoost | R | `install.packages("KBoost")` ([link](https://bioconductor.org/packages/KBoost)) | Iglesias-Martinez et al. 2021 |
| GRNBoost2 | Python | `pip install arboreto` ([link](https://github.com/aertslab/arboreto)) | Moerman et al. 2019 |
| 3DCEMA | Python | [GitHub](https://github.com/YueFan1014/3DCEMA) | Fan & Ma, 2021 |
| DeepRIG | Python | [GitHub](https://github.com/JChander/DeepRIG) | Wang et al. 2023 |
| IGEGRNs | Python (PyTorch) | [GitHub](https://github.com/DHUDBlab/IGEGRNS) | Gan et al. 2024 |

Additional R and Python packages are required for ensemble learning, network analysis, visualization, and functional enrichment. Please refer to the individual scripts for package dependencies.

## Directory Structure

```text
gsGNA/
├── code/
│   ├── gsGNA_Ensemble_Comparison.R
│   │   # Six fusion strategies, 57 subsets,
│   │   # TF-wise CV, nested CV, and model selection
│   │
│   ├── XGB_Biotic_Stress.R
│   │   # Cross-context biotic-stress analysis
│   │
│   ├── Network validation/
│   │   # TFBS and DEG validation
│   │
│   ├── RNAseq_TF_analysis.R
│   ├── hub_TF_identification.R
│   ├── GRN_module_identification.R
│   ├── Abiotic_TF_Enrichment.R
│   ├── core_pathway_TF_screening.R
│   └── Abiotic_Radar_Plot.R
│
├── data/
│   ├── Osj_TF_list.txt
│   ├── chip_seq.csv
│   ├── BioMart_data.txt
│   └── ...
│
├── .gitignore
├── LICENSE
└── README.md
```

## Data Availability

The RNA-seq datasets used in this study are publicly available from the NCBI Gene Expression Omnibus (GEO).

| Dataset | Stress Type | GEO Accession |
|---|---|---|
| Rice leaf | Drought stress | [GSE121303](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121303) |
| Rice leaf | Alkaline stress | [GSE104928](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104928) |
| Rice leaf | Cold stress | [GSE112547](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112547) |
| *Magnaporthe oryzae* | Biotic (external application) | [GSE113553](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113553) |
| *Nephotettix cincticeps* feeding | Biotic (external application) | [GSE157400](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157400) |
| Rice stripe virus | Biotic (external application) | [GSE176497](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176497) |

The biotic-stress datasets were used to explore the applicability of the ensemble framework in additional biological contexts and were not treated as independent abiotic-stress validation datasets.

## Citation

Manuscript under preparation. Citation information and DOI will be updated upon publication.

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

## Contact

For questions or issues, please open an issue on GitHub or contact the corresponding authors:

- Xuemei Wei (wxm201@dali.edu.cn)
- Junpeng Zhang (zjp@dali.edu.cn)
