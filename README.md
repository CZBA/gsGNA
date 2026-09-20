# gsGNA: A General Pipeline for Mining Plant Stress-Resistance Regulators

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R 4.0+](https://img.shields.io/badge/R-4.0%2B-blue.svg)](https://www.r-project.org/)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-green.svg)](https://www.python.org/)
[![DOI: pending](https://img.shields.io/badge/DOI-pending-lightgrey.svg)]()

## Overview

**gsGNA** (**G**eneral **S**tress-resistance **G**ene Regulatory **N**etwork **A**nalysis) is a universal computational pipeline for identifying plant stress-resistance regulators via multi-algorithm ensemble gene regulatory network (GRN) inference. The pipeline integrates six base GRN inference algorithms and six ensemble fusion strategies, selects the optimal strategy via multi-layer validation, and constructs high-confidence GRNs for downstream regulator mining.

The pipeline was developed and benchmarked on **rice** (*Oryza sativa*) leaf RNA-seq data under three abiotic stresses (drought, alkaline, cold), with external validation on three biotic stress datasets.

## Key features

- Integrates **six GRN inference algorithms**: GENIE3, KBoost, GRNBoost2, 3DCEMA, DeepRIG, IGEGRNS
- Compares **six ensemble fusion strategies**: AdaBoost, Bagging, Stacking, XGBoost, hard voting, quantile rank fusion
- **Multi-layer validation** using ChIP-seq (AUROC/AUPR), TFBS overlap, and TF overexpression/interference DEGs
- **Network topology analysis** (scale-free, small-world properties)
- **Functional module detection** (8 clustering algorithms) and GO/KEGG enrichment
- **Key TF identification** via degree centrality with permutation-based significance testing
- **Reverse screening** for novel regulatory TFs based on core pathway target genes

## Pipeline Workflow

<img width="693" height="853" alt="image" src="https://github.com/user-attachments/assets/22836793-a221-4567-800d-bfeb52d8e9c2" />

The gsGNA pipeline consists of four main stages:

1. **Input Data**: Gene expression matrices + TF list (from PlantTFDB)
2. **Base GRN Construction**: Six independent inference algorithms generate base GRNs
3. **GRN Fusion and Evaluation**: Six ensemble strategies are compared and validated using ChIP-seq, TFBS, and DEG data; XGBoost is selected as the optimal strategy
4. **Key TF Identification & Network Analysis**: Topological analysis, module detection, functional enrichment, and reverse screening for novel regulators

## Dependencies

### Base GRN inference algorithms

The six base GRN inference algorithms are **not bundled** in this repository. Please install them separately, run them on your expression matrix with the TF list, and place the output edge lists (columns: `TF`, `Target`, `Weight`) into the corresponding method folders before running the ensemble fusion step.

| Algorithm | Language | Installation / Source | Reference |
|---|---|---|---|
| GENIE3 | R | `BiocManager::install("GENIE3")` ([link](https://bioconductor.org/packages/GENIE3)) | Huynh-Thu et al. 2010 |
| KBoost | R | `install.packages("KBoost")` ([link](https://bioconductor.org/packages/KBoost)) | Iglesias-Martinez et al. 2021 |
| GRNBoost2 | Python | `pip install arboreto` ([link](https://github.com/aertslab/arboreto)) | Moerman et al. 2019 |
| 3DCEMA | Python | [github.com/YueFan1014/3DCEMA](https://github.com/YueFan1014/3DCEMA) | Fan & Ma. 2021 |
| DeepRIG | Python | [github.com/JChander/DeepRIG](https://github.com/JChander/DeepRIG) | Wang et al. 2023 |
| IGEGRNS | Python (PyTorch) | [github.com/DHUDBlab/IGEGRNS](https://github.com/DHUDBlab/IGEGRNS) | Gan et al. 2024 |


## Directory structure

```
gsGNA/
├── code/
│   ├── ensemble method/          # Six ensemble fusion strategies
│   │   ├── AdaBoost.R
│   │   ├── Bagging.R
│   │   ├── Stacking.R
│   │   ├── XGBoost.R             # Recommended optimal strategy
│   │   ├── Hard_Voting.R
│   │   ├── Quantile_Rank_Fusion_QRF.R
│   │   ├── Evaluate_*.R          # Evaluation scripts for each strategy
│   │   ├── Evaluation_QRF.R
│   │   └── XGB_Biotic_Stress.R   # External validation on biotic stress
│   ├── Network validation/       # TFBS and DEG overlap validation
│   ├── RNAseq_TF_analysis.R      # Data preprocessing and expression analysis
│   ├── hub_TF_identification.R   # Key TF identification via degree centrality
│   ├── GRN_module_identification.R  # Module detection (8 algorithms)
│   ├── Abiotic_TF_Enrichment.R   # GO/KEGG enrichment analysis
│   ├── core_pathway_TF_screening.R  # Reverse screening for novel regulators
│   └── Abiotic_Radar_Plot.R      # Visualization
├── data/
│   ├── Osj_TF_list.txt           # Rice TF annotation (PlantTFDB v5.0)
│   ├── chip_data.txt             # ChIP-seq validated interactions
│   ├── BioMart_data.txt          # Gene annotation
│   └── ...                       # Example expression datasets
├── .gitignore
├── LICENSE
└── README.md
```

## Data availability

The RNA-seq datasets used in this study are publicly available from the NCBI Gene Expression Omnibus (GEO):

| Dataset | Stress type | GEO accession |
|---|---|---|
| Drought | Abiotic | [GSE121303](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121303) |
| Alkaline | Abiotic | [GSE104928](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104928) |
| Cold | Abiotic | [GSE112547](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112547) |
| *Magnaporthe oryzae* | Biotic (external validation) | [GSE113553](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113553) |
| *Nephotettix cincticeps* feeding | Biotic (external validation) | [GSE157400](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157400) |
| Rice stripe virus | Biotic (external validation) | [GSE176497](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176497) |

## Citation

Manuscript under preparation. Citation information will be updated upon submission.

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

## Contact

For questions or issues, please open an issue on GitHub or contact the corresponding authors:

- Xuemei Wei (wxm201@dali.edu.cn)
- Junpeng Zhang (zjp@dali.edu.cn)
