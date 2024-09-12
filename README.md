# AWGE-ESPCA: An Edge Sparse PCA Model for Hermetia illucens Genomic Data Analysis



## Introduction

AWGE-ESPCA (Adaptive Noise Elimination Regularization and Weighted Gene Network - Edge Sparse PCA) is an unsupervised model designed for Hermetia illucens genome analysis. This model attempts to address three challenges in current research:

1. Limited available genomic data
2. Lack of AI feature selection models specific to insect genomes
3. Difficulty in selecting genes located in pathway enrichment regions

We hope this model may assist researchers in identifying potential biomarkers in insect genomic data analysis.





## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Simulation Experiment](#simulation-experiment)
- [File Structure](#file-structure)

## Installation

1. Clone the repository:

   ```
   git clone https://github.com/your-username/AWGE-ESPCA.git
   ```

2. Install required R packages:

   ```R
   install.packages(c("Matrix", "ggplot2", "pheatmap", "dplyr"))
   ```

## Usage

1. Set the working directory in the script:

   ```R
   setwd("/path/to/your/working/directory")
   ```

2. Run the main script:

   ```R
   source("Example.R")
   ```

3. The script will generate output files including `PC_data.csv` with principal component analysis results.

## Simulation Experiment

The Example.R script conducts a simulation experiment to demonstrate the potential advantages of AWGE-ESPCA over other PCA methods. It generates synthetic genomic data, constructs a gene network, and compares AWGE-ESPCA with PCA, SPCA, ESPCA, and DM-ESPCA. Results are saved in PC_data.csv for analysis.





## File Structure

```
AWGE-ESPCA/
│
├── Example.R              # Main script for data generation and analysis
├── functions/             # Directory containing function files
│   ├── fun_SPCA.R
│   ├── fun_ESPCA.R
│   ├── fun_DM_ESPCA.R
│   └── fun_AWGE_ESPCA.R
└── README.md
```

