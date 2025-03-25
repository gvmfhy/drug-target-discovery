# Drug Target Discovery Pipeline

An automated pipeline for identifying potential drug targets from microarray gene expression data, focusing on prostate cancer.

## Overview

This pipeline processes gene expression data from the Gene Expression Omnibus (GEO), maps microarray probes to genes, performs differential expression analysis, constructs gene co-expression networks, and validates potential targets against drug databases.

Key features:
- Probe-to-gene mapping using Bioconductor
- Differential expression analysis between cancer and normal samples
- Network-based prioritization of potential drug targets
- Integration with OpenTargets API for drug target validation

## Data

This repository includes:
- A prostate cancer dataset (GSE46602) with 36 cancer and 14 normal samples
- Results from running the pipeline on this dataset
- Visualizations of the network and top targets

## Installation

### Prerequisites
- Python 3.8+
- R 4.0+ with Bioconductor

### Setup
1. Clone this repository:
```
git clone https://github.com/gvmfhy/drug-target-discovery.git
cd drug-target-discovery
```

2. Create and activate a virtual environment:
```
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install Python dependencies:
```
pip install -r requirements.txt
```

4. Install R dependencies:
```
Rscript install_r_packages.R
```

## Usage

Run the pipeline on the included prostate cancer dataset:
```
python drug_target_gse46602/pipeline2 --matrix-file drug_target_gse46602/prostate_GSE46602_series_matrix.txt
```

Results will be saved to a timestamped directory in the `results/` folder.

## Results

The pipeline generates:
- Network visualization of gene co-expression
- Ranked list of potential drug targets
- Validation against known drug associations

### Key Findings
- DIPK1B emerged as the top network hub (score: 0.957) with zero known drugs
- CACNA1F had 324 known drugs but a lower network centrality (score: 0.667)
- Several interferon genes showed moderate centrality with a handful of drugs

See the [project narrative](reports/project_narrative.md) for a detailed discussion of findings.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

This project was inspired by the growing accessibility of computational approaches to drug discovery, as highlighted in a LinkedIn post noting that "right now, today, you can download from among 4M whole genome expression studies all those that relate to a given disease, analyze the data, and identify a drug target in your living room, on a plane."
