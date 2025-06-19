# GBM Drug Discovery Pipeline

A comprehensive multi-omic data integration and machine learning pipeline for discovering novel drug treatments for Glioblastoma (GBM), the most aggressive form of brain cancer.

## Project Overview

This pipeline addresses the critical challenge of finding effective treatments for GBM by:
- Integrating multiple omics datasets (genomics, transcriptomics, methylation, clinical data)
- Predicting drug response using machine learning models
- Identifying drug repurposing opportunities from existing FDA-approved drugs
- Discovering biomarkers for personalized treatment strategies

## What is Multi-Omics?

Multi-omics combines different types of biological data:
- Genomics: DNA mutations and copy number variations
- Transcriptomics: Gene expression levels (RNA-seq)
- Methylation: DNA methylation patterns
- Clinical: Patient demographics, survival, treatment response
- Drug Sensitivity: How cancer cells respond to drugs in vitro

## Key Capabilities

### 1. Data Integration
- Downloads and processes TCGA GBM datasets
- Integrates GDSC drug sensitivity data
- Handles large-scale genomic data efficiently

### 2. Machine Learning Models
- Drug response prediction models
- Patient stratification algorithms
- Biomarker discovery tools
- Multi-omic data fusion techniques

### 3. Knowledge Graph
- Neo4j-based knowledge graph for biological relationships
- Integrates genes, patients, drugs, and pathways
- Enables complex biological queries

### 4. Drug Repurposing
- Identifies existing drugs that might work against GBM
- Reduces drug development time from 10-15 years to 2-3 years
- Significantly lowers development costs

## Prediction Goals

1. Drug Response Prediction: Forecast which drugs will work for specific patients
2. Drug Repurposing: Find FDA-approved drugs effective against GBM
3. Biomarker Discovery: Identify molecular markers for drug sensitivity
4. Patient Stratification: Group patients by likely treatment response

## Installation & Setup

### Prerequisites
- Python 3.10+ (recommended for compatibility)
- Git
- Conda (for environment management)

### Quick Start

1. Clone the repository
   ```bash
   git clone https://github.com/mohamedessamhashim/GBM.git
   cd GBM
   ```

2. Create conda environment
   ```bash
   conda create -n gbm_env python=3.10 -y
   conda activate gbm_env
   ```

3. Install dependencies
   ```bash
   pip install -r requirements.txt
   ```

4. Download data (optional - will be downloaded automatically when needed)
   ```bash
   python run_pipeline.py --download-data
   ```

## Project Structure

```
GBM/
├── src/                          # Source code
│   ├── data_acquisition/         # Data download and parsing
│   ├── preprocessing/            # Data cleaning and normalization
│   ├── feature_engineering/      # Feature creation and selection
│   ├── models/                   # Machine learning models
│   ├── knowledge_graph/          # Neo4j knowledge graph
│   └── evaluation/               # Model evaluation and metrics
├── data/                         # Data storage
│   ├── raw/                      # Raw downloaded data
│   ├── processed/                # Cleaned and processed data
│   └── external/                 # External datasets
├── notebooks/                    # Jupyter notebooks for exploration
├── results/                      # Model outputs and visualizations
├── tests/                        # Unit tests
└── docs/                         # Documentation
```

## Usage Examples

### Basic Pipeline Execution
```bash
# Run the complete pipeline
python run_pipeline.py

# Run specific components
python run_pipeline.py --data-only
python run_pipeline.py --train-models
python run_pipeline.py --predict-drugs
```

### Jupyter Notebooks
```bash
# Start Jupyter Lab
jupyter lab

# Or start Jupyter Notebook
jupyter notebook
```

### Data Acquisition
```python
from src.data_acquisition.tcga_downloader import TCGADownloader
from src.data_acquisition.gdsc_parser import GDSCDataParser

# Download TCGA data
tcga = TCGADownloader()
tcga.download_gbm_data()

# Parse GDSC data
gdsc = GDSCDataParser()
gdsc.download_and_parse_all()
```

### Knowledge Graph
```python
from src.knowledge_graph.graph_builder import GBMKnowledgeGraph

# Initialize knowledge graph
kg = GBMKnowledgeGraph()

# Add patient data
kg.add_patient("TCGA-02-0003", mutations=["TP53", "EGFR"], survival_days=365)

# Query for drug candidates
drugs = kg.find_drug_candidates_for_mutations(["TP53", "EGFR"])
```

## Hardware Requirements

### CPU-Only Setup (Recommended for beginners)
- Minimum: 8GB RAM, 4 CPU cores
- Recommended: 16GB RAM, 8 CPU cores
- Storage: 50GB free space for data and models

### GPU Setup (For advanced users)
- When needed: Training large deep learning models
- Recommended: NVIDIA GPU with 8GB+ VRAM
- Frameworks: TensorFlow/PyTorch with GPU support

## Expected Outputs

### 1. Drug Response Predictions
- Probability scores for drug effectiveness
- Patient-specific treatment recommendations
- Confidence intervals for predictions

### 2. Drug Repurposing Candidates
- List of FDA-approved drugs with high GBM activity
- Mechanism of action explanations
- Clinical trial recommendations

### 3. Biomarker Discovery
- Gene signatures for drug sensitivity
- Patient stratification markers
- Prognostic indicators

### 4. Knowledge Graph Insights
- Biological pathway analysis
- Drug-target interactions
- Patient similarity networks

## Research Applications

### For Researchers
- Drug Discovery: Identify novel drug candidates
- Clinical Trials: Design patient stratification strategies
- Biomarker Research: Discover predictive markers
- Mechanism Studies: Understand drug resistance

### For Clinicians
- Personalized Medicine: Tailor treatments to patient profiles
- Treatment Planning: Optimize drug combinations
- Prognosis: Predict patient outcomes

## Contributing

We welcome contributions! Please see our contributing guidelines:

1. Fork the repository
2. Create a feature branch
3. Add comprehensive tests
4. Submit a pull request

### Development Setup
```bash
# Install development dependencies
pip install -r requirements-dev.txt

# Run tests
pytest tests/

# Run linting
flake8 src/
```

## Documentation

- API Documentation: See docstrings in source code
- Tutorials: Check `notebooks/` directory
- Research Papers: Referenced in `docs/`

## Clinical Disclaimer

This pipeline is for research purposes only. Results should be validated through:
- In vitro experiments
- Animal studies
- Clinical trials
- Expert review

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- TCGA (The Cancer Genome Atlas) for genomic data
- GDSC (Genomics of Drug Sensitivity in Cancer) for drug response data
- Open-source bioinformatics community

## Contact

- Repository: https://github.com/mohamedessamhashim/GBM
- Issues: Use GitHub Issues for bug reports and feature requests
- Discussions: Use GitHub Discussions for questions and collaboration

---

Note: This pipeline represents a significant step toward personalized medicine for GBM patients. By combining computational biology with clinical insights, we aim to accelerate the discovery of effective treatments for this devastating disease. 