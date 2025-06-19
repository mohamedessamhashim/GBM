# src/data_acquisition/tcga_downloader.py
import os
import requests
import tarfile
from pathlib import Path
import pandas as pd
import numpy as np
from tqdm import tqdm
import logging
import warnings
warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

"""
TCGA Downloader Module
----------------------
Handles downloading, extracting, and generating mock data for TCGA GBM datasets.
Used for multi-omic data integration in the GBM drug discovery pipeline.
"""

class TCGADownloader:
    def __init__(self, data_dir="data/raw"):
        """
        Initialize the downloader with the data directory.
        """
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
        # GDAC Firehose URLs (pre-processed TCGA data)
        self.base_url = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/"
        
        # GBM-specific datasets
        self.datasets = {
            'clinical': {
                'url': 'GBM/20160128/gdac.broadinstitute.org_GBM.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz',
                'description': 'Clinical data for GBM patients'
            },
            'expression': {
                'url': 'GBM/20160128/gdac.broadinstitute.org_GBM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz',
                'description': 'RNA-seq expression data'
            },
            'mutations': {
                'url': 'GBM/20160128/gdac.broadinstitute.org_GBM.Mutation_Packager_Calls.Level_3.2016012800.0.0.tar.gz',
                'description': 'Somatic mutation data'
            },
            'copy_number': {
                'url': 'GBM/20160128/gdac.broadinstitute.org_GBM.CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz',
                'description': 'Copy number variation data'
            },
            'methylation': {
                'url': 'GBM/20160128/gdac.broadinstitute.org_GBM.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0.tar.gz',
                'description': 'DNA methylation data'
            }
        }
        
    def download_gbm_data(self, force_download=False):
        """
        Download GBM datasets from GDAC Firehose.
        Pseudocode:
        1. For each data type (clinical, expression, etc.):
            a. Check if file exists.
            b. If not, download and extract.
            c. If download fails, create mock data.
        2. Return dictionary of downloaded files.
        """
        logger.info("Starting TCGA GBM data download...")
        
        downloaded_files = {}
        
        for data_type, dataset_info in self.datasets.items():
            url = self.base_url + dataset_info['url']
            filename = self.data_dir / f"tcga_gbm_{data_type}.tar.gz"
            
            if filename.exists() and not force_download:
                logger.info(f"✓ {data_type} data already exists: {filename}")
                downloaded_files[data_type] = filename
                continue
            
            try:
                logger.info(f"Downloading {data_type} data: {dataset_info['description']}")
                self._download_with_progress(url, filename)
                
                # Extract archive
                extract_dir = self.data_dir / f"tcga_gbm_{data_type}"
                extract_dir.mkdir(exist_ok=True)
                
                logger.info(f"Extracting {data_type} data...")
                with tarfile.open(filename, 'r:gz') as tar:
                    tar.extractall(extract_dir)
                    
                downloaded_files[data_type] = filename
                logger.info(f"✓ {data_type} data downloaded and extracted")
                
            except Exception as e:
                logger.error(f"Failed to download {data_type} data: {e}")
                logger.info(f"Will create mock data for {data_type}")
                downloaded_files[data_type] = self._create_mock_data(data_type)
        
        return downloaded_files
    
    def _download_with_progress(self, url, filename):
        """
        Download file with progress bar.
        Pseudocode:
        1. Stream file from URL.
        2. Write to disk in chunks, updating progress bar.
        3. Handle download errors.
        """
        try:
            response = requests.get(url, stream=True, timeout=30)
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            
            with open(filename, 'wb') as file, tqdm(
                desc=filename.name,
                total=total_size,
                unit='B',
                unit_scale=True,
                unit_divisor=1024,
            ) as pbar:
                for data in response.iter_content(chunk_size=8192):
                    size = file.write(data)
                    pbar.update(size)
                    
        except requests.exceptions.RequestException as e:
            logger.error(f"Download failed for {url}: {e}")
            raise
    
    def _create_mock_data(self, data_type):
        """
        Create realistic mock data for development.
        Pseudocode:
        1. Depending on data_type, call the appropriate mock data generator.
        2. Return path to mock data file.
        """
        logger.info(f"Creating mock {data_type} data...")
        
        if data_type == 'clinical':
            return self._create_mock_clinical_data()
        elif data_type == 'expression':
            return self._create_mock_expression_data()
        elif data_type == 'mutations':
            return self._create_mock_mutation_data()
        elif data_type == 'copy_number':
            return self._create_mock_cnv_data()
        elif data_type == 'methylation':
            return self._create_mock_methylation_data()
        else:
            logger.warning(f"Unknown data type: {data_type}")
            return None
    
    def _create_mock_clinical_data(self):
        """
        Create realistic mock clinical data.
        Pseudocode:
        1. Generate synthetic patient data with realistic distributions.
        2. Correlate clinical features with survival.
        3. Save to CSV and return path.
        """
        np.random.seed(42)
        
        n_patients = 200
        
        # Generate realistic GBM patient data
        clinical_data = {
            'patient_id': [f'TCGA-{i:02d}-{j:04d}' for i in range(1, 21) for j in range(1, 11)][:n_patients],
            'age_at_diagnosis': np.random.normal(60, 15, n_patients).astype(int),
            'gender': np.random.choice(['Male', 'Female'], n_patients, p=[0.6, 0.4]),
            'vital_status': np.random.choice(['Alive', 'Dead'], n_patients, p=[0.2, 0.8]),
            'overall_survival_days': np.random.exponential(365, n_patients).astype(int),
            'days_to_death': np.random.exponential(300, n_patients).astype(int),
            'days_to_last_followup': np.random.exponential(400, n_patients).astype(int),
            'mgmt_methylation_status': np.random.choice(['Methylated', 'Unmethylated', 'Unknown'], n_patients, p=[0.4, 0.5, 0.1]),
            'idh_mutation_status': np.random.choice(['Mutant', 'Wild-type', 'Unknown'], n_patients, p=[0.1, 0.8, 0.1]),
            'egfr_amplification': np.random.choice(['Amplified', 'Not Amplified', 'Unknown'], n_patients, p=[0.4, 0.5, 0.1]),
            'pten_mutation': np.random.choice(['Mutant', 'Wild-type', 'Unknown'], n_patients, p=[0.3, 0.6, 0.1]),
            'tp53_mutation': np.random.choice(['Mutant', 'Wild-type', 'Unknown'], n_patients, p=[0.3, 0.6, 0.1]),
            'treatment_response': np.random.choice(['Complete Response', 'Partial Response', 'Stable Disease', 'Progressive Disease'], n_patients),
            'karnofsky_performance_score': np.random.randint(30, 100, n_patients),
            'tumor_grade': np.random.choice(['Grade IV'], n_patients),  # GBM is always Grade IV
            'tumor_size_cm': np.random.uniform(1, 8, n_patients).round(1),
            'extent_of_resection': np.random.choice(['Gross Total', 'Subtotal', 'Biopsy Only'], n_patients, p=[0.4, 0.4, 0.2])
        }
        
        # Create correlations between variables
        # MGMT methylation correlates with better survival
        mgmt_methylated = clinical_data['mgmt_methylation_status'] == 'Methylated'
        clinical_data['overall_survival_days'][mgmt_methylated] = np.random.exponential(500, mgmt_methylated.sum()).astype(int)
        
        # IDH mutation correlates with better survival
        idh_mutant = clinical_data['idh_mutation_status'] == 'Mutant'
        clinical_data['overall_survival_days'][idh_mutant] = np.random.exponential(600, idh_mutant.sum()).astype(int)
        
        df = pd.DataFrame(clinical_data)
        
        # Save mock clinical data
        output_file = self.data_dir / "tcga_gbm_clinical_mock.csv"
        df.to_csv(output_file, index=False)
        logger.info(f"✓ Created mock clinical data: {output_file}")
        
        return output_file
    
    def _create_mock_expression_data(self):
        """
        Create realistic mock RNA-seq expression data.
        Pseudocode:
        1. Generate synthetic gene expression matrix.
        2. Adjust expression for clinical features.
        3. Save to CSV and return path.
        """
        np.random.seed(42)
        
        n_patients = 200
        n_genes = 20000
        
        # Generate realistic gene names
        gene_names = [f'GENE_{i:05d}' for i in range(n_genes)]
        
        # Generate log-normal expression values (typical for RNA-seq)
        # Most genes have low expression, few have high expression
        expression_values = np.random.lognormal(mean=3, sigma=2, size=(n_genes, n_patients))
        
        # Add some genes with differential expression based on clinical features
        # Load clinical data to get patient characteristics
        clinical_file = self.data_dir / "tcga_gbm_clinical_mock.csv"
        if clinical_file.exists():
            clinical_df = pd.read_csv(clinical_file)
            
            # MGMT-associated genes (first 100 genes)
            mgmt_status = clinical_df['mgmt_methylation_status'].values
            for i in range(100):
                methylated_mask = mgmt_status == 'Methylated'
                unmethylated_mask = mgmt_status == 'Unmethylated'
                
                if methylated_mask.sum() > 0:
                    expression_values[i, methylated_mask] *= np.random.uniform(1.5, 3.0)
                if unmethylated_mask.sum() > 0:
                    expression_values[i, unmethylated_mask] *= np.random.uniform(0.3, 0.7)
            
            # IDH-associated genes (next 100 genes)
            idh_status = clinical_df['idh_mutation_status'].values
            for i in range(100, 200):
                mutant_mask = idh_status == 'Mutant'
                wildtype_mask = idh_status == 'Wild-type'
                
                if mutant_mask.sum() > 0:
                    expression_values[i, mutant_mask] *= np.random.uniform(2.0, 4.0)
                if wildtype_mask.sum() > 0:
                    expression_values[i, wildtype_mask] *= np.random.uniform(0.2, 0.6)
        
        # Create DataFrame
        patient_ids = [f'TCGA-{i:02d}-{j:04d}' for i in range(1, 21) for j in range(1, 11)][:n_patients]
        
        expr_df = pd.DataFrame(
            expression_values,
            index=gene_names,
            columns=patient_ids
        )
        
        # Save mock expression data
        output_file = self.data_dir / "tcga_gbm_expression_mock.csv"
        expr_df.to_csv(output_file)
        logger.info(f"✓ Created mock expression data: {output_file}")
        
        return output_file
    
    def _create_mock_mutation_data(self):
        """
        Create realistic mock mutation data.
        Pseudocode:
        1. Generate synthetic mutation events for each patient.
        2. Save to CSV and return path.
        """
        np.random.seed(42)
        
        # Common GBM driver genes
        driver_genes = [
            'EGFR', 'TP53', 'PTEN', 'IDH1', 'ATRX', 'CIC', 'FUBP1', 'NF1', 
            'PIK3CA', 'PIK3R1', 'RB1', 'CDKN2A', 'CDKN2B', 'CDKN2C',
            'TERT', 'BRAF', 'KRAS', 'NRAS', 'HRAS', 'MAP2K1', 'MAP2K2'
        ]
        
        # Generate mutations
        mutations = []
        n_patients = 200
        patient_ids = [f'TCGA-{i:02d}-{j:04d}' for i in range(1, 21) for j in range(1, 11)][:n_patients]
        
        for patient_id in patient_ids:
            # Each patient has 5-15 mutations
            n_mutations = np.random.randint(5, 16)
            
            for _ in range(n_mutations):
                gene = np.random.choice(driver_genes)
                mutation_type = np.random.choice(['Missense', 'Nonsense', 'Frameshift', 'Splice_Site'])
                
                # Generate realistic mutation positions
                if gene == 'EGFR':
                    position = np.random.randint(1, 1210)
                elif gene == 'TP53':
                    position = np.random.randint(1, 393)
                elif gene == 'PTEN':
                    position = np.random.randint(1, 403)
                else:
                    position = np.random.randint(1, 1000)
                
                mutations.append({
                    'patient_id': patient_id,
                    'gene': gene,
                    'mutation_type': mutation_type,
                    'position': position,
                    'chromosome': np.random.randint(1, 23),
                    'start_position': position * 1000,
                    'end_position': position * 1000 + 1,
                    'reference_allele': np.random.choice(['A', 'T', 'C', 'G']),
                    'tumor_seq_allele1': np.random.choice(['A', 'T', 'C', 'G']),
                    'tumor_seq_allele2': np.random.choice(['A', 'T', 'C', 'G']),
                    'variant_classification': mutation_type,
                    'variant_type': 'SNP',
                    'impact': np.random.choice(['HIGH', 'MODERATE', 'LOW']),
                    'allele_frequency': np.random.uniform(0.1, 1.0)
                })
        
        df = pd.DataFrame(mutations)
        
        # Save mock mutation data
        output_file = self.data_dir / "tcga_gbm_mutations_mock.csv"
        df.to_csv(output_file, index=False)
        logger.info(f"✓ Created mock mutation data: {output_file}")
        
        return output_file
    
    def _create_mock_cnv_data(self):
        """
        Create realistic mock copy number variation data.
        Pseudocode:
        1. Generate synthetic CNV events for each patient.
        2. Save to CSV and return path.
        """
        np.random.seed(42)
        
        # Generate CNV data
        cnv_data = []
        n_patients = 200
        patient_ids = [f'TCGA-{i:02d}-{j:04d}' for i in range(1, 21) for j in range(1, 11)][:n_patients]
        
        # Common GBM genes with CNV
        cnv_genes = ['EGFR', 'CDKN2A', 'CDKN2B', 'PTEN', 'RB1', 'MDM2', 'MDM4', 'CDK4', 'CDK6']
        
        for patient_id in patient_ids:
            # Each patient has 2-8 CNV events
            n_cnv = np.random.randint(2, 9)
            
            for _ in range(n_cnv):
                gene = np.random.choice(cnv_genes)
                cnv_type = np.random.choice(['Amplification', 'Deletion'], p=[0.7, 0.3])
                
                cnv_data.append({
                    'patient_id': patient_id,
                    'gene': gene,
                    'cnv_type': cnv_type,
                    'copy_number': np.random.randint(0, 10) if cnv_type == 'Amplification' else np.random.randint(0, 2),
                    'chromosome': np.random.randint(1, 23),
                    'start_position': np.random.randint(1, 100000000),
                    'end_position': np.random.randint(1, 100000000),
                    'confidence': np.random.uniform(0.5, 1.0)
                })
        
        df = pd.DataFrame(cnv_data)
        
        # Save mock CNV data
        output_file = self.data_dir / "tcga_gbm_cnv_mock.csv"
        df.to_csv(output_file, index=False)
        logger.info(f"✓ Created mock CNV data: {output_file}")
        
        return output_file
    
    def _create_mock_methylation_data(self):
        """
        Create realistic mock methylation data.
        Pseudocode:
        1. Generate synthetic methylation beta values for each patient.
        2. Save to CSV and return path.
        """
        np.random.seed(42)
        
        # Generate methylation data
        methylation_data = []
        n_patients = 200
        patient_ids = [f'TCGA-{i:02d}-{j:04d}' for i in range(1, 21) for j in range(1, 11)][:n_patients]
        
        # Common methylation sites
        methylation_sites = [f'cg{i:07d}' for i in range(1000)]
        
        for patient_id in patient_ids:
            for site in methylation_sites:
                # Beta values range from 0 to 1
                beta_value = np.random.beta(2, 2)
                
                methylation_data.append({
                    'patient_id': patient_id,
                    'probe_id': site,
                    'beta_value': beta_value,
                    'chromosome': np.random.randint(1, 23),
                    'position': np.random.randint(1, 100000000),
                    'gene': np.random.choice(['EGFR', 'MGMT', 'PTEN', 'TP53', 'IDH1', 'Unknown'])
                })
        
        df = pd.DataFrame(methylation_data)
        
        # Save mock methylation data
        output_file = self.data_dir / "tcga_gbm_methylation_mock.csv"
        df.to_csv(output_file, index=False)
        logger.info(f"✓ Created mock methylation data: {output_file}")
        
        return output_file
    
    def load_processed_data(self):
        """
        Load all processed TCGA data.
        Pseudocode:
        1. Load each processed data file if it exists.
        2. Return a dictionary of DataFrames.
        """
        logger.info("Loading processed TCGA data...")
        
        data = {}
        
        # Load clinical data
        clinical_file = self.data_dir / "tcga_gbm_clinical_mock.csv"
        if clinical_file.exists():
            data['clinical'] = pd.read_csv(clinical_file)
            logger.info(f"✓ Loaded clinical data: {data['clinical'].shape}")
        
        # Load expression data
        expression_file = self.data_dir / "tcga_gbm_expression_mock.csv"
        if expression_file.exists():
            data['expression'] = pd.read_csv(expression_file, index_col=0)
            logger.info(f"✓ Loaded expression data: {data['expression'].shape}")
        
        # Load mutation data
        mutation_file = self.data_dir / "tcga_gbm_mutations_mock.csv"
        if mutation_file.exists():
            data['mutations'] = pd.read_csv(mutation_file)
            logger.info(f"✓ Loaded mutation data: {data['mutations'].shape}")
        
        # Load CNV data
        cnv_file = self.data_dir / "tcga_gbm_cnv_mock.csv"
        if cnv_file.exists():
            data['cnv'] = pd.read_csv(cnv_file)
            logger.info(f"✓ Loaded CNV data: {data['cnv'].shape}")
        
        # Load methylation data
        methylation_file = self.data_dir / "tcga_gbm_methylation_mock.csv"
        if methylation_file.exists():
            data['methylation'] = pd.read_csv(methylation_file)
            logger.info(f"✓ Loaded methylation data: {data['methylation'].shape}")
        
        return data

# Usage example
if __name__ == "__main__":
    downloader = TCGADownloader()
    
    # Download or create mock data
    downloaded_files = downloader.download_gbm_data()
    
    # Load processed data
    tcga_data = downloader.load_processed_data()
    
    print(f"Loaded {len(tcga_data)} datasets")
    for data_type, df in tcga_data.items():
        print(f"  {data_type}: {df.shape}") 