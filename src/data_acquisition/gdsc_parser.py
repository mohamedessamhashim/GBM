# src/data_acquisition/gdsc_parser.py
import pandas as pd
import numpy as np
from pathlib import Path
import requests
import warnings
from tqdm import tqdm
import logging
import zipfile
import io

warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

"""
GDSC Data Parser Module
----------------------
Handles downloading, parsing, and generating mock data for GDSC datasets.
Used for drug sensitivity and cell line data in the GBM drug discovery pipeline.
"""

class GDSCDataParser:
    def __init__(self, data_dir="data/raw/gdsc"):
        """
        Initialize the parser with the data directory.
        """
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
        # GDSC file URLs (these may change - check GDSC website for latest)
        self.urls = {
            'drug_response': 'https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/preprocessed/Cell_line_RMA_proc_basalExp.txt.zip',
            'drug_info': 'https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/suppData/TableS1E.xlsx',
            'cell_line_info': 'https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/suppData/TableS1F.xlsx',
            'sensitivity_data': 'https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/suppData/TableS4A.xlsx'
        }
        
        # FDA-approved drugs commonly used in cancer treatment
        self.fda_approved_drugs = [
            'Temozolomide', 'Bevacizumab', 'Lomustine', 'Carmustine', 'Procarbazine',
            'Vincristine', 'Irinotecan', 'Erlotinib', 'Gefitinib', 'Sorafenib',
            'Sunitinib', 'Pazopanib', 'Vandetanib', 'Cabozantinib', 'Regorafenib',
            'Dabrafenib', 'Vemurafenib', 'Trametinib', 'Cobimetinib', 'Binimetinib',
            'Olaparib', 'Rucaparib', 'Niraparib', 'Talazoparib', 'Veliparib',
            'Pembrolizumab', 'Nivolumab', 'Atezolizumab', 'Durvalumab', 'Avelumab',
            'Ipilimumab', 'Cisplatin', 'Carboplatin', 'Oxaliplatin', 'Doxorubicin',
            'Epirubicin', 'Paclitaxel', 'Docetaxel', 'Gemcitabine', '5-Fluorouracil',
            'Capecitabine', 'Methotrexate', 'Cyclophosphamide', 'Ifosfamide',
            'Etoposide', 'Topotecan', 'Irinotecan', 'Mitomycin', 'Bleomycin',
            'Vinblastine', 'Vinorelbine', 'Dactinomycin', 'Daunorubicin'
        ]
    
    def download_and_parse_all(self):
        """
        Download and parse all GDSC datasets.
        Pseudocode:
        1. Try to download and parse sensitivity data.
        2. If successful, filter for brain cell lines and save.
        3. If not, create mock data.
        4. Return parsed data dictionary.
        """
        logger.info("Starting GDSC data download and parsing...")
        
        parsed_data = {}
        
        try:
            # Download drug sensitivity data
            logger.info("Processing drug sensitivity data...")
            sensitivity_df = self._get_sensitivity_data()
            
            if sensitivity_df is not None:
                # Filter for CNS/brain cancer cell lines
                brain_lines = sensitivity_df[
                    sensitivity_df['GDSC Tissue descriptor 1'] == 'central_nervous_system'
                ]
                
                logger.info(f"Found {len(brain_lines)} brain cancer cell line experiments")
                
                # Get unique cell lines and drugs
                unique_cell_lines = brain_lines['Cell line name'].unique()
                unique_drugs = brain_lines['Drug name'].unique()
                
                logger.info(f"Unique cell lines: {len(unique_cell_lines)}")
                logger.info(f"Unique drugs tested: {len(unique_drugs)}")
                
                parsed_data['brain_sensitivity'] = brain_lines
                parsed_data['cell_lines'] = unique_cell_lines
                parsed_data['drugs'] = unique_drugs
                
                # Save processed data
                self._save_processed_data(parsed_data)
                
            else:
                logger.warning("Could not download real GDSC data, creating mock data")
                parsed_data = self._create_mock_gdsc_data()
                
        except Exception as e:
            logger.error(f"Error processing GDSC data: {e}")
            logger.info("Creating mock GDSC data for development")
            parsed_data = self._create_mock_gdsc_data()
        
        return parsed_data
    
    def _get_sensitivity_data(self):
        """
        Parse drug sensitivity data with error handling.
        Pseudocode:
        1. Try to load from local file.
        2. If not, download from URL and save locally.
        3. Return DataFrame or None on error.
        """
        try:
            # Try to read from local file first
            local_file = self.data_dir / "sensitivity_data.xlsx"
            if local_file.exists():
                logger.info("Loading local sensitivity data...")
                return pd.read_excel(local_file)
            
            # Download from URL
            logger.info("Downloading sensitivity data from GDSC...")
            url = self.urls['sensitivity_data']
            
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            # Read Excel file from memory
            df = pd.read_excel(io.BytesIO(response.content), engine='openpyxl')
            
            # Save locally for future use
            df.to_excel(local_file, index=False)
            logger.info(f"✓ Downloaded and saved sensitivity data: {df.shape}")
            
            return df
            
        except Exception as e:
            logger.error(f"Error downloading GDSC data: {e}")
            return None
    
    def _create_mock_gdsc_data(self):
        """
        Create comprehensive mock GDSC data for development.
        Pseudocode:
        1. Generate synthetic cell lines and drugs.
        2. Generate mock sensitivity, drug info, and expression data.
        3. Save all mock data.
        4. Return mock data dictionary.
        """
        logger.info("Creating mock GDSC data...")
        np.random.seed(42)
        
        # Mock cell lines (mix of real and fictional)
        cell_lines = [
            'U87MG', 'U251MG', 'T98G', 'LN229', 'A172', 'SF268', 'SF295', 'SF539',
            'SNB19', 'SNB75', 'GBM_CELL_LINE_001', 'GBM_CELL_LINE_002', 'GBM_CELL_LINE_003',
            'GBM_CELL_LINE_004', 'GBM_CELL_LINE_005', 'GBM_CELL_LINE_006', 'GBM_CELL_LINE_007',
            'GBM_CELL_LINE_008', 'GBM_CELL_LINE_009', 'GBM_CELL_LINE_010'
        ]
        
        # Use FDA-approved drugs
        drugs = self.fda_approved_drugs[:20]  # Use first 20 for mock data
        
        # Generate mock sensitivity data
        data = []
        
        for cell_line in cell_lines:
            for drug in drugs:
                # Generate realistic IC50 values (log-normal distribution)
                ic50 = np.random.lognormal(mean=1, sigma=1)
                
                # Generate AUC values (0-1, where lower is more sensitive)
                auc = np.random.beta(2, 3)  # Skewed towards lower values
                
                # Generate max concentration tested
                max_conc = np.random.choice([1, 5, 10, 20, 50, 100])
                
                # Generate RMSE (measure of fit quality)
                rmse = np.random.uniform(0.1, 0.5)
                
                # Generate tissue type (mostly CNS for GBM)
                tissue = np.random.choice(['central_nervous_system', 'brain'], p=[0.8, 0.2])
                
                data.append({
                    'Cell line name': cell_line,
                    'Drug name': drug,
                    'GDSC Tissue descriptor 1': tissue,
                    'IC50 (µM)': ic50,
                    'AUC': auc,
                    'Max conc tested (µM)': max_conc,
                    'RMSE': rmse,
                    'Experiment ID': f'EXP_{len(data):06d}',
                    'Drug ID': f'DRUG_{drugs.index(drug):03d}',
                    'Cell line ID': f'CELL_{cell_lines.index(cell_line):03d}'
                })
        
        sensitivity_df = pd.DataFrame(data)
        
        # Create additional mock datasets
        mock_data = {
            'brain_sensitivity': sensitivity_df,
            'cell_lines': np.array(cell_lines),
            'drugs': np.array(drugs),
            'drug_info': self._create_mock_drug_info(drugs),
            'cell_line_info': self._create_mock_cell_line_info(cell_lines),
            'expression_data': self._create_mock_expression_data(cell_lines)
        }
        
        # Save mock data
        self._save_processed_data(mock_data)
        
        return mock_data
    
    def _create_mock_drug_info(self, drugs):
        """
        Create mock drug information.
        Pseudocode:
        1. For each drug, assign category, target, and mechanism.
        2. Return DataFrame.
        """
        drug_info = []
        
        for i, drug in enumerate(drugs):
            # Drug categories
            categories = ['Alkylating Agent', 'Anti-angiogenic', 'Tyrosine Kinase Inhibitor', 
                         'PARP Inhibitor', 'Immune Checkpoint Inhibitor', 'Cytotoxic Agent']
            
            drug_info.append({
                'Drug name': drug,
                'Drug ID': f'DRUG_{i:03d}',
                'Category': np.random.choice(categories),
                'Target': f'TARGET_{i:03d}',
                'Mechanism': f'Mechanism for {drug}',
                'FDA_approved': True,
                'Year_approved': np.random.randint(1990, 2023),
                'Molecular_weight': np.random.uniform(200, 800),
                'LogP': np.random.uniform(-2, 5)
            })
        
        return pd.DataFrame(drug_info)
    
    def _create_mock_cell_line_info(self, cell_lines):
        """
        Create mock cell line information.
        Pseudocode:
        1. For each cell line, assign tissue, origin, and other metadata.
        2. Return DataFrame.
        """
        cell_line_info = []
        
        for i, cell_line in enumerate(cell_lines):
            cell_line_info.append({
                'Cell line name': cell_line,
                'Cell line ID': f'CELL_{i:03d}',
                'Tissue': 'central_nervous_system',
                'Cancer type': 'Glioblastoma',
                'Source': np.random.choice(['ATCC', 'DSMZ', 'ECACC']),
                'Gender': np.random.choice(['Male', 'Female', 'Unknown']),
                'Age': np.random.randint(20, 80),
                'Mutation_count': np.random.randint(50, 500),
                'Doubling_time_hours': np.random.uniform(20, 48)
            })
        
        return pd.DataFrame(cell_line_info)
    
    def _create_mock_expression_data(self, cell_lines):
        """
        Create mock expression data for cell lines.
        Pseudocode:
        1. Generate synthetic gene expression matrix for cell lines.
        2. Return DataFrame.
        """
        n_genes = 1000
        gene_names = [f'GENE_{i:04d}' for i in range(n_genes)]
        
        # Generate expression matrix
        expression_values = np.random.lognormal(mean=5, sigma=2, size=(n_genes, len(cell_lines)))
        
        # Add some drug resistance genes
        resistance_genes = ['MGMT', 'EGFR', 'IDH1', 'TP53', 'PTEN', 'ABCB1', 'ABCG2']
        for gene in resistance_genes:
            if gene not in gene_names:
                gene_names.append(gene)
                # Add high expression for resistance genes
                resistance_expr = np.random.lognormal(mean=7, sigma=1, size=len(cell_lines))
                expression_values = np.vstack([expression_values, resistance_expr])
        
        expr_df = pd.DataFrame(
            expression_values,
            index=gene_names,
            columns=cell_lines
        )
        
        return expr_df
    
    def _save_processed_data(self, data):
        """
        Save processed or mock data to disk.
        Pseudocode:
        1. For each key in data, save as CSV or Excel as appropriate.
        """
        logger.info("Saving processed GDSC data...")
        
        # Save sensitivity data
        if 'brain_sensitivity' in data:
            sensitivity_file = self.data_dir / "gdsc_brain_sensitivity.csv"
            data['brain_sensitivity'].to_csv(sensitivity_file, index=False)
            logger.info(f"✓ Saved sensitivity data: {sensitivity_file}")
        
        # Save drug info
        if 'drug_info' in data:
            drug_info_file = self.data_dir / "gdsc_drug_info.csv"
            data['drug_info'].to_csv(drug_info_file, index=False)
            logger.info(f"✓ Saved drug info: {drug_info_file}")
        
        # Save cell line info
        if 'cell_line_info' in data:
            cell_line_file = self.data_dir / "gdsc_cell_line_info.csv"
            data['cell_line_info'].to_csv(cell_line_file, index=False)
            logger.info(f"✓ Saved cell line info: {cell_line_file}")
        
        # Save expression data
        if 'expression_data' in data:
            expression_file = self.data_dir / "gdsc_expression.csv"
            data['expression_data'].to_csv(expression_file)
            logger.info(f"✓ Saved expression data: {expression_file}")
    
    def load_processed_data(self):
        """
        Load all processed GDSC data from disk.
        Pseudocode:
        1. Load each processed data file if it exists.
        2. Return a dictionary of DataFrames.
        """
        logger.info("Loading processed GDSC data...")
        
        data = {}
        
        # Load sensitivity data
        sensitivity_file = self.data_dir / "gdsc_brain_sensitivity.csv"
        if sensitivity_file.exists():
            data['sensitivity'] = pd.read_csv(sensitivity_file)
            logger.info(f"✓ Loaded sensitivity data: {data['sensitivity'].shape}")
        
        # Load drug info
        drug_info_file = self.data_dir / "gdsc_drug_info.csv"
        if drug_info_file.exists():
            data['drug_info'] = pd.read_csv(drug_info_file)
            logger.info(f"✓ Loaded drug info: {data['drug_info'].shape}")
        
        # Load cell line info
        cell_line_file = self.data_dir / "gdsc_cell_line_info.csv"
        if cell_line_file.exists():
            data['cell_line_info'] = pd.read_csv(cell_line_file)
            logger.info(f"✓ Loaded cell line info: {data['cell_line_info'].shape}")
        
        # Load expression data
        expression_file = self.data_dir / "gdsc_expression.csv"
        if expression_file.exists():
            data['expression'] = pd.read_csv(expression_file, index_col=0)
            logger.info(f"✓ Loaded expression data: {data['expression'].shape}")
        
        return data
    
    def analyze_drug_response(self, sensitivity_data=None):
        """
        Analyze drug response data for key patterns.
        Pseudocode:
        1. If no data provided, load processed sensitivity data.
        2. Compute summary statistics and visualizations.
        3. Return analysis results.
        """
        if sensitivity_data is None:
            sensitivity_data = self.load_processed_data().get('sensitivity')
        
        if sensitivity_data is None:
            logger.warning("No sensitivity data available for analysis")
            return None
        
        logger.info("Analyzing drug response patterns...")
        
        # Calculate drug response statistics
        drug_stats = sensitivity_data.groupby('Drug name').agg({
            'IC50 (µM)': ['mean', 'std', 'median'],
            'AUC': ['mean', 'std', 'median'],
            'Cell line name': 'count'
        }).round(3)
        
        # Flatten column names
        drug_stats.columns = ['_'.join(col).strip() for col in drug_stats.columns]
        drug_stats = drug_stats.reset_index()
        
        # Identify most and least effective drugs
        most_effective = drug_stats.nsmallest(10, 'AUC_mean')
        least_effective = drug_stats.nlargest(10, 'AUC_mean')
        
        # Calculate resistance patterns
        resistance_analysis = {
            'drug_stats': drug_stats,
            'most_effective': most_effective,
            'least_effective': least_effective,
            'total_experiments': len(sensitivity_data),
            'unique_drugs': sensitivity_data['Drug name'].nunique(),
            'unique_cell_lines': sensitivity_data['Cell line name'].nunique()
        }
        
        logger.info(f"✓ Analyzed {resistance_analysis['total_experiments']} drug response experiments")
        logger.info(f"✓ Found {resistance_analysis['unique_drugs']} unique drugs")
        logger.info(f"✓ Found {resistance_analysis['unique_cell_lines']} unique cell lines")
        
        return resistance_analysis

# Usage example
if __name__ == "__main__":
    gdsc_parser = GDSCDataParser()
    
    # Download and parse data
    gdsc_data = gdsc_parser.download_and_parse_all()
    
    # Load processed data
    processed_data = gdsc_parser.load_processed_data()
    
    # Analyze drug response
    analysis = gdsc_parser.analyze_drug_response()
    
    if analysis:
        print("\nTop 5 most effective drugs:")
        print(analysis['most_effective'][['Drug name', 'AUC_mean', 'IC50 (µM)_mean']].head())
        
        print("\nTop 5 least effective drugs:")
        print(analysis['least_effective'][['Drug name', 'AUC_mean', 'IC50 (µM)_mean']].head()) 