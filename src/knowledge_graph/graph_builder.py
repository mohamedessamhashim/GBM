# src/knowledge_graph/graph_builder.py
import neo4j
from neo4j import GraphDatabase
import pandas as pd
import numpy as np
from pathlib import Path
import json
import requests
from typing import Dict, List, Tuple, Any
import logging
from tqdm import tqdm
import time

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

"""
Knowledge Graph Builder Module
-----------------------------
Handles construction and querying of the Neo4j-based knowledge graph for GBM drug discovery.
Integrates genes, patients, drugs, pathways, and relationships.
"""

class GBMKnowledgeGraph:
    def __init__(self, uri="bolt://localhost:7687", user="neo4j", password="gbm_project_2024"):
        """
        Initialize connection to Neo4j database.
        Pseudocode:
        1. Set up driver and session.
        2. Create uniqueness constraints for key node types.
        """
        self.driver = GraphDatabase.driver(uri, auth=(user, password))
        self.session = None
        
        # Data sources
        self.data_sources = {
            'string_db': 'https://stringdb-downloads.org/download/protein.links.v11.5/9606.protein.links.v11.5.txt.gz',
            'drugbank': 'https://go.drugbank.com/releases/latest',
            'reactome': 'https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt',
            'oncokb': 'https://www.oncokb.org/api/v1/',
            'cosmic': 'https://cancer.sanger.ac.uk/cosmic'
        }
        
        # Initialize session
        self._init_session()
        
    def _init_session(self):
        """
        Initialize database session and create constraints.
        Pseudocode:
        1. Start session.
        2. Create uniqueness constraints for nodes.
        """
        try:
            self.session = self.driver.session()
            
            # Create uniqueness constraints
            constraints = [
                "CREATE CONSTRAINT gene_name_unique IF NOT EXISTS FOR (g:Gene) REQUIRE g.name IS UNIQUE",
                "CREATE CONSTRAINT drug_name_unique IF NOT EXISTS FOR (d:Drug) REQUIRE d.name IS UNIQUE", 
                "CREATE CONSTRAINT pathway_id_unique IF NOT EXISTS FOR (p:Pathway) REQUIRE p.id IS UNIQUE",
                "CREATE CONSTRAINT patient_id_unique IF NOT EXISTS FOR (pt:Patient) REQUIRE pt.id IS UNIQUE"
            ]
            
            for constraint in constraints:
                try:
                    self.session.run(constraint)
                    logger.info(f"✓ Created constraint: {constraint.split('(')[1].split(')')[0]}")
                except Exception as e:
                    logger.warning(f"Constraint might already exist: {e}")
                    
        except Exception as e:
            logger.error(f"Failed to initialize Neo4j session: {e}")
            raise
    
    def clear_database(self):
        """
        Clear all nodes and relationships - USE WITH CAUTION.
        Pseudocode:
        1. Run Cypher query to delete all nodes and relationships.
        """
        logger.warning("Clearing entire database...")
        self.session.run("MATCH (n) DETACH DELETE n")
        logger.info("✓ Database cleared")
    
    def load_genes_from_expression_data(self, expr_df: pd.DataFrame, 
                                      clinical_df: pd.DataFrame):
        """
        Load genes and patients from expression data.
        Pseudocode:
        1. For each gene, create a node with statistics.
        2. For each patient, create a node with clinical info.
        3. Create HAS_EXPRESSION relationships for highly expressed genes.
        """
        logger.info("Loading genes and patients into knowledge graph...")
        
        # Load genes
        genes_loaded = 0
        for gene_name in tqdm(expr_df.index, desc="Loading genes"):
            gene_expression = expr_df.loc[gene_name]
            
            # Calculate gene statistics
            mean_expr = float(gene_expression.mean())
            std_expr = float(gene_expression.std())
            max_expr = float(gene_expression.max())
            min_expr = float(gene_expression.min())
            
            # Create gene node
            gene_query = """
            MERGE (g:Gene {name: $gene_name})
            SET g.mean_expression = $mean_expr,
                g.std_expression = $std_expr,
                g.max_expression = $max_expr,
                g.min_expression = $min_expr,
                g.sample_count = $sample_count
            """
            
            self.session.run(gene_query, {
                'gene_name': gene_name,
                'mean_expr': mean_expr,
                'std_expr': std_expr,
                'max_expr': max_expr,
                'min_expr': min_expr,
                'sample_count': len(gene_expression)
            })
            genes_loaded += 1
        
        logger.info(f"✓ Loaded {genes_loaded} genes")
        
        # Load patients
        patients_loaded = 0
        for _, patient in tqdm(clinical_df.iterrows(), desc="Loading patients", total=len(clinical_df)):
            patient_query = """
            MERGE (p:Patient {id: $patient_id})
            SET p.age_at_diagnosis = $age,
                p.gender = $gender,
                p.vital_status = $vital_status,
                p.overall_survival_days = $survival_days,
                p.mgmt_methylation_status = $mgmt_status,
                p.idh_mutation_status = $idh_status,
                p.treatment_response = $treatment_response
            """
            
            self.session.run(patient_query, {
                'patient_id': patient['patient_id'],
                'age': int(patient['age_at_diagnosis']) if pd.notna(patient['age_at_diagnosis']) else None,
                'gender': patient['gender'] if pd.notna(patient['gender']) else None,
                'vital_status': patient['vital_status'] if pd.notna(patient['vital_status']) else None,
                'survival_days': float(patient['overall_survival_days']) if pd.notna(patient['overall_survival_days']) else None,
                'mgmt_status': patient['mgmt_methylation_status'] if pd.notna(patient['mgmt_methylation_status']) else None,
                'idh_status': patient['idh_mutation_status'] if pd.notna(patient['idh_mutation_status']) else None,
                'treatment_response': patient['treatment_response'] if pd.notna(patient['treatment_response']) else None
            })
            patients_loaded += 1
        
        logger.info(f"✓ Loaded {patients_loaded} patients")
        
        # Create expression relationships
        self._create_expression_relationships(expr_df, clinical_df)
    
    def _create_expression_relationships(self, expr_df: pd.DataFrame, clinical_df: pd.DataFrame):
        """
        Create HAS_EXPRESSION relationships between patients and genes.
        Pseudocode:
        1. For each patient, find highly expressed genes.
        2. Create relationships with expression level and z-score.
        """
        logger.info("Creating patient-gene expression relationships...")
        
        # Get common samples
        common_samples = set(clinical_df['patient_id']) & set(expr_df.columns)
        
        relationships_created = 0
        for patient_id in tqdm(common_samples, desc="Creating expression relationships"):
            patient_expressions = expr_df[patient_id]
            
            # Only create relationships for highly expressed genes (optimization)
            high_expr_genes = patient_expressions[patient_expressions > patient_expressions.quantile(0.75)]
            
            for gene_name, expression_level in high_expr_genes.items():
                expr_query = """
                MATCH (p:Patient {id: $patient_id})
                MATCH (g:Gene {name: $gene_name})
                MERGE (p)-[e:HAS_EXPRESSION]->(g)
                SET e.expression_level = $expr_level,
                    e.z_score = $z_score
                """
                
                # Calculate z-score
                gene_mean = expr_df.loc[gene_name].mean()
                gene_std = expr_df.loc[gene_name].std()
                z_score = (expression_level - gene_mean) / gene_std if gene_std > 0 else 0
                
                self.session.run(expr_query, {
                    'patient_id': patient_id,
                    'gene_name': gene_name,
                    'expr_level': float(expression_level),
                    'z_score': float(z_score)
                })
                relationships_created += 1
        
        logger.info(f"✓ Created {relationships_created} expression relationships")
    
    def load_drugs_and_targets(self, gdsc_data: Dict):
        """
        Load drugs and their targets from GDSC data.
        Pseudocode:
        1. For each drug, create a node with properties.
        2. Create TARGETS relationships to genes.
        """
        logger.info("Loading drugs and drug-target relationships...")
        
        if 'drug_info' not in gdsc_data:
            logger.warning("No drug info available, creating mock drug data")
            self._create_mock_drug_data()
            return
        
        drug_info = gdsc_data['drug_info']
        drugs_loaded = 0
        
        for _, drug in tqdm(drug_info.iterrows(), desc="Loading drugs", total=len(drug_info)):
            drug_query = """
            MERGE (d:Drug {name: $drug_name})
            SET d.category = $category,
                d.target = $target,
                d.mechanism = $mechanism,
                d.fda_approved = $fda_approved,
                d.year_approved = $year_approved,
                d.molecular_weight = $mol_weight,
                d.logp = $logp
            """
            
            self.session.run(drug_query, {
                'drug_name': drug['Drug name'],
                'category': drug.get('Category', 'Unknown'),
                'target': drug.get('Target', 'Unknown'),
                'mechanism': drug.get('Mechanism', 'Unknown'),
                'fda_approved': drug.get('FDA_approved', True),
                'year_approved': int(drug.get('Year_approved', 2000)),
                'mol_weight': float(drug.get('Molecular_weight', 500)),
                'logp': float(drug.get('LogP', 2.0))
            })
            drugs_loaded += 1
        
        logger.info(f"✓ Loaded {drugs_loaded} drugs")
        
        # Create drug-target relationships
        self._create_drug_target_relationships(gdsc_data)
    
    def _create_mock_drug_data(self):
        """
        Create mock drug data for development.
        Pseudocode:
        1. Create nodes for common GBM drugs with properties.
        """
        logger.info("Creating mock drug data...")
        
        # Common GBM drugs
        drugs = [
            {'name': 'Temozolomide', 'category': 'Alkylating Agent', 'target': 'DNA'},
            {'name': 'Bevacizumab', 'category': 'Anti-angiogenic', 'target': 'VEGF'},
            {'name': 'Lomustine', 'category': 'Alkylating Agent', 'target': 'DNA'},
            {'name': 'Carmustine', 'category': 'Alkylating Agent', 'target': 'DNA'},
            {'name': 'Erlotinib', 'category': 'Tyrosine Kinase Inhibitor', 'target': 'EGFR'},
            {'name': 'Gefitinib', 'category': 'Tyrosine Kinase Inhibitor', 'target': 'EGFR'},
            {'name': 'Sorafenib', 'category': 'Tyrosine Kinase Inhibitor', 'target': 'VEGFR'},
            {'name': 'Olaparib', 'category': 'PARP Inhibitor', 'target': 'PARP'},
            {'name': 'Pembrolizumab', 'category': 'Immune Checkpoint Inhibitor', 'target': 'PD-1'},
            {'name': 'Nivolumab', 'category': 'Immune Checkpoint Inhibitor', 'target': 'PD-1'}
        ]
        
        for drug in drugs:
            drug_query = """
            MERGE (d:Drug {name: $drug_name})
            SET d.category = $category,
                d.target = $target,
                d.fda_approved = true,
                d.year_approved = 2010
            """
            
            self.session.run(drug_query, {
                'drug_name': drug['name'],
                'category': drug['category'],
                'target': drug['target']
            })
        
        logger.info(f"✓ Created {len(drugs)} mock drugs")
    
    def _create_drug_target_relationships(self, gdsc_data: Dict):
        """
        Create TARGETS relationships between drugs and genes.
        Pseudocode:
        1. For each drug, find target gene and create relationship.
        """
        logger.info("Creating drug-target relationships...")
        
        if 'drug_info' not in gdsc_data:
            logger.warning("No drug info available for target relationships")
            return
        
        drug_info = gdsc_data['drug_info']
        relationships_created = 0
        
        for _, drug in drug_info.iterrows():
            target = drug.get('Target', 'Unknown')
            if target != 'Unknown':
                # Try to find matching gene
                target_query = """
                MATCH (d:Drug {name: $drug_name})
                MATCH (g:Gene) WHERE g.name CONTAINS $target OR g.name = $target
                MERGE (d)-[t:TARGETS]->(g)
                SET t.mechanism = $mechanism
                """
                
                try:
                    result = self.session.run(target_query, {
                        'drug_name': drug['Drug name'],
                        'target': target,
                        'mechanism': drug.get('Mechanism', 'Unknown')
                    })
                    
                    if result.single() is not None:
                        relationships_created += 1
                        
                except Exception as e:
                    # Skip if target gene not found
                    continue
        
        logger.info(f"✓ Created {relationships_created} drug-target relationships")
    
    def load_drug_response_data(self, gdsc_data: Dict):
        """
        Load drug response data from GDSC.
        Pseudocode:
        1. For each response, create RESPONDS_TO relationship between cell line and drug.
        """
        logger.info("Loading drug response data...")
        
        if 'sensitivity' not in gdsc_data:
            logger.warning("No sensitivity data available")
            return
        
        sensitivity_data = gdsc_data['sensitivity']
        responses_loaded = 0
        
        for _, response in tqdm(sensitivity_data.iterrows(), desc="Loading drug responses", total=len(sensitivity_data)):
            response_query = """
            MATCH (d:Drug {name: $drug_name})
            MERGE (c:CellLine {name: $cell_line})
            MERGE (c)-[r:RESPONDS_TO]->(d)
            SET r.ic50 = $ic50,
                r.auc = $auc,
                r.max_conc = $max_conc,
                r.rmse = $rmse,
                r.sensitive = $sensitive
            """
            
            # Determine if cell line is sensitive (AUC < 0.5)
            sensitive = response['AUC'] < 0.5
            
            self.session.run(response_query, {
                'drug_name': response['Drug name'],
                'cell_line': response['Cell line name'],
                'ic50': float(response['IC50 (µM)']),
                'auc': float(response['AUC']),
                'max_conc': float(response['Max conc tested (µM)']),
                'rmse': float(response['RMSE']),
                'sensitive': sensitive
            })
            responses_loaded += 1
        
        logger.info(f"✓ Loaded {responses_loaded} drug response relationships")
    
    def create_protein_interactions(self):
        """
        Create protein-protein interactions (mock for development).
        Pseudocode:
        1. Select genes from database.
        2. Randomly create INTERACTS_WITH relationships.
        """
        logger.info("Creating protein-protein interactions...")
        
        # Get some genes from the database
        gene_query = "MATCH (g:Gene) RETURN g.name as name LIMIT 100"
        result = self.session.run(gene_query)
        genes = [record['name'] for record in result]
        
        if len(genes) < 10:
            logger.warning("Not enough genes to create interactions")
            return
        
        # Create random interactions
        interactions_created = 0
        np.random.seed(42)
        
        for i in range(min(500, len(genes) * 2)):  # Create up to 500 interactions
            gene1, gene2 = np.random.choice(genes, 2, replace=False)
            confidence = np.random.randint(400, 1000)
            
            interaction_query = """
            MATCH (g1:Gene {name: $gene1})
            MATCH (g2:Gene {name: $gene2})
            MERGE (g1)-[i:INTERACTS_WITH]->(g2)
            SET i.confidence_score = $score,
                i.source = 'MOCK'
            """
            
            self.session.run(interaction_query, {
                'gene1': gene1,
                'gene2': gene2,
                'score': confidence
            })
            interactions_created += 1
        
        logger.info(f"✓ Created {interactions_created} protein interactions")
    
    def create_pathway_relationships(self):
        """
        Create pathway relationships (mock for development).
        Pseudocode:
        1. For each pathway, create node and BELONGS_TO relationships to genes.
        """
        logger.info("Creating pathway relationships...")
        
        # Common GBM pathways
        pathways = [
            {'name': 'PI3K/AKT/mTOR', 'genes': ['PIK3CA', 'PIK3R1', 'AKT1', 'MTOR']},
            {'name': 'RAS/RAF/MEK/ERK', 'genes': ['KRAS', 'BRAF', 'MAP2K1', 'MAPK1']},
            {'name': 'p53', 'genes': ['TP53', 'MDM2', 'MDM4']},
            {'name': 'RB', 'genes': ['RB1', 'CDKN2A', 'CDK4', 'CDK6']},
            {'name': 'EGFR', 'genes': ['EGFR', 'ERBB2', 'ERBB3', 'ERBB4']},
            {'name': 'DNA Repair', 'genes': ['MGMT', 'PARP1', 'BRCA1', 'BRCA2']}
        ]
        
        pathways_created = 0
        relationships_created = 0
        
        for pathway in pathways:
            # Create pathway node
            pathway_query = """
            MERGE (p:Pathway {name: $pathway_name})
            SET p.category = 'Signaling Pathway'
            """
            
            self.session.run(pathway_query, {
                'pathway_name': pathway['name']
            })
            pathways_created += 1
            
            # Create pathway-gene relationships
            for gene_name in pathway['genes']:
                gene_pathway_query = """
                MATCH (p:Pathway {name: $pathway_name})
                MATCH (g:Gene) WHERE g.name CONTAINS $gene_name OR g.name = $gene_name
                MERGE (g)-[r:BELONGS_TO]->(p)
                """
                
                try:
                    result = self.session.run(gene_pathway_query, {
                        'pathway_name': pathway['name'],
                        'gene_name': gene_name
                    })
                    
                    if result.single() is not None:
                        relationships_created += 1
                        
                except Exception as e:
                    # Skip if gene not found
                    continue
        
        logger.info(f"✓ Created {pathways_created} pathways")
        logger.info(f"✓ Created {relationships_created} pathway-gene relationships")
    
    def run_graph_queries(self):
        """
        Run example graph queries to test the database.
        Pseudocode:
        1. Run a set of predefined Cypher queries and log results.
        """
        logger.info("Running example graph queries...")
        
        queries = [
            {
                'name': 'Count all nodes by type',
                'query': 'MATCH (n) RETURN labels(n) as type, count(n) as count ORDER BY count DESC'
            },
            {
                'name': 'Count all relationships by type',
                'query': 'MATCH ()-[r]->() RETURN type(r) as type, count(r) as count ORDER BY count DESC'
            },
            {
                'name': 'Find drugs targeting EGFR',
                'query': '''
                MATCH (d:Drug)-[:TARGETS]->(g:Gene)
                WHERE g.name CONTAINS 'EGFR'
                RETURN d.name as drug, g.name as target
                '''
            },
            {
                'name': 'Find patients with high MGMT expression',
                'query': '''
                MATCH (p:Patient)-[e:HAS_EXPRESSION]->(g:Gene)
                WHERE g.name CONTAINS 'MGMT' AND e.z_score > 1
                RETURN p.id as patient, g.name as gene, e.z_score as z_score
                ORDER BY e.z_score DESC
                LIMIT 10
                '''
            },
            {
                'name': 'Find drug resistance patterns',
                'query': '''
                MATCH (c:CellLine)-[r:RESPONDS_TO]->(d:Drug)
                WHERE r.sensitive = false
                RETURN d.name as drug, count(r) as resistant_cell_lines
                ORDER BY resistant_cell_lines DESC
                LIMIT 10
                '''
            }
        ]
        
        for query_info in queries:
            try:
                result = self.session.run(query_info['query'])
                print(f"\n{query_info['name']}:")
                for record in result:
                    print(f"  {dict(record)}")
            except Exception as e:
                print(f"Query failed: {e}")
    
    def close(self):
        """
        Close the Neo4j session and driver.
        """
        if self.session:
            self.session.close()
        if self.driver:
            self.driver.close()
        logger.info("✓ Database connection closed")

# Usage example
if __name__ == "__main__":
    # Initialize knowledge graph
    kg = GBMKnowledgeGraph()
    
    try:
        # Clear database (optional)
        # kg.clear_database()
        
        # Load data (you would load your actual data here)
        # kg.load_genes_from_expression_data(expr_df, clinical_df)
        # kg.load_drugs_and_targets(gdsc_data)
        # kg.load_drug_response_data(gdsc_data)
        
        # Create additional relationships
        kg.create_protein_interactions()
        kg.create_pathway_relationships()
        
        # Run example queries
        kg.run_graph_queries()
        
    finally:
        kg.close() 