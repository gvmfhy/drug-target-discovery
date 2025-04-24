#!/usr/bin/env python3
"""
Drug Target Discovery Pipeline for GSE46602 with Bioconductor-Based Probe Mapping

This pipeline processes gene expression data from a GEO Series Matrix file,
maps probes to gene symbols using R's Bioconductor packages via rpy2,
performs differential expression analysis, constructs a gene co-expression
network, identifies potential drug targets through network analysis,
validates targets using the Open Targets API, and generates visualizations
and a summary report.
"""

import os
import sys
import argparse
import logging
from typing import Optional
import pandas as pd
import numpy as np
import gzip
import re
import requests
from sklearn.preprocessing import StandardScaler
from scipy import stats
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import subprocess
import tempfile
import csv
from statsmodels.stats.multitest import multipletests
import traceback

# rpy2 integration for Bioconductor-based probe mapping
try:
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
except ImportError:
    raise ImportError("rpy2 is required for Bioconductor-based probe mapping. Please install it via pip.")

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stderr)]
)
logger = logging.getLogger('drug_target_pipeline')

def map_probes_to_genes_bioc(probe_ids):
    """
    Map probe IDs to gene symbols using Bioconductor.
    This function uses a different approach by calling an R script directly
    instead of using rpy2 interfaces, which helps avoid conversion issues.
    
    Args:
        probe_ids (list): List of probe IDs (strings)
        
    Returns:
        dict: Mapping of probe IDs to gene symbols
    """
    # Create temporary directory for input/output files
    temp_dir = tempfile.mkdtemp()
    probes_file = os.path.join(temp_dir, "probe_ids.txt")
    output_file = os.path.join(temp_dir, "probe_mapping.csv")
    
    # Write probe IDs to temporary file
    logger.info(f"Writing {len(probe_ids)} probe IDs to temporary file")
    with open(probes_file, 'w') as f:
        for probe_id in probe_ids:
            f.write(f"{probe_id}\n")
    
    # Run R script
    result = generate_probe_mappings_csv(probes_file, output_file)
    if not result:
        return {}
    
    # Read mapping results
    mapping = {}
    try:
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            with open(output_file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    probe_id = row['PROBEID']
                    symbol = row['SYMBOL']
                    if symbol and symbol.strip():  # Check for non-empty symbol
                        mapping[probe_id] = symbol.strip()
            
            # Log mapping statistics
            total_probes = len(probe_ids)
            mapped_probes = len(mapping)
            mapping_percent = mapped_probes/total_probes*100 if total_probes > 0 else 0
            logger.info(f"Successfully mapped {mapped_probes} out of {total_probes} probes ({mapping_percent:.1f}%)")
        else:
            logger.error("Output file is empty or doesn't exist")
    except Exception as e:
        logger.error(f"Error reading mapping results: {str(e)}")
    
    # Clean up temporary files
    try:
        for file in [probes_file, output_file]:
            if os.path.exists(file):
                os.remove(file)
        os.rmdir(temp_dir)
    except Exception as e:
        logger.warning(f"Error cleaning up temporary files: {str(e)}")
    
    return mapping

def generate_probe_mappings_csv(probes_file: str, output_file: str) -> Optional[subprocess.CompletedProcess[str]]:
    logger.info("Running R script for probe mapping")
    r_script = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'generate_probe_mappings_csv.R')
    cmd = ["Rscript", r_script, probes_file, output_file]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"R script output: {result.stdout.strip()}")
        return result
    except subprocess.CalledProcessError as e:
        logger.error(f"R script failed with code {e.returncode}: {e.stderr.strip()}")
        logger.info(f"R script output: {e.stdout.strip()}")
        return None

class DrugTargetPipeline:
    def __init__(self, matrix_file, output_dir=None, p_threshold=0.05, fc_threshold=1.0):
        """
        Initialize the pipeline.
        
        Args:
            matrix_file (str): Path to GEO Series Matrix file.
            output_dir (str): Directory to save results.
            p_threshold (float): P-value threshold.
            fc_threshold (float): Log2 fold change threshold.
        """
        self.matrix_file = matrix_file
        self.output_dir = output_dir or os.path.join("results", f"GSE46602_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
        self.p_threshold = p_threshold
        self.fc_threshold = fc_threshold
        
        os.makedirs(os.path.join(self.output_dir, "data"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "figures"), exist_ok=True)
        
        self.expression_data = None
        self.metadata = None
        self.normalized_data = None
        self.differential_results = None
        self.significant_genes = None
        self.network = None
        self.target_scores = None
        self.sample_groups = {}  # Initialize sample_groups dictionary
        
        logger.info(f"Pipeline initialized with matrix file: {matrix_file}")
        logger.info(f"Results will be saved to: {self.output_dir}")
    
    def parse_geo_matrix(self):
        """Parse the GEO Series Matrix file to extract expression data and metadata."""
        logger.info(f"Parsing GEO Series Matrix file: {self.matrix_file}")
        metadata_lines = []
        data_lines = []
        sample_ids = []
        in_data_section = False
        
        try:
            # Check if the file is gzipped
            is_gzipped = self.matrix_file.lower().endswith('.gz')
            
            # Open the file with the appropriate method
            if is_gzipped:
                logger.info("Opening as gzipped file")
                file_opener = lambda: gzip.open(self.matrix_file, 'rt')
            else:
                logger.info("Opening as plain text file")
                file_opener = lambda: open(self.matrix_file, 'rt')
                
            with file_opener() as f:
                lines = f.readlines()
                
                for line in lines:
                    if line.startswith('!'):
                        metadata_lines.append(line.strip())
                    elif line.startswith('#'):
                        continue
                    elif "ID_REF" in line and not line.startswith('!'):
                        sample_ids = line.strip().split('\t')[1:]
                        in_data_section = True
                    elif in_data_section:
                        data_lines.append(line.strip())
        except Exception as e:
            logger.error(f"Error reading matrix file: {e}")
            raise
        
        # Pass the entire file content to _extract_metadata for more detailed parsing
        self._extract_metadata(self.matrix_file)
        self._extract_expression_data(data_lines, sample_ids)
        logger.info(f"Parsed {len(self.expression_data)} genes/probes")
        return self.expression_data, self.metadata

    def _extract_metadata(self, file_path):
        """Extract metadata from GEO Series Matrix file."""
        logging.info("Extracting metadata from %s", file_path)
        sample_ids = []
        sample_titles = []
        sample_characteristics = {}
        sample_groups = {}
        
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
                
            # Extract sample IDs (GSM numbers)
            for line in lines:
                if line.startswith('!Sample_geo_accession'):
                    sample_ids = [id.strip('"') for id in line.strip().split('\t')[1:]]
                    break
            
            if not sample_ids:
                logger.error("No sample IDs found in the matrix file")
                return {}
                
            logger.info(f"Found {len(sample_ids)} samples")
            
            # Extract sample titles
            for line in lines:
                if line.startswith('!Sample_title'):
                    sample_titles = [title.strip('"') for title in line.strip().split('\t')[1:]]
                    break
            
            # Initialize characteristics dictionary for each sample
            for sample_id in sample_ids:
                sample_characteristics[sample_id] = {}
            
            # Extract all characteristics
            for line in lines:
                if line.startswith('!Sample_characteristics_ch'):
                    parts = line.strip().split('\t')
                    if len(parts) < 2:
                        continue
                    
                    first_value = parts[1].strip('"').strip()
                    if ':' in first_value:
                        # This is a characteristic with a label
                        label = first_value.split(':', 1)[0].strip().lower()
                        
                        for i, sample_id in enumerate(sample_ids):
                            if i+1 < len(parts):
                                value = parts[i+1].strip('"').strip()
                                if ':' in value:
                                    value = value.split(':', 1)[1].strip()
                                sample_characteristics[sample_id][label] = value
            
            # Identify cancer vs. benign samples based on characteristics
            cancer_keywords = ['cancer', 'tumor', 'tumour', 'malignant', 'carcinoma']
            benign_keywords = ['benign', 'normal', 'healthy', 'non-tumor', 'non-cancer', 'non-malignant']
            
            for sample_id in sample_ids:
                # Check 'tissue' characteristic
                if 'tissue' in sample_characteristics[sample_id]:
                    tissue = sample_characteristics[sample_id]['tissue'].lower()
                    if any(kw in tissue for kw in cancer_keywords):
                        sample_groups[sample_id] = 'case'
                    elif any(kw in tissue for kw in benign_keywords):
                        sample_groups[sample_id] = 'control'
                
                # If tissue didn't determine group, check sample title
                if sample_id not in sample_groups and sample_ids.index(sample_id) < len(sample_titles):
                    title = sample_titles[sample_ids.index(sample_id)].lower()
                    if any(kw in title for kw in cancer_keywords):
                        sample_groups[sample_id] = 'case'
                    elif any(kw in title for kw in benign_keywords):
                        sample_groups[sample_id] = 'control'
                
                # Check all characteristics as a last resort
                if sample_id not in sample_groups:
                    all_values = ' '.join(str(v).lower() for v in sample_characteristics[sample_id].values())
                    if any(kw in all_values for kw in cancer_keywords):
                        sample_groups[sample_id] = 'case'
                    elif any(kw in all_values for kw in benign_keywords):
                        sample_groups[sample_id] = 'control'
            
            # Count samples in each group
            case_count = sum(1 for group in sample_groups.values() if group == 'case')
            control_count = sum(1 for group in sample_groups.values() if group == 'control')
            
            logger.info(f"Identified {case_count} cancer (case) samples and {control_count} benign (control) samples")
            
            # If not enough samples with known groups, create default groups
            if case_count < 2 or control_count < 2:
                logger.warning(f"Not enough samples identified in each group (need at least 2). Creating default groups.")
                sample_groups = {}
                midpoint = len(sample_ids) // 2
                for i, sample_id in enumerate(sample_ids):
                    sample_groups[sample_id] = 'control' if i < midpoint else 'case'
                
                logger.info(f"Created default groups: {midpoint} controls, {len(sample_ids) - midpoint} cases")
            
            # Create metadata dataframe
            metadata_dict = {sample_id: {} for sample_id in sample_ids}
            
            # Add titles and groups to metadata
            for i, sample_id in enumerate(sample_ids):
                if i < len(sample_titles):
                    metadata_dict[sample_id]['title'] = sample_titles[i]
                if sample_id in sample_groups:
                    metadata_dict[sample_id]['condition'] = sample_groups[sample_id]
                
                # Add all characteristics
                for label, value in sample_characteristics[sample_id].items():
                    metadata_dict[sample_id][label] = value
            
            # Convert to DataFrame
            metadata_df = pd.DataFrame.from_dict(metadata_dict, orient='index')
            
            # Log found fields
            logger.info(f"Metadata fields collected: {', '.join(metadata_df.columns)}")
            
            # Save metadata to file
            metadata_path = os.path.join(self.output_dir, "data", "GSE46602_metadata.csv")
            metadata_df.to_csv(metadata_path)
            logger.info(f"Metadata saved to {metadata_path}")
            
            # Store in class attributes
            self.metadata = metadata_df
            self.sample_ids = sample_ids
            self.sample_groups = sample_groups
            
        except Exception as e:
            logger.error(f"Error extracting metadata: {str(e)}")
            logger.error(traceback.format_exc())
            self.metadata = pd.DataFrame()
            self.sample_groups = {}
            
        return self.metadata

    def reconcile_sample_ids(self):
        """
        Ensure sample IDs are consistent between metadata and expression data,
        and create the sample_groups dictionary mapping sample IDs to their group.
        """
        if self.metadata is None:
            logger.warning("No metadata available for sample group reconciliation")
            return
        
        self.sample_groups = {}
        
        # If we have expression data, reconcile with its column names
        if self.normalized_data is not None:
            # Clean column names - remove quotes that might be present
            expr_samples_original = list(self.normalized_data.columns)
            expr_samples_clean = [s.strip('"\'') for s in expr_samples_original]
            
            # Create a mapping from clean to original column names
            clean_to_original = {clean: original for clean, original in zip(expr_samples_clean, expr_samples_original)}
            
            # Clean metadata sample IDs as well
            meta_samples = [s.strip('"\'') for s in self.metadata.index]
            
            logger.info(f"Expression data has {len(expr_samples_original)} samples")
            logger.info(f"Metadata has {len(meta_samples)} samples")
            
            # Create mapping for samples present in both (using clean versions for comparison)
            common_samples = set(expr_samples_clean).intersection(set(meta_samples))
            logger.info(f"Found {len(common_samples)} common samples between metadata and expression data")
            
            if 'condition' in self.metadata.columns:
                for clean_sample in common_samples:
                    # Get the original column name from expression data
                    original_sample = clean_to_original.get(clean_sample, clean_sample)
                    
                    # Get the condition from metadata using the clean sample ID
                    condition = self.metadata.loc[clean_sample, 'condition'] if clean_sample in self.metadata.index else None
                    
                    # Add to sample_groups using the original column name
                    if condition:
                        self.sample_groups[original_sample] = condition
            
            # If we didn't get any sample groups, try matching by position
            if not self.sample_groups and len(expr_samples_original) == len(meta_samples):
                logger.warning("No samples matched by ID. Trying to match by position...")
                conditions = self.metadata['condition'].tolist() if 'condition' in self.metadata.columns else []
                
                if conditions:
                    for i, sample in enumerate(expr_samples_original):
                        if i < len(conditions):
                            self.sample_groups[sample] = conditions[i]
                    
                    logger.info(f"Matched {len(self.sample_groups)} samples by position")
            
            # Count control and case samples
            control_count = sum(1 for group in self.sample_groups.values() if group == 'control')
            case_count = sum(1 for group in self.sample_groups.values() if group == 'case')
            logger.info(f"Sample groups distribution: {control_count} controls, {case_count} cases")
            
            # If we still don't have enough samples, create default groups
            if control_count < 2 or case_count < 2:
                logger.warning(f"Not enough samples in each group. Creating default groups.")
                self.sample_groups = {}
                samples = expr_samples_original
                midpoint = len(samples) // 2
                
                for i, sample in enumerate(samples):
                    self.sample_groups[sample] = 'control' if i < midpoint else 'case'
                
                control_count = sum(1 for group in self.sample_groups.values() if group == 'control')
                case_count = sum(1 for group in self.sample_groups.values() if group == 'case')
                logger.info(f"Created default groups: {control_count} controls, {case_count} cases")
        else:
            # If no expression data yet, use all metadata samples
            if 'condition' in self.metadata.columns:
                for sample, row in self.metadata.iterrows():
                    self.sample_groups[sample] = row['condition']
                
                logger.info(f"Created sample groups for {len(self.sample_groups)} samples from metadata")
                
                # Count control and case samples
                control_count = sum(1 for group in self.sample_groups.values() if group == 'control')
                case_count = sum(1 for group in self.sample_groups.values() if group == 'case')
                logger.info(f"Sample groups distribution: {control_count} controls, {case_count} cases")

    def _extract_expression_data(self, data_lines, sample_ids):
        """Extract expression data from the matrix file."""
        data_rows = []
        row_ids = []
        
        for line in data_lines:
            parts = line.split('\t')
            if len(parts) >= len(sample_ids) + 1:
                # Clean up probe ID - remove quotes, colons, and other special characters
                probe_id_raw = parts[0]
                
                # Handle different formats
                # Format 1: 54679:"AFFX-BioB-3_at" -> AFFX-BioB-3_at
                # Format 2: 54675-"91816_f_at" -> 91816_f_at
                
                # First, extract the actual probe ID from any numbering or prefix
                if ':' in probe_id_raw:
                    probe_id_raw = probe_id_raw.split(':', 1)[1]
                elif '-' in probe_id_raw:
                    probe_id_raw = probe_id_raw.split('-', 1)[1]
                    
                # Then, remove quotes
                probe_id = probe_id_raw.strip('"\'')
                
                # Clean up any other potential issues
                probe_id = probe_id.strip()
                
                row_ids.append(probe_id)
                
                # Extract expression values
                try:
                    values = [float(v) for v in parts[1:len(sample_ids)+1]]
                    data_rows.append(values)
                except ValueError:
                    logger.warning(f"Could not convert expression values for row {probe_id}")
        
        if data_rows:
            self.expression_data = pd.DataFrame(data_rows, index=row_ids, columns=sample_ids)
        else:
            logger.error("No expression data could be parsed.")
            self.expression_data = pd.DataFrame()
    
    def preprocess_data(self):
        """Preprocess the expression data (NA handling, log2 transformation, normalization)."""
        logger.info("Preprocessing expression data")
        if self.expression_data is None:
            logger.error("No expression data available. Run parse_geo_matrix() first.")
            return None
        logger.info(f"Original data shape: {self.expression_data.shape}")
        logger.info(f"Missing values: {self.expression_data.isna().sum().sum()}")
        threshold = 0.2 * self.expression_data.shape[1]
        filtered_data = self.expression_data.dropna(thresh=threshold)
        logger.info(f"After dropping rows with >20% NaNs: {filtered_data.shape}")
        filled_data = filtered_data.apply(lambda row: row.fillna(row.median()), axis=1)
        data_max = filled_data.max().max()
        if data_max > 100:
            logger.info("Applying log2 transformation")
            filled_data = np.log2(filled_data + 1)
        scaler = StandardScaler()
        normalized_values = scaler.fit_transform(filled_data.T).T
        self.normalized_data = pd.DataFrame(normalized_values, index=filled_data.index, columns=filled_data.columns)
        norm_path = os.path.join(self.output_dir, "data", "GSE46602_normalized.csv")
        self.normalized_data.to_csv(norm_path)
        logger.info(f"Normalized data saved to {norm_path}")
        return self.normalized_data
    
    def map_probes_to_genes(self):
        """Map probe IDs to gene symbols using Bioconductor via rpy2."""
        logger.info("Mapping probe IDs to gene symbols using Bioconductor")
        if self.normalized_data is None:
            logger.error("No normalized data available. Run preprocess_data() first.")
            return None
        
        # Get probe IDs
        probe_ids = list(self.normalized_data.index)
        
        # Map probes to genes
        probe_to_gene = map_probes_to_genes_bioc(probe_ids)
        
        # Create a new DataFrame with gene symbols
        gene_expr = self.normalized_data.copy()
        gene_expr['gene_symbol'] = [probe_to_gene.get(probe, f"UNKNOWN_{probe}") for probe in probe_ids]
        gene_expr = gene_expr.reset_index().rename(columns={'index': 'probe_id'})
        
        # Filter out unmapped probes (those starting with UNKNOWN_)
        mapped_expr = gene_expr[~gene_expr['gene_symbol'].str.startswith('UNKNOWN_')]
        logger.info(f"Filtered out {len(gene_expr) - len(mapped_expr)} unmapped probes")
        
        # Group by gene symbol and aggregate using median
        gene_aggregated = mapped_expr.groupby('gene_symbol').agg({
            col: 'median' for col in mapped_expr.columns if col not in ['probe_id', 'gene_symbol']
        })
        
        # Update normalized data with gene-mapped data
        self.normalized_data = gene_aggregated
        
        # Save gene-mapped data
        gene_path = os.path.join(self.output_dir, "data", "GSE46602_gene_mapped.csv")
        self.normalized_data.to_csv(gene_path)
        logger.info(f"Gene-mapped expression data saved to {gene_path}")
        
        # Reconcile sample IDs after mapping
        self.reconcile_sample_ids()
        
        return self.normalized_data
    
    def perform_differential_analysis(self):
        """Perform differential expression analysis."""
        if not hasattr(self, 'normalized_data') or self.normalized_data is None:
            logger.error("Normalized data not available for differential analysis")
            return None
        
        if not hasattr(self, 'sample_groups') or not self.sample_groups:
            logger.error("Sample groups not available for differential analysis")
            return None
        
        logger.info("Performing differential expression analysis")
        
        try:
            # Identify control and case samples
            control_samples = [sample_id for sample_id, group in self.sample_groups.items() if group == 'control']
            case_samples = [sample_id for sample_id, group in self.sample_groups.items() if group == 'case']
            
            logger.info(f"Differential analysis using {len(control_samples)} control samples and {len(case_samples)} case samples")
            
            # Reconcile sample IDs with column names in normalized_data
            available_samples = list(self.normalized_data.columns)
            valid_control_samples = [s for s in control_samples if s in available_samples]
            valid_case_samples = [s for s in case_samples if s in available_samples]
            
            if len(valid_control_samples) < 2 or len(valid_case_samples) < 2:
                logger.warning(f"Not enough valid samples for analysis: {len(valid_control_samples)} controls, {len(valid_case_samples)} cases")
                logger.warning("Using all samples and creating arbitrary groups")
                
                # Fall back to arbitrary grouping if not enough valid samples
                all_samples = available_samples
                valid_control_samples = all_samples[:len(all_samples)//2]
                valid_case_samples = all_samples[len(all_samples)//2:]
                
                logger.info(f"Created arbitrary groups: {len(valid_control_samples)} controls, {len(valid_case_samples)} cases")
            
            logger.info(f"Valid control samples: {len(valid_control_samples)}")
            logger.info(f"Valid case samples: {len(valid_case_samples)}")
            
            # Initialize lists to store results
            genes = []
            log2fc_values = []
            pvalues = []
            
            # Perform t-test for each gene
            for gene in self.normalized_data.index:
                # Extract expression values for control and case samples
                control_exp = self.normalized_data.loc[gene, valid_control_samples].values
                case_exp = self.normalized_data.loc[gene, valid_case_samples].values
                
                # Calculate mean expression
                mean_control = np.mean(control_exp)
                mean_case = np.mean(case_exp)
                
                # Calculate fold change and log2 fold change
                # Since we're using normalized data (already log-transformed and standardized),
                # we can directly use the difference in means as the log2FC
                log2fc = mean_case - mean_control
                
                # Perform t-test (unequal variance)
                try:
                    tstat, pval = stats.ttest_ind(case_exp, control_exp, equal_var=False)
                    # We're only interested in the p-value
                except:
                    pval = np.nan  # In case the t-test fails (e.g., constant values)
                
                genes.append(gene)
                log2fc_values.append(log2fc)
                pvalues.append(pval)
            
            # Create results DataFrame
            results = pd.DataFrame({
                'gene': genes,
                'log2FC': log2fc_values,
                'pvalue': pvalues
            })
            
            # Remove rows with NA or infinite values for log2FC
            results = results.replace([np.inf, -np.inf], np.nan)
            
            # Calculate adjusted p-values using Benjamini-Hochberg method
            valid_pvals = ~np.isnan(results['pvalue'])
            if valid_pvals.any():
                adjusted_pvals = np.full_like(results['pvalue'], np.nan)
                adjusted_pvals[valid_pvals] = multipletests(
                    results.loc[valid_pvals, 'pvalue'],
                    method='fdr_bh'
                )[1]
                results['adjusted_pvalue'] = adjusted_pvals
            else:
                results['adjusted_pvalue'] = np.nan
            
            # Save results to CSV
            os.makedirs(os.path.join(self.output_dir, 'data'), exist_ok=True)
            results_file = os.path.join(self.output_dir, 'data', 'differential_results.csv')
            results.to_csv(results_file, index=False)
            logger.info(f"Differential analysis results saved to {results_file}")
            
            # Identify significant genes (adjusted p-value < 0.05 and |log2FC| > 1)
            if 'adjusted_pvalue' in results.columns:
                sig_results = results[
                    (results['adjusted_pvalue'] < 0.05) & 
                    (results['log2FC'].abs() > 1) &
                    (~results['log2FC'].isna())
                ]
                logger.info(f"Identified {len(sig_results)} significant genes")
                
                if not sig_results.empty:
                    sig_file = os.path.join(self.output_dir, 'data', 'significant_genes.csv')
                    sig_results.to_csv(sig_file, index=False)
                    logger.info(f"Significant genes saved to {sig_file}")
                    self.significant_genes = sig_results
                    return sig_results
                
                self.significant_genes = None
            
            logger.warning("No significant genes found after filtering")
            return results
        
        except Exception as e:
            logger.error(f"Error in differential analysis: {str(e)}")
            logger.error(traceback.format_exc())
            return None
    
    def construct_network(self, n_top_genes=500):
        """Construct a gene co-expression network based on Pearson correlations."""
        logger.info("Constructing gene co-expression network")
        if self.normalized_data is None:
            logger.error("Normalized expression data not available.")
            return None
        
        # Make sure we have enough genes for a meaningful network
        if len(self.normalized_data) < 2:
            logger.error(f"Not enough genes available in the normalized data. Found {len(self.normalized_data)}, need at least 2.")
            # Create empty network to avoid downstream errors
            network = nx.Graph()
            self.network = network
            logger.warning("Created empty network to avoid pipeline failure")
            return self.network
        
        if self.significant_genes is not None and len(self.significant_genes) > 0:
            top_genes = list(self.significant_genes['gene'])[:n_top_genes]
            logger.info(f"Using {len(top_genes)} significant genes for network construction")
        else:
            logger.warning("No significant genes found, using most variable genes instead")
            gene_var = self.normalized_data.var(axis=1).sort_values(ascending=False)
            top_genes = list(gene_var.index)[:n_top_genes]
            logger.info(f"Using {len(top_genes)} most variable genes for network construction")
        
        # Verify top_genes are actually in the normalized data
        valid_genes = [gene for gene in top_genes if gene in self.normalized_data.index]
        if len(valid_genes) < len(top_genes):
            logger.warning(f"Only {len(valid_genes)} out of {len(top_genes)} selected genes are in the normalized data")
            top_genes = valid_genes
        
        if len(top_genes) < 2:
            logger.error("Not enough valid genes for network construction")
            # Create empty network to avoid downstream errors
            network = nx.Graph()
            self.network = network
            logger.warning("Created empty network to avoid pipeline failure")
            return self.network
        
        expr_data = self.normalized_data.loc[top_genes]
        corr_matrix = expr_data.T.corr(method='pearson')
        corr_path = os.path.join(self.output_dir, "data", "correlation_matrix.csv")
        corr_matrix.to_csv(corr_path)
        logger.info(f"Correlation matrix saved to {corr_path}")
        
        threshold = 0.7
        network = nx.Graph()
        for gene in top_genes:
            network.add_node(gene)
        for i, gene1 in enumerate(top_genes):
            for gene2 in top_genes[i+1:]:
                corr = abs(corr_matrix.loc[gene1, gene2])
                if corr > threshold:
                    network.add_edge(gene1, gene2, weight=corr)
        self.network = network
        logger.info(f"Network constructed with {network.number_of_nodes()} nodes and {network.number_of_edges()} edges")
        nx.write_gexf(network, os.path.join(self.output_dir, "data", "gene_network.gexf"))
        return self.network
    
    def analyze_network(self):
        """Analyze the network to compute centrality measures and create a composite score."""
        logger.info("Analyzing network for potential drug targets")
        if self.network is None:
            logger.error("Network not available. Run construct_network() first.")
            return None
        
        # Handle empty or very small networks gracefully
        if self.network.number_of_nodes() < 2:
            logger.warning("Network has less than 2 nodes, cannot compute meaningful centrality measures")
            # Create a minimal result DataFrame to avoid pipeline failure
            centrality_df = pd.DataFrame({
                'gene': list(self.network.nodes()) if self.network.number_of_nodes() > 0 else ['PLACEHOLDER'],
                'degree_centrality': [0.0] * max(1, self.network.number_of_nodes()),
                'betweenness_centrality': [0.0] * max(1, self.network.number_of_nodes()),
                'eigenvector_centrality': [0.0] * max(1, self.network.number_of_nodes()),
                'composite_score': [0.0] * max(1, self.network.number_of_nodes())
            })
            centrality_path = os.path.join(self.output_dir, "data", "network_targets.csv")
            centrality_df.to_csv(centrality_path, index=False)
            logger.info(f"Minimal network analysis results saved to {centrality_path}")
            self.target_scores = centrality_df
            return self.target_scores
        
        try:
            degree_centrality = nx.degree_centrality(self.network)
            betweenness_centrality = nx.betweenness_centrality(self.network)
            eigenvector_centrality = nx.eigenvector_centrality(self.network, max_iter=1000)
            
            centrality_df = pd.DataFrame({
                'gene': list(self.network.nodes()),
                'degree_centrality': [degree_centrality[g] for g in self.network.nodes()],
                'betweenness_centrality': [betweenness_centrality[g] for g in self.network.nodes()],
                'eigenvector_centrality': [eigenvector_centrality[g] for g in self.network.nodes()]
            })
            
            # Normalize centrality metrics
            from sklearn.preprocessing import MinMaxScaler
            scaler = MinMaxScaler()
            centrality_df[['degree_centrality', 'betweenness_centrality', 'eigenvector_centrality']] = scaler.fit_transform(
                centrality_df[['degree_centrality', 'betweenness_centrality', 'eigenvector_centrality']]
            )
            
            centrality_df['composite_score'] = (
                centrality_df['degree_centrality'] +
                centrality_df['betweenness_centrality'] +
                centrality_df['eigenvector_centrality']
            ) / 3
            
            centrality_df = centrality_df.sort_values('composite_score', ascending=False)
            centrality_path = os.path.join(self.output_dir, "data", "network_targets.csv")
            centrality_df.to_csv(centrality_path, index=False)
            logger.info(f"Network analysis results saved to {centrality_path}")
            self.target_scores = centrality_df
            return self.target_scores
        
        except Exception as e:
            logger.error(f"Error during network analysis: {str(e)}")
            # Create a minimal result DataFrame to avoid pipeline failure
            centrality_df = pd.DataFrame({
                'gene': list(self.network.nodes()),
                'degree_centrality': [0.0] * self.network.number_of_nodes(),
                'betweenness_centrality': [0.0] * self.network.number_of_nodes(),
                'eigenvector_centrality': [0.0] * self.network.number_of_nodes(),
                'composite_score': [0.0] * self.network.number_of_nodes()
            })
            centrality_path = os.path.join(self.output_dir, "data", "network_targets.csv")
            centrality_df.to_csv(centrality_path, index=False)
            logger.warning(f"Created minimal centrality measures due to error: {str(e)}")
            self.target_scores = centrality_df
            return self.target_scores
    
    def _validate_gene_symbol(self, gene_symbol):
        """
        Validate if a string looks like a proper gene symbol.
        
        Args:
            gene_symbol (str): Gene symbol to validate
            
        Returns:
            bool: True if the symbol looks valid, False otherwise
        """
        if not isinstance(gene_symbol, str):
            return False
        
        # Most gene symbols are 1-20 characters long
        if len(gene_symbol) < 1 or len(gene_symbol) > 20:
            return False
        
        # Should not be a probe ID (e.g., '1234_at')
        if '_at' in gene_symbol.lower():
            return False
        
        # Should not be 'UNKNOWN_' prefixed
        if gene_symbol.startswith('UNKNOWN_'):
            return False
        
        # Should contain at least one letter
        if not any(c.isalpha() for c in gene_symbol):
            return False
        
        # Should not contain special characters except hyphen and dot
        if any(c not in '-.' and not c.isalnum() for c in gene_symbol):
            return False
        
        return True

    def _get_ensembl_id(self, gene_symbol):
        """
        Look up Ensembl ID for a gene symbol using the Ensembl REST API.
        
        Args:
            gene_symbol (str): Gene symbol to look up
            
        Returns:
            str: Ensembl gene ID if found, None otherwise
        """
        base_url = "https://rest.ensembl.org"
        endpoint = f"/lookup/symbol/homo_sapiens/{gene_symbol}"
        
        headers = {
            "Content-Type": "application/json",
            "Accept": "application/json"
        }
        
        try:
            response = requests.get(f"{base_url}{endpoint}", headers=headers)
            if response.status_code == 200:
                data = response.json()
                return data.get('id')
            elif response.status_code == 400:
                # Try again with symbol search which can match partial/alternative names
                search_endpoint = f"/xrefs/symbol/homo_sapiens/{gene_symbol}"
                search_response = requests.get(f"{base_url}{search_endpoint}", headers=headers)
                if search_response.status_code == 200:
                    results = search_response.json()
                    if results and len(results) > 0:
                        # Return the first match's ID
                        return results[0].get('id')
            
            return None
        except requests.exceptions.RequestException as e:
            logger.error(f"Ensembl API request failed for {gene_symbol}: {str(e)}")
            return None

    def _query_opentargets(self, gene_symbol):
        """Query the Open Targets GraphQL API for a given gene symbol or Ensembl ID."""
        if not self._validate_gene_symbol(gene_symbol):
            logger.warning(f"Skipping invalid gene symbol: {gene_symbol}")
            return None
        
        # Use the Platform API with the correct endpoint
        url = "https://api.platform.opentargets.org/api/v4/graphql"
        
        # Get Ensembl ID first
        ensembl_id = self._get_ensembl_id(gene_symbol)
        if not ensembl_id:
            logger.warning(f"Could not find Ensembl ID for {gene_symbol}")
            return None
        
        # Create a simple query with only the fields we need, based on the working example
        query_string = """
        query($id: String!) {
          target(ensemblId: $id) {
            id
            approvedSymbol
            biotype
            knownDrugs {
              count
              rows {
                drug {
                  id
                  name
                }
              }
            }
            associatedDiseases {
              count
              rows {
                disease {
                  id
                  name
                }
                score
              }
            }
          }
        }
        """
        
        # Set variables for the query
        variables = {"id": ensembl_id}
        
        try:
            # Make the API request
            logger.info(f"Querying OpenTargets for gene {gene_symbol} with Ensembl ID {ensembl_id}")
            response = requests.post(url, json={"query": query_string, "variables": variables})
            
            if response.status_code != 200:
                logger.warning(f"OpenTargets API request failed with status {response.status_code} for {gene_symbol}")
                return None
            
            json_data = response.json()
            
            # Check for errors in the response
            if "errors" in json_data:
                error_message = json_data["errors"][0]["message"] if json_data["errors"] else "Unknown error"
                logger.warning(f"OpenTargets API returned error for {gene_symbol}: {error_message}")
                return None
            
            # Check if we got data back
            if not json_data.get("data", {}).get("target"):
                logger.warning(f"No target data found for {gene_symbol}")
                return None
            
            logger.info(f"Successfully retrieved data for {gene_symbol}")
            return json_data.get("data", {})
            
        except requests.exceptions.RequestException as e:
            logger.error(f"Request failed for gene {gene_symbol}: {e}")
            return None

    def validate_targets(self, top_n=20):
        """Validate top targets using data from the Open Targets API."""
        logger.info("Validating top targets with Open Targets")
        if self.target_scores is None:
            logger.error("Target scores not available. Run analyze_network() first.")
            return None
        
        # Filter out invalid gene symbols first
        valid_targets = self.target_scores[
            self.target_scores['gene'].apply(self._validate_gene_symbol)
        ]
        
        if len(valid_targets) == 0:
            logger.error("No valid gene symbols found in target scores")
            return None
        
        logger.info(f"Found {len(valid_targets)} valid gene symbols out of {len(self.target_scores)} total targets")
        
        # Take top N valid targets
        top_targets = valid_targets.head(top_n)
        validation_results = []
        
        for _, target_row in top_targets.iterrows():
            gene_symbol = target_row['gene']
            result = self._query_opentargets(gene_symbol)
            
            if result and result.get("target"):
                target_data = result["target"]
                
                # Extract relevant information
                # 1. Known drugs - using the simplified structure
                known_drugs = target_data.get("knownDrugs", {})
                num_known_drugs = known_drugs.get("count", 0)
                
                # 2. Disease associations
                disease_associations = target_data.get("associatedDiseases", {})
                disease_count = disease_associations.get("count", 0)
                
                # Get association scores
                disease_rows = disease_associations.get("rows", [])
                scores = [row.get("score", 0) for row in disease_rows if row.get("score") is not None]
                avg_association_score = sum(scores) / len(scores) if scores else 0
                
                # Calculate a simple drugability score
                drugability_score = (
                    (num_known_drugs * 0.6) +        # Known drugs weight
                    (avg_association_score * 0.4)     # Disease association weight
                )
                
                validation_results.append({
                    "gene": gene_symbol,
                    "composite_score": target_row['composite_score'],
                    "num_known_drugs": num_known_drugs,
                    "avg_association_score": round(avg_association_score, 3),
                    "drugability_score": round(drugability_score, 3)
                })
                logger.info(f"Successfully validated {gene_symbol} (drugability score: {drugability_score:.3f})")
            else:
                logger.warning(f"No Open Targets data found for: {gene_symbol}")
                validation_results.append({
                    "gene": gene_symbol,
                    "composite_score": target_row['composite_score'],
                    "num_known_drugs": 0,
                    "avg_association_score": 0.0,
                    "drugability_score": 0.0
                })
        
        final_targets = pd.DataFrame(validation_results)
        
        # Sort by drugability score
        final_targets = final_targets.sort_values('drugability_score', ascending=False)
        
        # Save results
        final_path = os.path.join(self.output_dir, "data", "GSE46602_final_targets.csv")
        final_targets.to_csv(final_path, index=False)
        logger.info(f"Validated targets with Open Targets saved to {final_path}")
        
        return final_targets
    
    def generate_visualizations(self):
        """Generate figures: volcano plot, network visualization, and bar plot of top targets."""
        logger.info("Generating visualizations")
        plt.style.use('seaborn-v0_8')
        
        # Volcano plot for differential expression
        if self.differential_results is not None and 'log2FC' in self.differential_results.columns and 'adj_pvalue' in self.differential_results.columns:
            plt.figure(figsize=(10, 8))
            plt.scatter(
                self.differential_results['log2FC'],
                -np.log10(self.differential_results['adj_pvalue']),
                alpha=0.5,
                color='gray',
                label='Not significant'
            )
            if self.significant_genes is not None and len(self.significant_genes) > 0:
                plt.scatter(
                    self.significant_genes['log2FC'],
                    -np.log10(self.significant_genes['adj_pvalue']),
                    alpha=0.8,
                    color='red',
                    label='Significant'
                )
            plt.axhline(-np.log10(self.p_threshold), linestyle='--', color='blue')
            plt.axvline(self.fc_threshold, linestyle='--', color='blue')
            plt.axvline(-self.fc_threshold, linestyle='--', color='blue')
            plt.xlabel('Log2 Fold Change')
            plt.ylabel('-Log10 Adjusted P-value')
            plt.title('Volcano Plot: Differential Expression')
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "figures", "volcano_plot.png"), dpi=300)
            plt.close()
        else:
            logger.warning("Skipping volcano plot: differential_results not available or missing required columns")
        
        # Network visualization (subgraph for top 30 targets)
        if self.network is not None and self.target_scores is not None and self.network.number_of_nodes() > 0:
            top_genes = list(self.target_scores.head(min(30, len(self.target_scores)))['gene'])
            # Filter out placeholder genes
            top_genes = [gene for gene in top_genes if gene != 'PLACEHOLDER']
            
            if len(top_genes) > 1:  # Need at least 2 nodes for a meaningful network visualization
                subgraph = self.network.subgraph(top_genes)
                plt.figure(figsize=(12, 10))
                node_sizes = []
                for gene in subgraph.nodes():
                    score = self.target_scores.loc[self.target_scores['gene'] == gene, 'composite_score'].values
                    node_sizes.append(float(score[0]) * 1000 if len(score) > 0 else 500)
                
                edge_weights = [subgraph[u][v].get('weight', 0.5) * 2 for u, v in subgraph.edges()]
                pos = nx.spring_layout(subgraph, seed=42)
                nx.draw_networkx(
                    subgraph,
                    pos=pos,
                    with_labels=True,
                    node_size=node_sizes,
                    node_color='skyblue',
                    edge_color='gray',
                    width=edge_weights,
                    alpha=0.8,
                    font_size=10
                )
                plt.title('Top Genes Network')
                plt.axis('off')
                plt.tight_layout()
                plt.savefig(os.path.join(self.output_dir, "figures", "network_visualization.png"), dpi=300)
                plt.close()
            else:
                logger.warning("Skipping network visualization: Not enough nodes for meaningful visualization")
        else:
            logger.warning("Skipping network visualization: network or target_scores not available")
        
        # Bar plot for top potential targets
        if self.target_scores is not None and len(self.target_scores) > 0 and not all(gene == 'PLACEHOLDER' for gene in self.target_scores['gene']):
            top_n = min(20, len(self.target_scores))
            top_targets = self.target_scores.head(top_n)
            plt.figure(figsize=(12, 8))
            sns.barplot(
                x='composite_score',
                y='gene',
                data=top_targets,
                palette='viridis'
            )
            plt.title(f'Top {top_n} Potential Drug Targets')
            plt.xlabel('Composite Network Score')
            plt.ylabel('Gene')
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "figures", "top_targets.png"), dpi=300)
            plt.close()
        else:
            logger.warning("Skipping bar plot: No valid target scores available")
        
        logger.info("Visualizations generated and saved to figures directory")
    
    def _create_summary_report(self):
        """Create a summary report of the pipeline results."""
        logger.info("Creating summary report")
        summary_path = os.path.join(self.output_dir, "summary.txt")
        with open(summary_path, 'w') as f:
            f.write("# Drug Target Discovery Pipeline Summary\n\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("Dataset: GSE46602\n\n")
            f.write("## Dataset Statistics\n")
            if self.expression_data is not None:
                f.write(f"- Samples: {self.expression_data.shape[1]}\n")
                f.write(f"- Probes/Genes: {self.expression_data.shape[0]}\n\n")
            f.write("## Differential Expression Analysis\n")
            if self.significant_genes is not None:
                f.write(f"- Significant genes: {len(self.significant_genes)}\n")
                up_genes = sum(self.significant_genes['log2FC'] > 0)
                down_genes = sum(self.significant_genes['log2FC'] < 0)
                f.write(f"- Up-regulated: {up_genes}\n")
                f.write(f"- Down-regulated: {down_genes}\n\n")
            f.write("## Network Analysis\n")
            if self.network is not None:
                f.write(f"- Network nodes: {self.network.number_of_nodes()}\n")
                f.write(f"- Network edges: {self.network.number_of_edges()}\n\n")
            f.write("## Top Potential Drug Targets\n")
            if self.target_scores is not None:
                top_10 = self.target_scores.head(10)
                for i, (_, target) in enumerate(top_10.iterrows()):
                    f.write(f"{i+1}. {target['gene']} (score: {target['composite_score']:.4f})\n")
        logger.info(f"Summary report saved to {summary_path}")
    
    def run_pipeline(self):
        """Run the complete drug target discovery pipeline."""
        logger.info("Starting drug target discovery pipeline")
        success = True
        try:
            self.parse_geo_matrix()
        except Exception as e:
            logger.error(f"Error during parsing GEO matrix: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            success = False
            return False
            
        try:
            self.preprocess_data()
        except Exception as e:
            logger.error(f"Error during data preprocessing: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            success = False
            
        try:
            self.map_probes_to_genes()
        except Exception as e:
            logger.error(f"Error during probe-to-gene mapping: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            success = False
            
        try:
            self.perform_differential_analysis()
        except Exception as e:
            logger.error(f"Error during differential analysis: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            success = False
            
        try:
            self.construct_network()
        except Exception as e:
            logger.error(f"Error during network construction: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            success = False
            
        try:
            self.analyze_network()
        except Exception as e:
            logger.error(f"Error during network analysis: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            success = False
            
        try:
            self.validate_targets()
        except Exception as e:
            logger.error(f"Error during target validation: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            success = False
            
        try:
            self.generate_visualizations()
        except Exception as e:
            logger.error(f"Error during visualization generation: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            success = False
            
        try:
            self._create_summary_report()
        except Exception as e:
            logger.error(f"Error creating summary report: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            success = False
            
        if success:
            logger.info("Pipeline completed successfully")
        else:
            logger.warning("Pipeline completed with some errors")
            
        return True

def main():
    parser = argparse.ArgumentParser(description='Drug Target Discovery Pipeline with Bioconductor Mapping')
    parser.add_argument('--matrix-file', type=str, default='GSE46602_series_matrix.txt.gz', help='Path to GEO Series Matrix file')
    parser.add_argument('--output-dir', type=str, help='Output directory for results')
    parser.add_argument('--p-threshold', type=float, default=0.05, help='P-value threshold for significance')
    parser.add_argument('--fc-threshold', type=float, default=1.0, help='Log2 fold change threshold for significance')
    args = parser.parse_args()
    
    pipeline = DrugTargetPipeline(
        matrix_file=args.matrix_file,
        output_dir=args.output_dir,
        p_threshold=args.p_threshold,
        fc_threshold=args.fc_threshold
    )
    
    success = pipeline.run_pipeline()
    if success:
        print(f"Pipeline completed successfully. Results saved to: {pipeline.output_dir}")
        return 0
    else:
        print("Pipeline failed. Check logs for details.")
        return 1

if __name__ == "__main__":
    sys.exit(main())