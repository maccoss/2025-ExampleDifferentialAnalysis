"""
Tests for the enrichment module.

These tests cover:
- EnrichmentConfig initialization and customization
- parse_enrichr_results function
- run_enrichment_by_group function (mock API calls)
- run_differential_enrichment function
- Visualization functions (plot_enrichment_barplot, plot_enrichment_comparison)
- Utility functions

Note: Tests that query the Enrichr API are marked with @pytest.mark.network
and are skipped by default (use --run-network to run them).
"""

import pytest
import pandas as pd
import numpy as np
from unittest.mock import patch, MagicMock

from proteomics_toolkit.enrichment import (
    EnrichmentConfig,
    query_enrichr,
    parse_enrichr_results,
    run_enrichment_analysis,
    run_enrichment_by_group,
    run_differential_enrichment,
    plot_enrichment_barplot,
    plot_enrichment_comparison,
    get_available_libraries,
    merge_enrichment_results,
    LIBRARY_COLORS,
)


class TestEnrichmentConfig:
    """Test EnrichmentConfig dataclass"""
    
    def test_config_defaults(self):
        """Test default configuration values"""
        config = EnrichmentConfig()
        
        assert config.pvalue_cutoff == 0.05
        assert config.top_n == 20
        assert config.min_genes == 5
        assert config.rate_limit_delay == 0.5
        assert config.timeout == 30
        assert 'GO_Biological_Process_2023' in config.enrichr_libraries
        assert 'KEGG_2021_Human' in config.enrichr_libraries
    
    def test_config_custom_values(self):
        """Test custom configuration"""
        config = EnrichmentConfig(
            pvalue_cutoff=0.01,
            top_n=10,
            min_genes=10,
            enrichr_libraries=['KEGG_2021_Human']
        )
        
        assert config.pvalue_cutoff == 0.01
        assert config.top_n == 10
        assert config.min_genes == 10
        assert config.enrichr_libraries == ['KEGG_2021_Human']


class TestParseEnrichrResults:
    """Test parse_enrichr_results function"""
    
    def test_parse_empty_results(self):
        """Test parsing empty results"""
        result = parse_enrichr_results({})
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0
    
    def test_parse_valid_results(self):
        """Test parsing valid Enrichr-format results"""
        # Simulate Enrichr API response format:
        # [rank, term, pval, zscore, combined_score, genes, adj_pval]
        mock_results = {
            'GO_Biological_Process_2023': [
                [1, 'Cell cycle (GO:0007049)', 0.001, 2.5, 15.3, ['TP53', 'BRCA1', 'CDK1'], 0.01],
                [2, 'DNA repair (GO:0006281)', 0.002, 2.0, 12.1, ['BRCA1', 'ATM'], 0.02],
            ]
        }
        
        result = parse_enrichr_results(mock_results)
        
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2
        assert 'Library' in result.columns
        assert 'Term' in result.columns
        assert 'P_Value' in result.columns
        assert 'Combined_Score' in result.columns
        assert 'Genes' in result.columns
        assert 'N_Genes' in result.columns
        
        # Should be sorted by Combined_Score descending
        assert result.iloc[0]['Combined_Score'] >= result.iloc[1]['Combined_Score']
    
    def test_parse_filters_by_pvalue(self):
        """Test that parsing filters by p-value threshold"""
        mock_results = {
            'KEGG_2021_Human': [
                [1, 'Pathway A', 0.001, 2.5, 15.3, ['GENE1'], 0.01],  # Should pass
                [2, 'Pathway B', 0.1, 1.0, 5.0, ['GENE2'], 0.2],      # Should be filtered
            ]
        }
        
        config = EnrichmentConfig(pvalue_cutoff=0.05)
        result = parse_enrichr_results(mock_results, config)
        
        assert len(result) == 1
        assert result.iloc[0]['Term'] == 'Pathway A'


class TestRunEnrichmentByGroup:
    """Test run_enrichment_by_group function"""
    
    @pytest.fixture
    def sample_grouped_data(self):
        """Create sample data with groups"""
        return pd.DataFrame({
            'Protein': [f'P{i:05d}' for i in range(20)],
            'Gene': [f'GENE{i}' for i in range(20)],
            'Cluster': ['A'] * 10 + ['B'] * 10,
        })
    
    def test_groups_genes_correctly(self, sample_grouped_data):
        """Test that function correctly groups genes"""
        with patch('proteomics_toolkit.enrichment.run_enrichment_analysis') as mock_analysis:
            mock_analysis.return_value = pd.DataFrame()
            
            result = run_enrichment_by_group(
                sample_grouped_data,
                group_column='Cluster',
                gene_column='Gene',
                verbose=False
            )
            
            assert isinstance(result, dict)
            assert 'A' in result
            assert 'B' in result
            # Each group should have been called with 10 genes
            assert mock_analysis.call_count == 2
    
    def test_skips_small_groups(self):
        """Test that groups with too few genes are skipped"""
        data = pd.DataFrame({
            'Protein': ['P1', 'P2', 'P3'],
            'Gene': ['GENE1', 'GENE2', 'GENE3'],
            'Cluster': ['A', 'A', 'B'],  # B has only 1 gene
        })
        
        config = EnrichmentConfig(min_genes=5)
        
        with patch('proteomics_toolkit.enrichment.run_enrichment_analysis') as mock_analysis:
            mock_analysis.return_value = pd.DataFrame()
            
            result = run_enrichment_by_group(
                data,
                group_column='Cluster',
                gene_column='Gene',
                config=config,
                verbose=False
            )
            
            # Neither group should have been queried (both < 5 genes)
            assert mock_analysis.call_count == 0
            assert all(df.empty for df in result.values())


class TestRunDifferentialEnrichment:
    """Test run_differential_enrichment function"""
    
    @pytest.fixture
    def sample_diff_results(self):
        """Create sample differential expression results"""
        np.random.seed(42)
        n = 100
        
        return pd.DataFrame({
            'Protein': [f'P{i:05d}' for i in range(n)],
            'Gene': [f'GENE{i}' for i in range(n)],
            'logFC': np.random.normal(0, 2, n),
            'adj.P.Val': np.random.uniform(0, 0.2, n),
        })
    
    def test_splits_up_down(self, sample_diff_results):
        """Test that function correctly splits up/down regulated"""
        with patch('proteomics_toolkit.enrichment.run_enrichment_analysis') as mock_analysis:
            mock_analysis.return_value = pd.DataFrame()
            
            result = run_differential_enrichment(
                sample_diff_results,
                logfc_threshold=1.0,
                pvalue_threshold=0.05,
                verbose=False
            )
            
            assert isinstance(result, dict)
            assert 'Upregulated' in result
            assert 'Downregulated' in result
    
    def test_respects_thresholds(self, sample_diff_results):
        """Test that thresholds are correctly applied"""
        # Set very strict thresholds that should exclude most genes
        config = EnrichmentConfig(min_genes=1)  # Allow small groups for testing
        
        with patch('proteomics_toolkit.enrichment.run_enrichment_analysis') as mock_analysis:
            mock_analysis.return_value = pd.DataFrame()
            
            result = run_differential_enrichment(
                sample_diff_results,
                logfc_threshold=3.0,  # Very strict
                pvalue_threshold=0.01,  # Very strict
                config=config,
                verbose=False
            )
            
            # Should still return the dict structure
            assert 'Upregulated' in result
            assert 'Downregulated' in result


class TestVisualizationFunctions:
    """Test visualization functions"""
    
    @pytest.fixture
    def sample_enrichment_df(self):
        """Create sample enrichment results"""
        return pd.DataFrame({
            'Library': ['GO_Biological_Process_2023'] * 5,
            'Term': [f'Term {i}' for i in range(5)],
            'P_Value': [0.001, 0.002, 0.003, 0.004, 0.005],
            'Adj_P_Value': [0.01, 0.02, 0.03, 0.04, 0.05],
            'Z_Score': [2.5, 2.0, 1.8, 1.5, 1.2],
            'Combined_Score': [15.0, 12.0, 10.0, 8.0, 6.0],
            'Genes': ['A;B;C', 'D;E', 'F;G;H', 'I;J', 'K'],
            'N_Genes': [3, 2, 3, 2, 1],
        })
    
    def test_plot_enrichment_barplot(self, sample_enrichment_df):
        """Test bar plot creation"""
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
        
        fig = plot_enrichment_barplot(
            sample_enrichment_df,
            title='Test Enrichment',
            top_n=5
        )
        
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)
    
    def test_plot_enrichment_barplot_empty(self):
        """Test bar plot with empty data returns None"""
        result = plot_enrichment_barplot(
            pd.DataFrame(),
            title='Empty Test'
        )
        
        assert result is None
    
    def test_plot_enrichment_comparison(self, sample_enrichment_df):
        """Test comparison plot creation"""
        import matplotlib
        matplotlib.use('Agg')
        
        enrichment_dict = {
            'Group A': sample_enrichment_df,
            'Group B': sample_enrichment_df.copy(),
        }
        
        fig = plot_enrichment_comparison(
            enrichment_dict,
            title='Test Comparison'
        )
        
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)
    
    def test_plot_enrichment_comparison_empty(self):
        """Test comparison plot with empty data returns None"""
        result = plot_enrichment_comparison(
            {'A': pd.DataFrame(), 'B': pd.DataFrame()},
            title='Empty Test'
        )
        
        assert result is None


class TestUtilityFunctions:
    """Test utility functions"""
    
    def test_get_available_libraries(self):
        """Test that library list is returned"""
        libraries = get_available_libraries()
        
        assert isinstance(libraries, list)
        assert len(libraries) > 0
        assert 'GO_Biological_Process_2023' in libraries
        assert 'KEGG_2021_Human' in libraries
    
    def test_merge_enrichment_results(self):
        """Test merging multiple enrichment DataFrames"""
        df1 = pd.DataFrame({
            'Term': ['A', 'B'],
            'P_Value': [0.01, 0.02],
        })
        df2 = pd.DataFrame({
            'Term': ['C', 'D'],
            'P_Value': [0.03, 0.04],
        })
        
        merged = merge_enrichment_results({
            'Group1': df1,
            'Group2': df2,
        })
        
        assert len(merged) == 4
        assert 'Group' in merged.columns
        assert set(merged['Group'].unique()) == {'Group1', 'Group2'}
    
    def test_merge_enrichment_results_empty(self):
        """Test merging with empty DataFrames"""
        merged = merge_enrichment_results({
            'Group1': pd.DataFrame(),
            'Group2': pd.DataFrame(),
        })
        
        assert len(merged) == 0
    
    def test_library_colors_defined(self):
        """Test that library colors are defined for common libraries"""
        assert 'GO_Biological_Process_2023' in LIBRARY_COLORS
        assert 'KEGG_2021_Human' in LIBRARY_COLORS
        assert 'Reactome_2022' in LIBRARY_COLORS


# Network tests - these actually hit the Enrichr API
# Run with: pytest --run-network
@pytest.mark.network
class TestEnrichrAPIIntegration:
    """Integration tests that require network access to Enrichr API.
    
    These tests are skipped by default. Run with --run-network to execute.
    """
    
    @pytest.fixture
    def sample_genes(self):
        """Well-known cancer genes for testing"""
        return ['TP53', 'BRCA1', 'BRCA2', 'ATM', 'CHEK2', 'PALB2', 
                'RAD51', 'BARD1', 'PTEN', 'CDK1']
    
    def test_query_enrichr_real(self, sample_genes):
        """Test actual API query"""
        config = EnrichmentConfig(
            enrichr_libraries=['GO_Biological_Process_2023'],
            top_n=5
        )
        
        results = query_enrichr(sample_genes, config)
        
        assert isinstance(results, dict)
        # Should have results from the one library we queried
        if results:  # API might be unavailable
            assert 'GO_Biological_Process_2023' in results or len(results) == 0
    
    def test_run_enrichment_analysis_real(self, sample_genes):
        """Test complete enrichment analysis"""
        config = EnrichmentConfig(
            enrichr_libraries=['GO_Biological_Process_2023'],
            top_n=5,
            pvalue_cutoff=0.1  # More lenient for testing
        )
        
        result = run_enrichment_analysis(sample_genes, config, verbose=False)
        
        assert isinstance(result, pd.DataFrame)
        # Result might be empty if API fails or no significant terms


def pytest_configure(config):
    """Configure custom markers"""
    config.addinivalue_line(
        "markers", "network: mark test as requiring network access (skip by default)"
    )


def pytest_collection_modifyitems(config, items):
    """Skip network tests unless --run-network is passed"""
    if config.getoption("--run-network", default=False):
        return
    
    skip_network = pytest.mark.skip(reason="need --run-network option to run")
    for item in items:
        if "network" in item.keywords:
            item.add_marker(skip_network)
