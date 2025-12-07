"""
Tests for proteomics_toolkit.temporal_clustering module

Tests cover:
- TemporalClusteringConfig configuration
- Temporal mean calculation
- Cluster detection and optimization
- Pattern classification
- Integration with statistical results
"""

import pandas as pd
import numpy as np
import pytest
import warnings

from proteomics_toolkit.temporal_clustering import (
    TemporalClusteringConfig,
    calculate_temporal_means,
    get_week_columns,
    determine_optimal_clusters,
    cluster_temporal_trends,
    name_clusters_by_pattern,
    classify_trend_pattern,
    merge_with_statistics,
    filter_significant_proteins,
)


class TestTemporalClusteringConfig:
    """Test the TemporalClusteringConfig class"""

    def test_config_initialization(self):
        """Test configuration initialization with defaults"""
        config = TemporalClusteringConfig()

        assert config.n_clusters == 4
        assert config.auto_detect_clusters is True
        assert config.min_clusters == 2
        assert config.max_clusters == 8
        assert config.clustering_method == 'kmeans'
        assert config.random_seed == 42
        assert config.p_value_threshold == 0.05
        assert config.use_adjusted_pvalue is False
        assert config.min_fold_change == 0.0
        assert len(config.enrichr_libraries) > 0
        assert 'GO_Biological_Process_2023' in config.enrichr_libraries

    def test_config_custom_values(self):
        """Test configuration with custom values"""
        config = TemporalClusteringConfig(
            n_clusters=6,
            auto_detect_clusters=False,
            p_value_threshold=0.01,
            use_adjusted_pvalue=True,
            min_fold_change=0.5,
        )

        assert config.n_clusters == 6
        assert config.auto_detect_clusters is False
        assert config.p_value_threshold == 0.01
        assert config.use_adjusted_pvalue is True
        assert config.min_fold_change == 0.5


class TestCalculateTemporalMeans:
    """Test temporal mean calculation"""

    @pytest.fixture
    def temporal_data(self):
        """Create sample temporal data"""
        np.random.seed(42)
        
        # Create protein data with sample columns
        data = pd.DataFrame({
            'Protein': ['P00001', 'P00002', 'P00003'],
            'Gene': ['GENE1', 'GENE2', 'GENE3'],
            # Subject 1 samples
            'S1_W0': [100.0, 200.0, 300.0],
            'S1_W2': [110.0, 180.0, 320.0],
            'S1_W4': [120.0, 160.0, 340.0],
            # Subject 2 samples
            'S2_W0': [90.0, 210.0, 280.0],
            'S2_W2': [100.0, 190.0, 300.0],
            'S2_W4': [110.0, 170.0, 320.0],
        })
        
        # Create metadata dictionary
        metadata = {
            'S1_W0': {'Subject': 'S1', 'Week': 0},
            'S1_W2': {'Subject': 'S1', 'Week': 2},
            'S1_W4': {'Subject': 'S1', 'Week': 4},
            'S2_W0': {'Subject': 'S2', 'Week': 0},
            'S2_W2': {'Subject': 'S2', 'Week': 2},
            'S2_W4': {'Subject': 'S2', 'Week': 4},
        }
        
        return data, metadata

    def test_calculate_temporal_means_basic(self, temporal_data):
        """Test basic temporal mean calculation"""
        data, metadata = temporal_data
        
        result, unique_weeks = calculate_temporal_means(
            data_df=data,
            metadata_dict=metadata,
            week_column='Week',
            subject_column='Subject'
        )
        
        assert isinstance(result, pd.DataFrame)
        assert 'Protein' in result.columns
        assert 'Gene' in result.columns
        assert len(unique_weeks) == 3
        assert set(unique_weeks) == {0, 2, 4}
        assert len(result) == 3  # 3 proteins

    def test_calculate_temporal_means_week_columns(self, temporal_data):
        """Test that week columns are created correctly"""
        data, metadata = temporal_data
        
        result, unique_weeks = calculate_temporal_means(
            data_df=data,
            metadata_dict=metadata,
            week_column='Week',
            subject_column='Subject'
        )
        
        week_cols = get_week_columns(result)
        assert len(week_cols) == 3
        assert 'Week_0' in week_cols
        assert 'Week_2' in week_cols
        assert 'Week_4' in week_cols


class TestDetermineOptimalClusters:
    """Test optimal cluster detection"""

    @pytest.fixture
    def clusterable_data(self):
        """Create data with clear cluster structure"""
        np.random.seed(42)
        n_proteins = 30
        n_timepoints = 4
        
        # Create 3 distinct patterns
        pattern1 = np.linspace(0, 1, n_timepoints)  # Increasing
        pattern2 = np.linspace(1, 0, n_timepoints)  # Decreasing
        pattern3 = np.array([0, 1, 1, 0])  # Transient
        
        # Generate data with noise
        X = np.vstack([
            pattern1 + np.random.normal(0, 0.1, (10, n_timepoints)),
            pattern2 + np.random.normal(0, 0.1, (10, n_timepoints)),
            pattern3 + np.random.normal(0, 0.1, (10, n_timepoints)),
        ])
        
        return X

    def test_determine_optimal_clusters_basic(self, clusterable_data):
        """Test basic cluster detection"""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            optimal_k, fig, scores = determine_optimal_clusters(
                X_scaled=clusterable_data,
                k_range=(2, 6),
                random_seed=42,
                plot=False
            )
        
        assert isinstance(optimal_k, int)
        assert 2 <= optimal_k <= 6
        assert 'k_values' in scores
        assert 'silhouette_scores' in scores
        assert 'inertias' in scores
        assert 'selected_k' in scores
        assert scores['selected_k'] == optimal_k

    def test_determine_optimal_clusters_with_plot(self, clusterable_data):
        """Test cluster detection with plot generation"""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            optimal_k, fig, scores = determine_optimal_clusters(
                X_scaled=clusterable_data,
                k_range=(2, 5),
                random_seed=42,
                plot=True
            )
        
        assert fig is not None
        # Figure should have 2 subplots (elbow and silhouette)


class TestClusterTemporalTrends:
    """Test clustering of temporal trends"""

    @pytest.fixture
    def temporal_means_df(self):
        """Create temporal means dataframe"""
        np.random.seed(42)
        return pd.DataFrame({
            'Protein': [f'P{i:05d}' for i in range(20)],
            'Gene': [f'GENE{i}' for i in range(20)],
            'Week_0': np.random.randn(20),
            'Week_2': np.random.randn(20),
            'Week_4': np.random.randn(20),
            'Week_6': np.random.randn(20),
        })

    def test_cluster_temporal_trends_kmeans(self, temporal_means_df):
        """Test K-means clustering"""
        config = TemporalClusteringConfig(
            n_clusters=3,
            auto_detect_clusters=False,
            clustering_method='kmeans',
            random_seed=42
        )
        
        week_cols = get_week_columns(temporal_means_df)
        labels, X_scaled, model, fig = cluster_temporal_trends(
            temporal_df=temporal_means_df,
            week_columns=week_cols,
            config=config
        )
        
        assert isinstance(labels, np.ndarray)
        assert len(np.unique(labels)) == 3
        assert len(labels) == 20

    def test_cluster_temporal_trends_auto_detect(self, temporal_means_df):
        """Test clustering with auto cluster detection"""
        config = TemporalClusteringConfig(
            auto_detect_clusters=True,
            min_clusters=2,
            max_clusters=5,
            random_seed=42
        )
        
        week_cols = get_week_columns(temporal_means_df)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            labels, X_scaled, model, fig = cluster_temporal_trends(
                temporal_df=temporal_means_df,
                week_columns=week_cols,
                config=config
            )
        
        assert isinstance(labels, np.ndarray)
        n_clusters = len(np.unique(labels))
        assert 2 <= n_clusters <= 5


class TestPatternClassification:
    """Test pattern classification functions"""

    def test_classify_trend_pattern_increasing(self):
        """Test classification of increasing trend"""
        # Strong increasing trend
        trend = np.array([-1.0, -0.3, 0.3, 1.0])
        pattern_name, confidence = classify_trend_pattern(trend)
        
        # classify_trend_pattern returns (name, confidence) tuple
        # Pattern names include: 'Up & Stay Up', 'Down & Stay Down', 'Up then Down', etc.
        assert isinstance(pattern_name, str)
        assert isinstance(confidence, (int, float))
        assert confidence >= 0
        # Increasing patterns tend to have "Up" in the name
        assert 'Up' in pattern_name or 'Increase' in pattern_name.lower()

    def test_classify_trend_pattern_decreasing(self):
        """Test classification of decreasing trend"""
        # Strong decreasing trend
        trend = np.array([1.0, 0.3, -0.3, -1.0])
        pattern_name, confidence = classify_trend_pattern(trend)
        
        # Pattern names include: 'Down & Stay Down', etc.
        assert isinstance(pattern_name, str)
        assert isinstance(confidence, (int, float))
        assert confidence >= 0
        # Decreasing patterns tend to have "Down" in the name
        assert 'Down' in pattern_name or 'Decrease' in pattern_name.lower()

    def test_classify_trend_pattern_transient(self):
        """Test classification of transient pattern"""
        # Transient pattern (up then down)
        trend = np.array([-0.5, 0.8, 0.8, -0.5])
        pattern_name, confidence = classify_trend_pattern(trend)
        
        # Should return a tuple (name, confidence)
        assert isinstance(pattern_name, str)
        assert isinstance(confidence, (int, float))
        # Transient patterns may include "Up then Down" or "Down then Up"
        assert len(pattern_name) > 0

    def test_name_clusters_by_pattern(self):
        """Test cluster naming based on patterns"""
        np.random.seed(42)
        
        # Create Z-scored data matrix (proteins x timepoints)
        X_scaled = np.array([
            [-1.0, -0.5, 0.5, 1.0],   # Increasing
            [1.0, 0.5, -0.5, -1.0],   # Decreasing
            [-0.5, 0.8, 0.8, -0.5],   # Transient
            [0.5, 0.0, -0.5, 0.0],    # Mixed
        ])
        
        labels = np.array([0, 1, 2, 3])
        week_cols = ['Week_0', 'Week_2', 'Week_4', 'Week_6']
        
        cluster_names = name_clusters_by_pattern(X_scaled, labels, week_cols)
        
        assert isinstance(cluster_names, dict)
        assert len(cluster_names) == 4


class TestMergeWithStatistics:
    """Test merging temporal data with statistical results"""

    @pytest.fixture
    def temporal_df(self):
        """Create temporal dataframe"""
        return pd.DataFrame({
            'Protein': ['P00001', 'P00002', 'P00003'],
            'Gene': ['GENE1', 'GENE2', 'GENE3'],
            'Week_0': [0.0, 0.0, 0.0],
            'Week_4': [1.0, -1.0, 0.5],
        })

    @pytest.fixture
    def stats_df(self):
        """Create statistics dataframe"""
        return pd.DataFrame({
            'Protein': ['P00001', 'P00002', 'P00003', 'P00004'],
            'Gene': ['GENE1', 'GENE2', 'GENE3', 'GENE4'],
            'logFC': [1.5, -1.2, 0.3, 0.8],
            'P.Value': [0.001, 0.01, 0.5, 0.03],
            'adj.P.Val': [0.01, 0.05, 0.8, 0.1],
        })

    def test_merge_with_statistics(self, temporal_df, stats_df):
        """Test basic merge"""
        cluster_labels = np.array([0, 1, 0])
        cluster_names = {0: 'Cluster 1', 1: 'Cluster 2'}
        
        merged = merge_with_statistics(temporal_df, stats_df, cluster_labels, cluster_names)
        
        assert 'P.Value' in merged.columns
        assert 'logFC' in merged.columns
        assert 'Cluster' in merged.columns
        assert 'Cluster_Name' in merged.columns
        assert len(merged) == 3  # 3 proteins in temporal_df

    def test_filter_significant_proteins(self, temporal_df, stats_df):
        """Test filtering by significance"""
        cluster_labels = np.array([0, 1, 0])
        cluster_names = {0: 'Cluster 1', 1: 'Cluster 2'}
        
        merged = merge_with_statistics(temporal_df, stats_df, cluster_labels, cluster_names)
        
        config = TemporalClusteringConfig(
            p_value_threshold=0.05,
            use_adjusted_pvalue=False
        )
        
        filtered = filter_significant_proteins(merged, config)
        
        # Should only keep proteins with P.Value < 0.05
        assert len(filtered) <= len(merged)
        if len(filtered) > 0:
            assert all(filtered['P.Value'] < 0.05)


class TestIntegration:
    """Integration tests for the temporal clustering pipeline"""

    @pytest.fixture
    def complete_test_data(self):
        """Create complete test dataset"""
        np.random.seed(42)
        n_proteins = 50
        
        # Create protein data
        data = pd.DataFrame({
            'Protein': [f'P{i:05d}' for i in range(n_proteins)],
            'Gene': [f'GENE{i}' for i in range(n_proteins)],
        })
        
        # Add sample columns for 3 subjects, 4 timepoints
        subjects = ['S1', 'S2', 'S3']
        weeks = [0, 2, 4, 8]
        
        metadata = {}
        for subj in subjects:
            for week in weeks:
                col_name = f'{subj}_W{week}'
                # Generate realistic-looking abundance data
                data[col_name] = np.random.lognormal(mean=10, sigma=1, size=n_proteins)
                metadata[col_name] = {'Subject': subj, 'Week': week}
        
        # Create stats dataframe
        stats = pd.DataFrame({
            'Protein': data['Protein'],
            'Gene': data['Gene'],
            'logFC': np.random.normal(0, 1, n_proteins),
            'P.Value': np.random.uniform(0, 1, n_proteins),
            'adj.P.Val': np.random.uniform(0, 1, n_proteins),
        })
        
        return data, metadata, stats

    def test_full_pipeline(self, complete_test_data):
        """Test the full temporal clustering pipeline"""
        data, metadata, stats = complete_test_data
        
        # Calculate temporal means
        temporal_means, weeks = calculate_temporal_means(
            data_df=data,
            metadata_dict=metadata,
            week_column='Week',
            subject_column='Subject'
        )
        
        assert len(weeks) == 4
        assert len(temporal_means) == 50
        
        # Configure and run clustering
        config = TemporalClusteringConfig(
            n_clusters=3,
            auto_detect_clusters=False,
            random_seed=42
        )
        
        week_cols = get_week_columns(temporal_means)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            labels, X_scaled, model, fig = cluster_temporal_trends(
                temporal_means, week_cols, config
            )
        
        assert isinstance(labels, np.ndarray)
        assert len(labels) == 50
        
        # Name clusters
        cluster_names = name_clusters_by_pattern(X_scaled, labels, week_cols)
        assert isinstance(cluster_names, dict)
        
        # Merge with statistics
        merged = merge_with_statistics(temporal_means, stats, labels, cluster_names)
        assert 'P.Value' in merged.columns
        assert 'Cluster' in merged.columns
        
        # Filter significant
        config_sig = TemporalClusteringConfig(p_value_threshold=0.5)  # High threshold for test
        filtered = filter_significant_proteins(merged, config_sig)
        
        assert len(filtered) > 0


class TestEdgeCases:
    """Test edge cases and error handling"""

    def test_empty_metadata(self):
        """Test with empty metadata"""
        data = pd.DataFrame({
            'Protein': ['P00001'],
            'Gene': ['GENE1'],
            'Sample1': [100.0],
        })
        
        temporal_means, weeks = calculate_temporal_means(
            data_df=data,
            metadata_dict={},
            week_column='Week',
            subject_column='Subject'
        )
        
        # Should handle gracefully
        assert len(weeks) == 0 or len(temporal_means) == 1

    def test_single_timepoint(self):
        """Test with single timepoint"""
        data = pd.DataFrame({
            'Protein': ['P00001', 'P00002'],
            'Gene': ['GENE1', 'GENE2'],
            'S1_W0': [100.0, 200.0],
            'S2_W0': [110.0, 210.0],
        })
        
        metadata = {
            'S1_W0': {'Subject': 'S1', 'Week': 0},
            'S2_W0': {'Subject': 'S2', 'Week': 0},
        }
        
        temporal_means, weeks = calculate_temporal_means(
            data_df=data,
            metadata_dict=metadata,
            week_column='Week',
            subject_column='Subject'
        )
        
        assert len(weeks) == 1
        assert weeks[0] == 0

    def test_missing_subject_data(self):
        """Test handling of incomplete subject data"""
        data = pd.DataFrame({
            'Protein': ['P00001'],
            'Gene': ['GENE1'],
            'S1_W0': [100.0],
            'S1_W4': [120.0],
            # S2 only has one timepoint
            'S2_W0': [90.0],
        })
        
        metadata = {
            'S1_W0': {'Subject': 'S1', 'Week': 0},
            'S1_W4': {'Subject': 'S1', 'Week': 4},
            'S2_W0': {'Subject': 'S2', 'Week': 0},
        }
        
        temporal_means, weeks = calculate_temporal_means(
            data_df=data,
            metadata_dict=metadata,
            week_column='Week',
            subject_column='Subject'
        )
        
        # Should still work with partial data
        assert len(weeks) == 2
