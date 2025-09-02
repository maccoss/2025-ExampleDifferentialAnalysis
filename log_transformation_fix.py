# FIX FOR LOG TRANSFORMATION ISSUE
# Add this cell to your notebook to fix the large logFC values

# The issue was that the normalization method wasn't being passed to the statistical config
# This prevented the auto log transformation from working properly

print("FIXING LOG TRANSFORMATION ISSUE...")
print(f"Current normalization method: {normalization_method}")

# Re-create the config with proper normalization method
config = ptk.statistical_analysis.StatisticalConfig()

# CRITICAL FIX: Set the normalization method so auto log transformation can work
config.normalization_method = normalization_method.lower()  # Convert to lowercase for consistency

print(f"Config normalization_method set to: {config.normalization_method}")
print(f"Log transformation mode: {config.log_transform_before_stats}")
print(f"Log base: {config.log_base}")

# Re-run the statistical analysis with the fixed config
print("\nRe-running statistical analysis with proper log transformation...")
differential_results = ptk.statistical_analysis.run_comprehensive_statistical_analysis(
    normalized_data=normalized_data,
    sample_metadata=sample_metadata,
    config=config
)

print(f"\n=== FIXED RESULTS ===")
print(f"Number of proteins: {len(differential_results)}")
print(f"logFC range: {differential_results['logFC'].min():.4f} to {differential_results['logFC'].max():.4f}")
print(f"logFC mean: {differential_results['logFC'].mean():.4f}")
print(f"logFC standard deviation: {differential_results['logFC'].std():.4f}")

print("\nFirst few results (should now have reasonable logFC values):")
print(differential_results[['logFC', 'pvalue', 'padj']].head(10))

# Verify the fix worked
max_abs_logfc = abs(differential_results['logFC']).max()
if max_abs_logfc < 5.0:
    print(f"\n✅ SUCCESS: Maximum |logFC| = {max_abs_logfc:.4f} (reasonable for log-transformed data)")
else:
    print(f"\n❌ STILL ISSUE: Maximum |logFC| = {max_abs_logfc:.4f} (still too large)")
