import pandas as pd
import numpy as np
from scipy import stats

# 1. Load the merged PCA dataset
df = pd.read_csv("pca/pca_with_metadata.csv")

# Identify our PCA coordinate coordinates and potential metadata features
pc_cols = [f"PC{i}" for i in range(1, 6)]
metadata_cols = [col for col in df.columns if col not in pc_cols and col not in ["Sample_Full_ID", "BQC ID"]]

results = []

print("Scanning all clinical features for optimal group splits...")

for col in metadata_cols:
    # Drop rows missing this specific clinical data point
    clean_df = df[[col, "PC1", "PC2"]].dropna()

    # Skip columns that don't have enough data variability
    if len(clean_df) < 10 or clean_df[col].nunique() < 2:
        continue

    # Check data type to determine the appropriate statistical model
    is_numeric = pd.api.types.is_numeric_dtype(clean_df[col])
    unique_vals = clean_df[col].nunique()

    # Treat as Categorical if it's text or has very few distinct integer states (like checkboxes/groups)
    if not is_numeric or unique_vals <= 5:
        # Group PC1 values by the category labels
        groups_pc1 = [group["PC1"].values for name, group in clean_df.groupby(col)]
        groups_pc2 = [group["PC2"].values for name, group in clean_df.groupby(col)]

        # Protect against empty categories
        if any(len(g) < 2 for g in groups_pc1): continue

        # Calculate ANOVA F-statistic (Higher F means cleaner, more separate groups)
        f_stat_pc1, p_val_pc1 = stats.f_oneway(*groups_pc1)
        f_stat_pc2, p_val_pc2 = stats.f_oneway(*groups_pc2)

        results.append({
            "Clinical_Feature": col,
            "Type": "Categorical/Discrete",
            "Unique_Groups": unique_vals,
            "PC1_Metric_Value": f_stat_pc1,
            "PC1_p_value": p_val_pc1,
            "PC2_Metric_Value": f_stat_pc2,
            "PC2_p_value": p_val_pc2,
            "Ranking_Score": f_stat_pc1  # Primary rank based on group separation power
        })

    else:
        # Treat as Continuous: Compute Pearson correlation coefficient
        r_pc1, p_val_pc1 = stats.pearsonr(clean_df[col], clean_df["PC1"])
        r_pc2, p_val_pc2 = stats.pearsonr(clean_df[col], clean_df["PC2"])

        # Convert correlation strength to an absolute scale matching F-statistic intent
        results.append({
            "Clinical_Feature": col,
            "Type": "Continuous Metric",
            "Unique_Groups": "N/A (Continuous)",
            "PC1_Metric_Value": r_pc1,
            "PC1_p_value": p_val_pc1,
            "PC2_Metric_Value": r_pc2,
            "PC2_p_value": p_val_pc2,
            "Ranking_Score": abs(r_pc1) * 100  # Scaled variance score for uniform sorting
        })

# 2. Compile results and sort by the strongest splitting power on PC1
summary_df = pd.DataFrame(results)
summary_df = summary_df.sort_values(by="Ranking_Score", ascending=False)

# Save the complete statistical breakdown
summary_df.to_csv("pca/clinical_split_rankings.tsv", sep="\t", index=False)

# 3. Print the top 10 best columns for splitting your PCA to screen
print("\n" + "=" * 85)
print("🏆 TOP 10 CLINICAL COLUMNS THAT SPLIT YOUR GENOMIC GROUPS BEST ON PC1 🏆")
print("=" * 85)
for idx, row in summary_df.head(10).iterrows():
    metric_label = "F-Stat" if row["Type"] == "Categorical/Discrete" else "Corr (r)"
    print(f"📍 Feature: {row['Clinical_Feature']}")
    print(f"   - Variable Type : {row['Type']}")
    print(f"   - PC1 Separation: {metric_label} = {row['PC1_Metric_Value']:.4f} (p = {row['PC1_p_value']:.2e})")
    print(f"   - PC2 Separation: {metric_label} = {row['PC2_Metric_Value']:.4f} (p = {row['PC2_p_value']:.2e})")
    print("-" * 85)
