import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load the integrated dataset
df = pd.read_csv("pca/pca_with_metadata.csv")

# Standardize string representations to avoid whitespace errors
who_col = "If a screening test for SARS-CoV-2 by PCR was performed, what is the most severe severity level (according to WHO) achieved?"
diabetes_col = "Diabetes ?"

for c in [who_col, diabetes_col]:
    if c in df.columns:
        df[c] = df[c].astype(str).str.strip()

# Clean missing records and filter out generic strings
clean_df = df.dropna(subset=[diabetes_col, who_col])
clean_df = clean_df[~clean_df[who_col].str.lower().isin(["nan", "unknown", "missing"])]
clean_df = clean_df[~clean_df[diabetes_col].str.lower().isin(["nan", "unknown", "missing"])]

# Fix a uniform category sorting for the WHO Levels
who_categories = sorted(clean_df[who_col].unique())

# Set up a 2x2 grid canvas layout
fig, axes = plt.subplots(2, 2, figsize=(24, 18), dpi=300)

# =====================================================================
# ROW 1: THE OVERALL COHORT PLOTS
# =====================================================================

# Panel A: Baseline Comorbidity (Diabetes)
sns.scatterplot(
    data=clean_df,
    x="PC1", y="PC2",
    hue=diabetes_col,
    palette="Set1",
    alpha=0.6, s=25, edgecolor=None,
    ax=axes[0, 0]
)
axes[0, 0].set_title("PCA Split by Baseline Comorbidity\n(Diabetes ?)", fontsize=13, weight='bold')
axes[0, 0].set_xlabel("PC1 (37.62% Var) - Tracks Chronic Health Status", fontsize=11)
axes[0, 0].set_ylabel("PC2 (10.38% Var)", fontsize=11)
axes[0, 0].legend(title="Diabetes Status", loc='upper right')

# Panel B: Overall Acute Presentation (WHO Severity Level)
sns.scatterplot(
    data=clean_df,
    x="PC1", y="PC2",
    hue=who_col,
    hue_order=who_categories,
    palette="viridis",
    alpha=0.6, s=25, edgecolor=None,
    ax=axes[0, 1]
)
axes[0, 1].set_title("SARS-CoV-2 severe severity level (according to WHO)", fontsize=13, weight='bold')
axes[0, 1].set_xlabel("PC1 (37.62% Var)", fontsize=11)
axes[0, 1].set_ylabel("PC2 (10.38% Var) - Tracks Infection Severity", fontsize=11)
axes[0, 1].legend(title="WHO Level", bbox_to_anchor=(1.02, 1), loc='upper left')


# =====================================================================
# ROW 2: THE STRATIFIED COHORT PLOTS
# =====================================================================

# Panel C: Patients WITH Diabetes (Diabetes Status = Yes)
df_diab_yes = clean_df[clean_df[diabetes_col].str.lower() == "yes"]
sns.scatterplot(
    data=df_diab_yes,
    x="PC1", y="PC2",
    hue=who_col,
    hue_order=who_categories,
    palette="viridis",
    alpha=0.7, s=30, edgecolor=None,
    ax=axes[1, 0]
)
axes[1, 0].set_title("SARS-CoV-2 severe severity level (according to WHO)\n[Filtered for Diabetes Status = Yes]", fontsize=13, weight='bold')
axes[1, 0].set_xlabel("PC1 (37.62% Variance Explained)", fontsize=11)
axes[1, 0].set_ylabel("PC2 (10.38% Variance Explained)", fontsize=11)
axes[1, 0].legend().remove() # Unified by the panel to its right

# Panel D: Patients WITHOUT Diabetes (Diabetes Status = No)
df_diab_no = clean_df[clean_df[diabetes_col].str.lower() == "no"]
sns.scatterplot(
    data=df_diab_no,
    x="PC1", y="PC2",
    hue=who_col,
    hue_order=who_categories,
    palette="viridis",
    alpha=0.6, s=25, edgecolor=None,
    ax=axes[1, 1]
)
axes[1, 1].set_title("SARS-CoV-2 severe severity level (according to WHO)\n[Filtered for Diabetes Status = No]", fontsize=13, weight='bold')
axes[1, 1].set_xlabel("PC1 (37.62% Variance Explained)", fontsize=11)
axes[1, 1].set_ylabel("PC2 (10.38% Variance Explained)", fontsize=11)
axes[1, 1].legend(title="WHO Severity Level", bbox_to_anchor=(1.02, 1), loc='upper left', markerscale=1.5)


# Optimize layout spacing and save
plt.tight_layout()
output_image = "pca/genomic_comorbidity_analysis_4panel.png"
plt.savefig(output_image, bbox_inches='tight')

print(f"Success! 4-panel integrated figure generated.")
print(f"File saved to location: {output_image}")
