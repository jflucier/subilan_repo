import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load the fully aligned dataset
df = pd.read_csv("pca/pca_with_metadata.csv")

# Clean up any potential white-space strings in the target variables
for c in df.columns:
    if df[c].dtype == 'object':
        df[c] = df[c].astype(str).str.strip()

# Set up a side-by-side canvas
fig, axes = plt.subplots(1, 2, figsize=(20, 8), dpi=300)

# --- PLOT 1: Baseline Comorbidity (Diabetes) ---
df_diabetes = df.dropna(subset=["Diabetes ?"])
sns.scatterplot(
    data=df_diabetes,
    x="PC1", y="PC2",
    hue="Diabetes ?",
    palette="Set1",
    alpha=0.6, s=20, edgecolor=None,
    ax=axes[0]
)
axes[0].set_title("PCA Split by Baseline Comorbidity\n(Diabetes ?)", fontsize=14, weight='bold')
axes[0].set_xlabel("PC1 (37.62% Var) - Tracks Chronic Health Status", fontsize=11)
axes[0].set_ylabel("PC2 (10.38% Var)", fontsize=11)
axes[0].legend(title="Diabetes Status", loc='upper right')

# --- PLOT 2: Acute Presentation (WHO Severity Level) ---
who_col = "If a screening test for SARS-CoV-2 by PCR was performed, what is the most severe severity level (according to WHO) achieved?"
df_who = df.dropna(subset=[who_col])
# Exclude generic text strings like "nan" or "Unknown" if present
df_who = df_who[~df_who[who_col].str.lower().isin(["nan", "unknown", "missing"])]

# Sort categories logically
categories = sorted(df_who[who_col].unique())

sns.scatterplot(
    data=df_who,
    x="PC1", y="PC2",
    hue=who_col,
    hue_order=categories,
    palette="viridis",
    alpha=0.6, s=20, edgecolor=None,
    ax=axes[1]
)
# --- UPDATED SUBPLOT TITLE ---
axes[1].set_title("SARS-CoV-2 severe severity level (according to WHO)", fontsize=14, weight='bold')
axes[1].set_xlabel("PC1 (37.62% Var)", fontsize=11)
axes[1].set_ylabel("PC2 (10.38% Var) - Tracks Infection Severity", fontsize=11)
axes[1].legend(title="WHO Level", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
output_image = "pca/genomic_split_dual_plot.png"
plt.savefig(output_image, bbox_inches='tight')
print(f"Success! High-resolution comparative plot saved to: {output_image}")
