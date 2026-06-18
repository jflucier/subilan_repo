import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load the integrated dataset
df = pd.read_csv("pca/pca_with_metadata.csv")

who_col = "If a screening test for SARS-CoV-2 by PCR was performed, what is the most severe severity level (according to WHO) achieved?"
diabetes_col = "Diabetes ?"

# Standardize text configurations
for c in [who_col, diabetes_col]:
    if c in df.columns:
        df[c] = df[c].astype(str).str.strip()

# Drop rows missing critical metadata fields
clean_df = df.dropna(subset=[diabetes_col, who_col])
clean_df = clean_df[~clean_df[who_col].str.lower().isin(["nan", "unknown", "missing"])]
clean_df = clean_df[~clean_df[diabetes_col].str.lower().isin(["nan", "unknown", "missing"])]

# --- THREE-BIN SEVERITY FUNCTION ---
def bin_severity(val):
    v = str(val).strip()
    if v in ["Mild", "Uninfected"]:
        return "Mild / Uninfected"
    elif v == "Moderate":
        return "Moderate"
    elif v in ["Severe", "Dead"]:
        return "Severe / Dead"
    return "Mild / Uninfected"  # Fallback for unexpected cases

# Apply the 3-bin grouping transformation
binned_col = "Binned WHO Severity"
clean_df[binned_col] = clean_df[who_col].apply(bin_severity)

# Set up a 2x2 grid layout canvas
fig, axes = plt.subplots(2, 2, figsize=(24, 18), dpi=300)

# Custom color profiles for consistency across subplots (Blue -> Yellow/Orange -> Red gradient)
hue_order = ["Mild / Uninfected", "Moderate", "Severe / Dead"]
palette_colors = {"Mild / Uninfected": "#1f77b4", "Moderate": "#ff7f0e", "Severe / Dead": "#d62728"}

# =====================================================================
# ROW 1: THE OVERALL COHORT PLOTS
# =====================================================================

# Panel A: Baseline Comorbidity (Diabetes)
diab_counts = clean_df[diabetes_col].value_counts()
diab_labels = {k: f"{k} (n={v})" for k, v in diab_counts.items()}
df_panel_a = clean_df.copy()
df_panel_a[diabetes_col] = df_panel_a[diabetes_col].map(diab_labels)

sns.scatterplot(
    data=df_panel_a,
    x="PC1", y="PC2",
    hue=diabetes_col,
    hue_order=[diab_labels.get(k, k) for k in ["Yes", "No"] if k in diab_labels],
    palette="Set1",
    alpha=0.6, s=25, edgecolor=None,
    ax=axes[0, 0]
)
axes[0, 0].set_title("PCA Split by Baseline Comorbidity\n(Diabetes ?)", fontsize=13, weight='bold')
axes[0, 0].set_xlabel("PC1 (37.62% Var) - Tracks Chronic Health Status", fontsize=11)
axes[0, 0].set_ylabel("PC2 (10.38% Var)", fontsize=11)
axes[0, 0].legend(title="Diabetes Status", loc='upper right')

# Panel B: Overall Acute Presentation (3-Binned Severity)
overall_counts = clean_df[binned_col].value_counts()
overall_labels = {k: f"{k} (n={v})" for k, v in overall_counts.items()}
df_panel_b = clean_df.copy()
df_panel_b[binned_col] = df_panel_b[binned_col].map(overall_labels)

sns.scatterplot(
    data=df_panel_b,
    x="PC1", y="PC2",
    hue=binned_col,
    hue_order=[overall_labels[k] for k in hue_order if k in overall_labels],
    palette={overall_labels[k]: palette_colors[k] for k in hue_order if k in overall_labels},
    alpha=0.6, s=25, edgecolor=None,
    ax=axes[0, 1]
)
axes[0, 1].set_title("SARS-CoV-2 severe severity level (according to WHO)\n[Overall Cohort Binned]", fontsize=13, weight='bold')
axes[0, 1].set_xlabel("PC1 (37.62% Var)", fontsize=11)
axes[0, 1].set_ylabel("PC2 (10.38% Var) - Tracks Infection Severity", fontsize=11)
axes[0, 1].legend(title="Clinical State", bbox_to_anchor=(1.02, 1), loc='upper left')


# =====================================================================
# ROW 2: THE STRATIFIED COHORT PLOTS
# =====================================================================

# Panel C: Patients WITH Diabetes (Diabetes Status = Yes)
df_diab_yes = clean_df[clean_df[diabetes_col].str.lower() == "yes"].copy()
yes_counts = df_diab_yes[binned_col].value_counts()
yes_labels = {k: f"{k} (n={v})" for k, v in yes_counts.items()}
df_diab_yes[binned_col] = df_diab_yes[binned_col].map(yes_labels)

sns.scatterplot(
    data=df_diab_yes,
    x="PC1", y="PC2",
    hue=binned_col,
    hue_order=[yes_labels[k] for k in hue_order if k in yes_labels],
    palette={yes_labels[k]: palette_colors[k] for k in hue_order if k in yes_labels},
    alpha=0.7, s=30, edgecolor=None,
    ax=axes[1, 0]
)
axes[1, 0].set_title("SARS-CoV-2 severe severity level (according to WHO)\n[Filtered for Diabetes Status = Yes]", fontsize=13, weight='bold')
axes[1, 0].set_xlabel("PC1 (37.62% Variance Explained)", fontsize=11)
axes[1, 0].set_ylabel("PC2 (10.38% Variance Explained)", fontsize=11)
axes[1, 0].legend(title="Clinical State (Diabetes=Yes)", loc='upper right')

# Panel D: Patients WITHOUT Diabetes (Diabetes Status = No)
df_diab_no = clean_df[clean_df[diabetes_col].str.lower() == "no"].copy()
no_counts = df_diab_no[binned_col].value_counts()
no_labels = {k: f"{k} (n={v})" for k, v in no_counts.items()}
df_diab_no[binned_col] = df_diab_no[binned_col].map(no_labels)

sns.scatterplot(
    data=df_diab_no,
    x="PC1", y="PC2",
    hue=binned_col,
    hue_order=[no_labels[k] for k in hue_order if k in no_labels],
    palette={no_labels[k]: palette_colors[k] for k in hue_order if k in no_labels},
    alpha=0.6, s=25, edgecolor=None,
    ax=axes[1, 1]
)
axes[1, 1].set_title("SARS-CoV-2 severe severity level (according to WHO)\n[Filtered for Diabetes Status = No]", fontsize=13, weight='bold')
axes[1, 1].set_xlabel("PC1 (37.62% Variance Explained)", fontsize=11)
axes[1, 1].set_ylabel("PC2 (10.38% Variance Explained)", fontsize=11)
axes[1, 1].legend(title="Clinical State (Diabetes=No)", bbox_to_anchor=(1.02, 1), loc='upper left', markerscale=1.5)


# Optimize spacing parameters and export image file
plt.tight_layout()
output_image = "pca/genomic_comorbidity_binned_3groups.png"
plt.savefig(output_image, bbox_inches='tight')

print("Success! 3-bin integrated layout figure generated.")
print(f"File saved to location: {output_image}")
