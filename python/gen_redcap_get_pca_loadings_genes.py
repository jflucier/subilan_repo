import pandas as pd
import numpy as np

# 1. Reconstruct your PCA transformation steps to grab components
# We pull the inputs exactly as processed inside your filter_redcap_fpkm_matrix.py script
print("Loading expression profile matrix...")
fpkm_df = pd.read_csv("fpkm_matrix_cleaned.tsv", sep="\t", index_col=0)

# Re-apply identical transformations
pca_input = fpkm_df.T
pca_input = pca_input.loc[:, pca_input.var() > 0.1]
pca_input_log = np.log2(pca_input + 1)

# Fit PCA to extract feature components
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()
scaled_data = scaler.fit_transform(pca_input_log)

pca = PCA(n_components=2)
pca.fit(scaled_data)

# 2. Extract loading arrays (weights of individual genes per component)
loadings = pd.DataFrame(
    pca.components_.T,
    columns=["PC1_Loading", "PC2_Loading"],
    index=pca_input.columns
)

# 3. Export top positive and negative driving variables
print("\n" + "=" * 70)
print("🧬 TOP DRIVING GENES RESPONSIBLE FOR HISTOGRAM SHIFTS 🧬")
print("=" * 70)

for pc in ["PC1", "PC2"]:
    col = f"{pc}_Loading"

    # Highest absolute weights indicate highest directional variance influence
    top_genes = loadings[col].abs().sort_values(ascending=False).head(20)
    top_signed_genes = loadings.loc[top_genes.index, [col]]

    print(f"\n🚀 Top 20 Most Influential Driver Genes along {pc}:")
    print("-" * 50)
    for gene, weight in top_signed_genes[col].items():
        direction = "📈 Positive Shift" if weight > 0 else "📉 Negative Shift"
        print(f"  - {gene:<20} | Weight = {weight:+.4f} ({direction})")

    # Save coordinate results to a clean output text file
    top_signed_genes.to_csv(f"pca/{pc}_top_driving_genes.tsv", sep="\t")

print("\nLoadings files successfully saved to 'pca/PC1_top_driving_genes.tsv' and 'pca/PC2_top_driving_genes.tsv'!")
