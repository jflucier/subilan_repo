import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import mygene
import math

# 1. Reconstruct your PCA transformation steps
print("Loading expression profile matrix...")
fpkm_df = pd.read_csv("fpkm_matrix_cleaned.tsv", sep="\t", index_col=0)

# Re-apply identical transformations
pca_input = fpkm_df.T
pca_input = pca_input.loc[:, pca_input.var() > 0.1]
pca_input_log = np.log2(pca_input + 1)

# --- DYNAMICALLY DETERMINE PRECISE GENE COUNT ---
num_genes = pca_input.shape[1]
uniform_baseline = 1.0 / math.sqrt(num_genes)
# Set target threshold at 2.5 times background noise to extract true signals
method1_threshold = 2.5 * uniform_baseline

print("\n" + "=" * 75)
print("📊 DYNAMIC MATHEMATICAL INTEGRITY AUDIT")
print("=" * 75)
print(f"  - Precise Number of Post-Variance Genes (N) : {num_genes}")
print(f"  - Uniform Distribution Background Baseline   : {uniform_baseline:.5f}")
print(f"  - Method 1 Mathematical Cutoff (2.5x Noise) : {method1_threshold:.5f}")
print("=" * 75 + "\n")

# Fit PCA to extract feature components
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

# 3. Process components and apply Method 1 filter
mg = mygene.MyGeneInfo()

for pc in ["PC1", "PC2"]:
    col = f"{pc}_Loading"

    # --- METHOD 1 DYNAMIC FILTER: Capture all genes exceeding the threshold ---
    all_passing_genes = loadings[loadings[col].abs() >= method1_threshold].copy()

    # Sort by absolute strength so highest drivers print first
    all_passing_genes["Absolute_Weight"] = all_passing_genes[col].abs()
    all_passing_genes = all_passing_genes.sort_values(by="Absolute_Weight", ascending=False)

    total_found = len(all_passing_genes)
    print(f"⚙️ Component {pc}: Found {total_found} genes passing the Method 1 threshold (>= {method1_threshold:.5f})")

    if total_found == 0:
        print(f"  ⚠️ Warning: No genes passed the threshold for {pc}.")
        continue

    # Batch query BioMart for symbols and descriptions via mygene API
    ensembl_ids = all_passing_genes.index.tolist()
    query_results = mg.querymany(
        ensembl_ids,
        scopes="ensembl.gene",
        fields="symbol,name",
        species="human",
        verbose=False
    )

    # Map query elements into local mappers
    symbol_map = {}
    name_map = {}
    for res in query_results:
        ens_id = res.get("query")
        if ens_id:
            symbol_map[ens_id] = res.get("symbol", "Uncharacterized")
            name_map[ens_id] = res.get("name", "Long non-coding RNA / Novel Feature")

    # Add annotations to dataframe
    all_passing_genes["Gene_Symbol"] = all_passing_genes.index.map(symbol_map)
    all_passing_genes["Gene_Name"] = all_passing_genes.index.map(name_map)

    # 4. Print clean annotated tabular display
    print("-" * 115)
    print(f"{'Ensembl ID':<20} | {'Symbol':<15} | {'Weight':<10} | {'Direction':<12} | {'Full Gene Name'}")
    print("-" * 115)

    for gene, row in all_passing_genes.iterrows():
        weight = row[col]
        direction = "📈 Positive" if weight > 0 else "📉 Negative"
        print(f"{gene:<20} | {row['Gene_Symbol']:<15} | {weight:+.4f} | {direction:<12} | {row['Gene_Name']}")

    # Save all passing annotated results to a clean output text file
    output_tsv = f"pca/{pc}_method1_passing_genes.tsv"
    # Clean up tracking column before exporting
    all_passing_genes.drop(columns=["Absolute_Weight"]).to_csv(output_tsv, sep="\t")
    print("-" * 115)
    print(f"Saved {total_found} thresholded variables to: {output_tsv}\n")
