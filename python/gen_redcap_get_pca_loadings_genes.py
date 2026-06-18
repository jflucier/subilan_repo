import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import mygene

# 1. Reconstruct your PCA transformation steps
print("Loading expression profile matrix...")
fpkm_df = pd.read_csv("fpkm_matrix_cleaned.tsv", sep="\t", index_col=0)

# Re-apply identical transformations
pca_input = fpkm_df.T
pca_input = pca_input.loc[:, pca_input.var() > 0.1]
pca_input_log = np.log2(pca_input + 1)

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

# 3. Export top driving variables with full BioMart annotations
print("\n" + "=" * 115)
print("🧬 TOP DRIVING GENES WITH BIOMART ANNOTATIONS 🧬")
print("=" * 115)

mg = mygene.MyGeneInfo()

for pc in ["PC1", "PC2"]:
    col = f"{pc}_Loading"

    # Highest absolute weights indicate highest directional variance influence
    top_genes = loadings[col].abs().sort_values(ascending=False).head(20)
    top_signed_genes = loadings.loc[top_genes.index, [col]].copy()

    # Batch query BioMart for Symbols and Names via mygene
    ensembl_ids = top_signed_genes.index.tolist()
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
    top_signed_genes["Gene_Symbol"] = top_signed_genes.index.map(symbol_map)
    top_signed_genes["Gene_Name"] = top_signed_genes.index.map(name_map)

    print(f"\n🚀 Top 20 Most Influential Driver Genes along {pc}:")
    print("-" * 115)
    print(f"{'Ensembl ID':<20} | {'Symbol':<15} | {'Weight':<10} | {'Direction':<12} | {'Full Gene Name'}")
    print("-" * 115)

    for gene, row in top_signed_genes.iterrows():
        weight = row[col]
        direction = "📈 Positive" if weight > 0 else "📉 Negative"
        print(f"{gene:<20} | {row['Gene_Symbol']:<15} | {weight:+.4f} | {direction:<12} | {row['Gene_Name']}")

    # Save complete annotated results to a clean output text file
    top_signed_genes.to_csv(f"pca/{pc}_top_driving_genes_annotated.tsv", sep="\t")

print("\n" + "-" * 115)
print(
    "Loadings files successfully saved to 'pca/PC1_top_driving_genes_annotated.tsv' and 'pca/PC2_top_driving_genes_annotated.tsv'!")
