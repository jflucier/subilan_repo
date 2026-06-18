import pandas as pd
import numpy as np
import gseapy as gp
import os

print("Preparing global transcript ranking from PCA model weights...")

# 1. Reconstruct loading weights for ALL post-variance genes from your matrix
fpkm_df = pd.read_csv("fpkm_matrix_cleaned.tsv", sep="\t", index_col=0)
pca_input = fpkm_df.T
pca_input = pca_input.loc[:, pca_input.var() > 0.1]
pca_input_log = np.log2(pca_input + 1)

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()
scaled_data = scaler.fit_transform(pca_input_log)

pca = PCA(n_components=2)
pca.fit(scaled_data)

# Extract ALL genes into a single ranking frame
rank_df = pd.DataFrame(pca.components_[1], index=pca_input.columns, columns=["Weight"])

# 2. Annotate ALL genes to fetch maximum possible human symbols
import mygene

mg = mygene.MyGeneInfo()
print(f"Batch-annotating all {len(rank_df)} genes for total pathway recall...")
query_results = mg.querymany(rank_df.index.tolist(), scopes="ensembl.gene", fields="symbol", species="human",
                             verbose=False)

symbol_map = {res.get("query"): res.get("symbol") for res in query_results if res.get("query") and "symbol" in res}
rank_df["Symbol"] = rank_df.index.map(symbol_map)

# Remove uncharacterized markers and sort continuously from top positive to lowest negative weight
rank_df = rank_df.dropna(subset=["Symbol"])
rank_df = rank_df.sort_values(by="Weight", ascending=False)

# Structure exactly as expected by the gseapy preranked API engine (2 columns: Symbol, Score)
rnk_series = rank_df[["Symbol", "Weight"]]

print(f"Submitting entire continuous gradient ({len(rnk_series)} ranked genes) to Preranked GSEA engine...")

# 3. Execute the Preranked Pathway Analysis
try:
    pre_res = gp.prerank(
        rnk=rnk_series,
        gene_sets=['Reactome_2022', 'KEGG_2021_Human'],
        threads=4,
        min_size=5,
        max_size=500,
        permutation_num=1000,  # Runs 1000 background permutations for robust statistical validation
        outdir='pca/gsea_preranked_results',
        seed=42
    )

    # 4. Filter and display top results
    gsea_results = pre_res.res2d.sort_values(by="NES", key=abs, ascending=False)  # Rank by Normalized Enrichment Score

    print("\n" + "=" * 115)
    print("🏆 TOP BIOLOGICAL PATHWAYS DRIVING YOUR WHO SEVERITY GRADIENT (PRERANKED GSEA) 🏆")
    print("=" * 115)
    if not gsea_results.empty:
        # Pull out pathways that are significantly enriched on either the severe or mild ends
        sig_gsea = gsea_results[gsea_results["FDR q-val"] <= 0.25]
        if not sig_gsea.empty:
            print(sig_gsea[['Term', 'ES', 'NES', 'NOM p-val', 'FDR q-val']].head(12).to_string(index=False))
        else:
            print("No pathways hit FDR <= 0.25. Displaying top trending structural pathways:")
            print(gsea_results[['Term', 'ES', 'NES', 'NOM p-val', 'FDR q-val']].head(10).to_string(index=False))
    else:
        print("No pathways returned from the calculation matrix.")
    print("-" * 115)

except Exception as e:
    print(f"❌ Preranked GSEA tracking run failed: {str(e)}")
