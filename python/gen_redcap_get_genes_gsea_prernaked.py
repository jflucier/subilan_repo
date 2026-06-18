import pandas as pd
import numpy as np
import gseapy as gp
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import mygene
import math
import os

print("Preparing global transcript ranking from PCA model weights (Targeting PC2)...")

# 1. Reconstruct loading weights for ALL post-variance genes from your matrix
fpkm_df = pd.read_csv("fpkm_matrix_cleaned.tsv", sep="\t", index_col=0)
pca_input = fpkm_df.T
pca_input = pca_input.loc[:, pca_input.var() > 0.1]
pca_input_log = np.log2(pca_input + 1)

scaler = StandardScaler()
scaled_data = scaler.fit_transform(pca_input_log)

pca = PCA(n_components=2)
pca.fit(scaled_data)

# --- FIX: Explicitly grab row index 1 to isolate PC2 (Infection Severity) Weights ---
pc2_weights = pca.components_[1]
rank_df = pd.DataFrame(pc2_weights, index=pca_input.columns, columns=["Weight"])

# 2. Annotate ALL genes to fetch maximum possible human symbols
mg = mygene.MyGeneInfo()
print(f"Batch-annotating all {len(rank_df)} genes for total pathway recall...")
query_results = mg.querymany(rank_df.index.tolist(), scopes="ensembl.gene", fields="symbol", species="human",
                             verbose=False)

symbol_map = {res.get("query"): res.get("symbol") for res in query_results if res.get("query") and "symbol" in res}
rank_df["Symbol"] = rank_df.index.map(symbol_map)

# Remove uncharacterized markers and sort continuously from top positive to lowest negative weight
rank_df = rank_df.dropna(subset=["Symbol"])
# Deduplicate symbols if any duplicates exist by keeping the one with highest absolute loading weight
rank_df["Abs_Weight"] = rank_df["Weight"].abs()
rank_df = rank_df.sort_values("Abs_Weight", ascending=False).drop_duplicates(subset=["Symbol"])
rank_df = rank_df.sort_values(by="Weight", ascending=False)

# Structure exactly as a 2-column DataFrame expected by gseapy (Symbol, Score)
rnk_series = rank_df[["Symbol", "Weight"]]

print(f"Submitting entire continuous gradient ({len(rnk_series)} ranked genes) to Preranked GSEA engine...")

# 3. Execute the Preranked Pathway Analysis
try:
    pre_res = gp.prerank(
        rnk=rnk_series,
        gene_sets=['Reactome_2022', 'KEGG_2021_Human'],
        processes=12,  # --- FIX: Changed from 'threads' to 'processes' ---
        min_size=5,
        max_size=500,
        permutation_num=1000,  # Runs 1000 background permutations for statistical validation
        outdir='pca/gsea_preranked_results',
        seed=42
    )

    # 4. Filter and display top results
    gsea_results = pre_res.res2d

    print("\n" + "=" * 115)
    print("🏆 TOP BIOLOGICAL PATHWAYS DRIVING YOUR WHO SEVERITY GRADIENT (PRERANKED GSEA) 🏆")
    print("=" * 115)
    if not gsea_results.empty:
        # Sort by absolute Normalized Enrichment Score to find strongest pathways on both ends
        gsea_results["Abs_NES"] = gsea_results["NES"].abs()
        gsea_results = gsea_results.sort_values(by="Abs_NES", ascending=False)

        # In GSEA Preranked, an FDR (q-val) <= 0.25 is the standard discovery threshold
        sig_gsea = gsea_results[gsea_results["FDR q-val"] <= 0.25]
        cols_to_print = ['Term', 'ES', 'NES', 'NOM p-val', 'FDR q-val']

        if not sig_gsea.empty:
            print("✅ Statistically Significant Pathways found (FDR <= 0.25):")
            print(sig_gsea[cols_to_print].head(15).to_string(index=False))
        else:
            print("ℹ️ No pathways hit the broad FDR <= 0.25 threshold. Displaying top trending pathways:")
            print(gsea_results[cols_to_print].head(15).to_string(index=False))
    else:
        print("No pathways returned from the calculation matrix.")
    print("-" * 115)

except Exception as e:
    print(f"❌ Preranked GSEA tracking run failed: {str(e)}")
