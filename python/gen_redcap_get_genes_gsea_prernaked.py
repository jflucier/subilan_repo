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

# Isolate PC2 (Infection Severity) Weights
pc2_weights = pca.components_[1]  # Index 1 is exactly PC2
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
rank_df["Abs_Weight"] = rank_df["Weight"].abs()
rank_df = rank_df.sort_values("Abs_Weight", ascending=False).drop_duplicates(subset=["Symbol"])
rank_df = rank_df.sort_values(by="Weight", ascending=False)

# Structure exactly as a 2-column DataFrame expected by gseapy (Symbol, Score)
rnk_series = rank_df[["Symbol", "Weight"]]

print(f"Submitting entire continuous gradient ({len(rnk_series)} ranked genes) to Preranked GSEA engine...")

# 3. Loop through individual gene sets sequentially to prevent parsing crashes
gene_libraries = ['Reactome_2022', 'KEGG_2021_Human']
compiled_results = []

for gset in gene_libraries:
    print(f"\n🚀 Running 1,000 permutations over pathway library: {gset}...")
    try:
        out_directory = f'pca/gsea_preranked_results/{gset}'
        pre_res = gp.prerank(
            rnk=rnk_series,
            gene_sets=gset,  # Pass exactly one string library name per run iteration
            processes=4,
            min_size=5,
            max_size=500,
            permutation_num=1000,
            outdir=out_directory,
            seed=42
        )

        # Capture internal result data frame layer
        res_df = pre_res.res2d
        if not res_df.empty:
            res_df["Library"] = gset
            # Force columns to lowercase to standardize across various gseapy versions
            res_df.columns = [c.lower() for c in res_df.columns]
            compiled_results.append(res_df)

    except Exception as e:
        print(f"   ❌ Preranked execution failed for library {gset}: {str(e)}")

# 4. Compile and print top comprehensive pathway trends cleanly
print("\n" + "=" * 115)
print("🏆 TOP BIOLOGICAL PATHWAYS DRIVING YOUR WHO SEVERITY GRADIENT (PRERANKED GSEA) 🏆")
print("=" * 115)

if compiled_results:
    gsea_results = pd.concat(compiled_results, ignore_index=True)

    # Save a raw copy to disk immediately so your data is secure
    gsea_results.to_csv("pca/gsea_preranked_compiled_results.tsv", sep="\t", index=False)

    # Dynamically find case-insensitive columns
    cols = {c.lower(): c for c in gsea_results.columns}

    nes_col = cols.get("nes")
    fdr_col = cols.get("fdr") or cols.get("fdr q-val") or cols.get("padj")
    pval_col = cols.get("pval") or cols.get("nom p-val")
    term_col = cols.get("term") or cols.get("name") or gsea_results.columns[0]
    lib_col = cols.get("library") or "Library"

    if nes_col and fdr_col:
        # Sort by absolute NES to show strongest results
        gsea_results["abs_nes"] = gsea_results[nes_col].abs()
        gsea_results = gsea_results.sort_values(by="abs_nes", ascending=False)

        # Standard GSEA significance filter (FDR <= 0.25)
        sig_gsea = gsea_results[gsea_results[fdr_col] <= 0.25]

        # Build strict dynamic list of present columns to print
        cols_to_print = [lib_col, term_col, nes_col, pval_col, fdr_col]
        cols_to_print = [c for c in cols_to_print if c in gsea_results.columns]

        if not sig_gsea.empty:
            print(f"✅ Statistically Significant Pathways found (FDR <= 0.25):")
            print("-" * 115)
            print(sig_gsea[cols_to_print].head(15).to_string(index=False))
        else:
            print("ℹ️ No pathways hit the broad FDR <= 0.25 threshold. Displaying top trending pathways:")
            print("-" * 115)
            print(gsea_results[cols_to_print].head(15).to_string(index=False))
    else:
        print("⚠️ GSEA finished but structure differs. Printing raw dataframe sample slice:")
        print(gsea_results.head(10).to_string())
else:
    print("❌ No comprehensive pathways returned from the calculation matrix loops.")
print("-" * 115)
