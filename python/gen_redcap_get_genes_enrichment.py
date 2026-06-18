import pandas as pd
import requests
import json
import time


def run_raw_enrichr(pc_name):
    tsv_path = f"pca/{pc_name}_method1_passing_genes.tsv"
    print(f"\nParsing drivers for {pc_name} from '{tsv_path}'...")

    # 1. Load your thresholded file
    df = pd.read_csv(tsv_path, sep="\t")
    gene_list = df[df["Gene_Symbol"] != "Uncharacterized"]["Gene_Symbol"].dropna().tolist()

    print(f"🧬 Submitting {len(gene_list)} characterized genes directly to Enrichr REST API...")
    if len(gene_list) < 5:
        print(f"⚠️ Skipping {pc_name}: Too few gene symbols for enrichment analysis.")
        return

    # Join gene symbols by a newline character for the API payload
    genes_str = "\n".join(gene_list)
    description = f"{pc_name} Thresholded Drivers"

    # 2. Add the list of genes to Enrichr
    ENRICHR_URL = 'https://maayanlab.cloud'
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        print("❌ Failed to upload gene list to Enrichr.")
        return

    user_list_id = response.json()['userListId']

    # 3. Query the specific gene sets (Reactome and KEGG)
    gene_sets = ['Reactome_2022', 'KEGG_2021_Human']
    QUERY_URL = 'https://maayanlab.cloud'

    all_results = []

    for gset in gene_sets:
        print(f"   Querying pathway library: {gset}...")
        params = {'userListId': user_list_id, 'backgroundType': gset}
        resp = requests.get(QUERY_URL, params=params)

        if not resp.ok:
            print(f"   ⚠️ Query failed for library {gset}")
            continue

        data = resp.json()[gset]

        # Enrichr returns results as a list of lists:
        # [rank, term_name, p-value, z-score, combined_score, overlapping_genes, adjusted_p-value]
        for item in data:
            all_results.append({
                'Gene_Set': gset,
                'Term': item[1],
                'P-value': item[2],
                'Adjusted_P_value': item[6],
                'Combined_Score': item[4],
                'Genes': ", ".join(item[5])
            })

    # 4. Convert to DataFrame and rank
    if not all_results:
        print("❌ No pathway enrichment results recovered.")
        return

    results_df = pd.DataFrame(all_results)
    results_df = results_df.sort_values(by="Adjusted_P_value", ascending=True)

    # Save complete table
    output_path = f"pca/{pc_name}_enrichment_api_results.tsv"
    results_df.to_csv(output_path, sep="\t", index=False)

    # 5. Print a beautiful ranked output summary
    print("\n" + "=" * 115)
    print(f"🏆 TOP ENRICHED PATHWAYS FOR {pc_name} (ACUTE SEVERITY) VIA DIRECT API 🏆")
    print("=" * 115)
    print(f"{'Gene Set':<15} | {'Term':<55} | {'Adj P-value':<12} | {'Genes Overlap'}")
    print("-" * 115)

    for idx, row in results_df.head(10).iterrows():
        print(
            f"{row['Gene_Set']:<15} | {row['Term'][:55]:<55} | {row['Adjusted_P_value']:.2e} | {row['Genes'][:30]}...")
    print("-" * 115)
    print(f"Full results successfully saved to: {output_path}")


# Execute the direct API pipeline over your post-variance lists
run_raw_enrichr("PC2")
