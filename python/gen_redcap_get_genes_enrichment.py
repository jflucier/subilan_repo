import pandas as pd
import gseapy as gp
import os

# Create an output directory for the results
os.makedirs("pca/enrichment_results", exist_ok=True)


def run_enrichment_for_component(pc_name):
    tsv_path = f"pca/{pc_name}_method1_passing_genes.tsv"
    print(f"\nParsing drivers for {pc_name} from '{tsv_path}'...")

    # 1. Load your thresholded matrix file
    df = pd.read_csv(tsv_path, sep="\t")

    # Filter out uncharacterized lncRNAs to build a clean list of characterized symbols
    gene_list = df[df["Gene_Symbol"] != "Uncharacterized"]["Gene_Symbol"].dropna().tolist()

    print(f"🧬 Submitting {len(gene_list)} characterized genes to the Enrichr API engine...")
    if len(gene_list) < 5:
        print(f"⚠️ Skipping {pc_name}: Too few characterized gene symbols to run significant enrichment.")
        return

    # 2. Run over human functional pathways
    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=['Reactome_2022', 'KEGG_2021_Human'],
            organism='human',
            outdir=f'pca/enrichment_results/{pc_name}',
            cutoff=0.05
        )

        # 3. Print a ranked summary table directly to your terminal screen
        results_df = enr.results
        if not results_df.empty:
            # Filter for statistically significant entries
            sig_results = results_df[results_df["Adjusted P-value"] <= 0.05]
            print(f"🏆 TOP SIGNIFICANT PATHWAYS FOR {pc_name}:")
            print("-" * 115)
            if not sig_results.empty:
                print(sig_results[['Gene_Set', 'Term', 'Adjusted P-value', 'Genes']].head(10).to_string(index=False))
            else:
                print("No single pathway passed the strict Adjusted P-value <= 0.05 FDR threshold.")
                print("Top suggestive terms instead:")
                print(results_df[['Gene_Set', 'Term', 'P-value', 'Genes']].head(3).to_string(index=False))
            print("-" * 115)
        else:
            print("No matching records returned from the enrichment database.")

    except Exception as e:
        print(f"❌ Enrichment execution failed for {pc_name}: {str(e)}")


# Execute the analysis over your post-variance lists
run_enrichment_for_component("PC1")
run_enrichment_for_component("PC2")
