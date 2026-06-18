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

    # 2. Run over human functional pathways (Setting a relaxed cutoff of 1.0 to capture suggestive trends)
    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=['Reactome_2022', 'KEGG_2021_Human'],
            organism='human',
            outdir=f'pca/enrichment_results/{pc_name}',
            cutoff=1.0  # Captures all suggestive pathways for exploratory view
        )

        results_df = enr.results

        # 3. Print summaries safely by verifying index footprints first
        print(f"\n🏆 FUNCTIONAL PATHWAY PROFILE FOR {pc_name} (ACUTE SEVERITY):")
        print("-" * 115)

        if not results_df.empty and "Term" in results_df.columns:
            # Sort strictly by nominal P-value to pull strongest biological trends to the top
            results_df = results_df.sort_values(by="P-value", ascending=True)

            # Separate into significant vs suggestive blocks
            if "Adjusted P-value" in results_df.columns:
                sig_results = results_df[results_df["Adjusted P-value"] <= 0.05]
            else:
                sig_results = pd.DataFrame()

            if not sig_results.empty:
                print("✅ Statistically Significant Terms found (FDR <= 0.05):")
                cols_to_print = [c for c in ['Gene_Set', 'Term', 'Adjusted P-value', 'Genes'] if
                                 c in results_df.columns]
                print(sig_results[cols_to_print].head(10).to_string(index=False))
            else:
                print("ℹ️ No terms passed the strict Adjusted P-value <= 0.05 FDR threshold.")
                print("   Displaying Top 10 Suggestive Trends (Ranked by raw nominal P-value):")
                print("." * 115)
                cols_to_print = [c for c in ['Gene_Set', 'Term', 'P-value', 'Adjusted P-value', 'Genes'] if
                                 c in results_df.columns]
                print(results_df[cols_to_print].head(10).to_string(index=False))
        else:
            print("❌ No matching records returned from the enrichment database libraries.")

        print("-" * 115)

    except Exception as e:
        print(f"❌ Enrichment execution failed for {pc_name}: {str(e)}")


# Execute the analysis over your post-variance lists
run_enrichment_for_component("PC1")
run_enrichment_for_component("PC2")
