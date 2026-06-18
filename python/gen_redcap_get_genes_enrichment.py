import pandas as pd
import gseapy as gp

# Load the annotated drivers we just saved
df = pd.read_csv("pca/PC2_top_driving_genes_annotated.tsv", sep="\t")

# Filter out uncharacterized genes to get clean symbols for pathway mapping
gene_symbols = df[df["Gene_Symbol"] != "Uncharacterized"]["Gene_Symbol"].dropna().tolist()

print(f"Submitting {len(gene_symbols)} characterized genes to Enrichr...")

# Run enrichment against Reactome and KEGG pathways
enr = gp.enrichr(
    gene_list=gene_symbols,
    gene_sets=['Reactome_2022', 'KEGG_2021_Human'],
    organism='human',
    outdir='pca/enrichr_results',
    cutoff=0.05
)

# View top 5 enriched pathways
results = enr.results
if not results.empty:
    print("\n🏆 TOP ENRICHED PATHWAYS FOR PC2 (ACUTE SEVERITY):")
    print(results[['Gene_Set', 'Term', 'Adjusted P-value', 'Genes']].head(5).to_string(index=False))
else:
    print("No statistically significant pathways found with this subset size.")
