# --- REVISED MULTIPLE FOR PATHWAY ENRICHMENT ---
# Changing to 1.5x background noise yields a robust list for pathway tools
method1_threshold = 1.5 * uniform_baseline

print("\n" + "=" * 75)
print("📊 DYNAMIC MATHEMATICAL INTEGRITY AUDIT (ADJUSTED FOR PATHWAYS)")
print("=" * 75)
print(f"  - Precise Number of Post-Variance Genes (N) : {num_genes}")
print(f"  - Uniform Distribution Background Baseline   : {uniform_baseline:.5f}")
print(f"  - Method 1 Mathematical Cutoff (1.5x Noise) : {method1_threshold:.5f}")
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

    # Capture all genes exceeding the revised 1.5x background threshold
    all_passing_genes = loadings[loadings[col].abs() >= method1_threshold].copy()

    # Sort by absolute strength so highest drivers print first
    all_passing_genes["Absolute_Weight"] = all_passing_genes[col].abs()
    all_passing_genes = all_passing_genes.sort_values(by="Absolute_Weight", ascending=False)

    total_found = len(all_passing_genes)
    print(f"⚙️ Component {pc}: Found {total_found} genes passing the adjusted threshold (>= {method1_threshold:.5f})")

    if total_found == 0:
        print(f"  ⚠️ Warning: No genes passed the threshold for {pc}.")
        continue

    # Slice output display to top 25 in terminal to avoid console flooding, but save ALL to TSV
    display_genes = all_passing_genes.head(25)

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

    # 4. Print clean annotated tabular display (showing first 25 preview rows)
    print("-" * 115)
    print(f"Displaying top 25 / {total_found} total passing genes for {pc}:")
    print("-" * 115)
    print(f"{'Ensembl ID':<20} | {'Symbol':<15} | {'Weight':<10} | {'Direction':<12} | {'Full Gene Name'}")
    print("-" * 115)

    for gene, row in all_passing_genes.head(25).iterrows():
        weight = row[col]
        direction = "📈 Positive" if weight > 0 else "📉 Negative"
        print(f"{gene:<20} | {row['Gene_Symbol']:<15} | {weight:+.4f} | {direction:<12} | {row['Gene_Name']}")

    # Save all passing annotated results to a clean output text file
    output_tsv = f"pca/{pc}_method1_passing_genes.tsv"
    all_passing_genes.drop(columns=["Absolute_Weight"]).to_csv(output_tsv, sep="\t")
    print("-" * 115)
    print(f"Saved ALL {total_found} thresholded variables to: {output_tsv}\n")
