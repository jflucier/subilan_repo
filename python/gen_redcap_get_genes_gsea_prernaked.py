# 4. Compile and print top comprehensive pathway trends cleanly
print("\n" + "=" * 115)
print("🏆 TOP BIOLOGICAL PATHWAYS DRIVING YOUR WHO SEVERITY GRADIENT (PRERANKED GSEA) 🏆")
print("=" * 115)

if compiled_results:
    gsea_results = pd.concat(compiled_results, ignore_index=True)

    # --- MATCH LOWERCASE COLUMNS RETURNED BY YOUR ENVIRONMENT ---
    if "nes" in gsea_results.columns:
        gsea_results["Abs_NES"] = gsea_results["nes"].abs()
        gsea_results = gsea_results.sort_values(by="Abs_NES", ascending=False)

        # Save the complete integrated file to disk cleanly
        gsea_results.drop(columns=["Abs_NES"]).to_csv("pca/gsea_preranked_compiled_results.tsv", sep="\t", index=False)

        # Standard biological discovery cutoff threshold for Preranked GSEA is False Discovery Rate <= 0.25
        # Lowercase 'fdr' used here to match environment profile
        sig_gsea = gsea_results[gsea_results["fdr"] <= 0.25]
        cols_to_print = ['Library', 'Term', 'nes', 'pval', 'fdr']

        if not sig_gsea.empty:
            print("✅ Statistically Significant Pathways found (FDR <= 0.25):")
            print(sig_gsea[cols_to_print].head(15).to_string(index=False))
        else:
            print("ℹ️ No pathways hit the broad FDR <= 0.25 threshold. Displaying top trending pathways:")
            print(gsea_results[cols_to_print].head(15).to_string(index=False))
    else:
        print("⚠️ Warning: The GSEA run completed, but structural tables lack 'nes' column.")
        print(f"Available columns: {gsea_results.columns.tolist()}")
else:
    print("❌ No comprehensive pathways returned from the calculation matrix loops.")
print("-" * 115)
