import pandas as pd
import numpy as np
import os
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

print("Loading raw counts and clinical metadata...")
counts_df = pd.read_csv("star/raw_counts_matrix.tsv", sep="\t", index_col=0)
# --- ADD THIS LINE TO COLLAPSE DUPLICATES ---
counts_df = counts_df.groupby(level=0).sum()
# --------------------------------------------
clinical_df = pd.read_csv("redcap.filtered_output.csv")

# Clean matrix column headers: 'star/alignment/BQC19719.VAR15535-23_trimmed_clean' -> 'BQC19719.VAR15535-23_trimmed_clean'
counts_df.columns = [os.path.basename(col) for col in counts_df.columns]

# --- 1. Filter out Bad Samples ---
print("Filtering out bad samples...")
with open("star/bad_samples.txt", "r") as f:
    bad_samples = [line.strip() for line in f if line.strip()]

# Drop the columns present in bad_samples
counts_df = counts_df.drop(columns=bad_samples, errors='ignore')
print(f"Remaining samples in matrix after filtering: {counts_df.shape[1]}")

# --- 2. Compute Gene Lengths from GTF for FPKM ---
print("Extracting gene lengths from GTF file...")
gtf_path = "/fast2/jflucier/star/Homo_sapiens.GRCh38.p14.Ensembl115.gtf"

gene_lengths = {}
with open(gtf_path, "r") as f:
    for line in f:
        if line.startswith("#"): continue
        parts = line.split("\t")
        if parts[2] == "exon":
            attrs = parts[8]
            # Fast parse of gene_id from attributes string
            if 'gene_id "' in attrs:
                gene_id = attrs.split('gene_id "')[1].split('"')[0]
                exon_len = int(parts[4]) - int(parts[3]) + 1
                gene_lengths[gene_id] = gene_lengths.get(gene_id, 0) + exon_len

lengths_series = pd.Series(gene_lengths).reindex(counts_df.index).dropna()
counts_df = counts_df.loc[lengths_series.index]

# --- 3. Calculate FPKM Matrix ---
print("Calculating FPKM...")
total_reads_m = counts_df.sum() / 1e6
gene_lengths_kb = lengths_series / 1e3

fpkm_df = counts_df.div(gene_lengths_kb, axis=0).div(total_reads_m, axis=1)
fpkm_df.to_csv("fpkm_matrix_cleaned.tsv", sep="\t")
print("Cleaned FPKM Matrix saved.")

# --- 4. Run PCA ---
print("Preparing data for PCA...")
pca_input = fpkm_df.T
pca_input = pca_input.loc[:, pca_input.var() > 0.1] # Drop low-variance genes

pca_input_log = np.log2(pca_input + 1)
scaler = StandardScaler()
scaled_data = scaler.fit_transform(pca_input_log)

print("Running PCA modeling...")
pca = PCA(n_components=5)
pca_coords = pca.fit_transform(scaled_data)

# Create PCA coordinates dataframe
pca_df = pd.DataFrame(
    pca_coords,
    columns=[f"PC{i}" for i in range(1, 6)],
    index=pca_input.index
).reset_index().rename(columns={"index": "Sample_Full_ID"})

# Extract 'BQC ID' from 'Sample_Full_ID' (e.g., 'BQC19719.VAR15535-23_trimmed_clean' -> 'BQC19719')
pca_df["BQC ID"] = pca_df["Sample_Full_ID"].str.split('.').str[0]

# --- 5. Merge with Filtered Metadata Columns ---
print("Merging PCA coordinates with clinical metadata...")
pca_df["BQC ID"] = pca_df["BQC ID"].astype(str)
clinical_df["BQC ID"] = clinical_df["BQC ID"].astype(str)

# inner join keeps samples present in both the expression profiling and the clinical file
final_merged_pca = pd.merge(pca_df, clinical_df, on="BQC ID", how="inner")
final_merged_pca.to_csv("pca/pca_with_metadata.csv", index=False)

print(f"Success! Merged {len(final_merged_pca)} total rows into 'pca_with_metadata.csv'.")
print(f"Explained Variance Ratio per PC: {pca.explained_variance_ratio_}")
