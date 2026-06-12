import pandas as pd

# Load the list of columns to keep and remove trailing whitespaces/newlines
with open("redcap_tokeep.txt", "r", encoding="utf-8") as f:
    columns_to_keep = [line.strip() for line in f if line.strip()]

csv_filename = "/jbod2/def-gimap5/20250923_BQC19/data/clinical/redcap_clinical_data_raw_2024-04-11.csv"
header = pd.read_csv(csv_filename, nrows=0).columns.tolist()

# Match standard columns from your text file
valid_columns = [col for col in columns_to_keep if col in header]

# --- FIX: Dynamically catch BQCID and force it to the front ---
if "BQCID" in header and "BQCID" not in valid_columns:
    valid_columns.insert(0, "BQCID")
elif "BQC ID" in header and "BQC ID" not in valid_columns:
    valid_columns.insert(0, "BQC ID")

# Read and filter the CSV using chunks
chunk_size = 10000
output_filename = "filtered_output.csv"

first_chunk = True
# Set low_memory=False to suppress the mixed-type warning
for chunk in pd.read_csv(csv_filename, usecols=valid_columns, chunksize=chunk_size, low_memory=False):
    chunk = chunk[valid_columns]

    # --- FIX: Standardize column header name to 'BQC ID' for the PCA script ---
    if "BQCID" in chunk.columns:
        chunk = chunk.rename(columns={"BQCID": "BQC ID"})

    # Write to new CSV
    chunk.to_csv(output_filename, mode='a', index=False, header=first_chunk)
    first_chunk = False

print(f"Done! Filtered file saved as {output_filename}")
