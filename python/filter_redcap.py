import pandas as pd
import os

# 1. Load the target clinical features
with open("redcap_tokeep.txt", "r", encoding="utf-8") as f:
    columns_to_keep = [line.strip() for line in f if line.strip()]

csv_filename = "/jbod2/def-gimap5/20250923_BQC19/data/clinical/redcap_clinical_data_label_2024-04-11.csv"
header = pd.read_csv(csv_filename, nrows=0).columns.tolist()


# 2. Build a normalized mapping dictionary
def normalize(text):
    if not isinstance(text, str): return ""
    text = text.lower().strip().replace(":", "").replace("?", "").replace(" ", "").replace("_", "")
    text = text.replace("âge", "age").replace("âg", "age").replace("age", "age")
    text = text.replace("≥38.0", "380").replace(">=38.0", "380")
    return text


raw_mapping = {normalize(col): col for col in header}

# 3. Match columns dynamically
valid_columns = []
for col in columns_to_keep:
    norm_col = normalize(col)
    if col in header:
        valid_columns.append(col)
    elif norm_col in raw_mapping:
        valid_columns.append(raw_mapping[norm_col])

# Force ensure BQC ID variant is at the front
id_col = "BQC ID" if "BQC ID" in header else ("BQCID" if "BQCID" in header else header)
if id_col not in valid_columns:
    valid_columns.insert(0, id_col)

valid_columns = list(dict.fromkeys(valid_columns))
print(f"Successfully matched {len(valid_columns)} columns.")

# 4. Stream and filter using chunks into a temporary unflattened file
chunk_size = 20000
temp_filename = "filtered_output_raw.tmp"
final_filename = "filtered_output.csv"

# Remove any old runs hanging around to prevent cross-contamination
for f in [temp_filename, final_filename]:
    if os.path.exists(f): os.remove(f)

print("Streaming and extracting matched columns from raw clinical file...")
first_chunk = True
for chunk in pd.read_csv(csv_filename, usecols=valid_columns, chunksize=chunk_size, low_memory=False):
    chunk = chunk[valid_columns]

    # Standardize the ID header name to match downstream expectations
    if id_col in chunk.columns:
        chunk = chunk.rename(columns={id_col: "BQC ID"})

    chunk.to_csv(temp_filename, mode='a', index=False, header=first_chunk)
    first_chunk = False

# 5. Load extracted data, flatten longitudinal rows, and save
print("Flattening repeated longitudinal records to one row per patient...")
raw_extracted_df = pd.read_csv(temp_filename, low_memory=False)

# Group by 'BQC ID' and grab the first populated (non-null) entry for each column
flattened_df = raw_extracted_df.groupby("BQC ID", as_index=False).first()

# Save clean flattened matrix
flattened_df.to_csv(final_filename, index=False)

# Clean up temp file
if os.path.exists(temp_filename):
    os.remove(temp_filename)

print(f"Done! Reduced from {len(raw_extracted_df)} records down to {len(flattened_df)} unique patient rows.")
print(f"Final clean file saved as {final_filename}")
