import pandas as pd
import os
import csv

# 1. Load the target clinical features
with open("redcap_tokeep.txt", "r", encoding="utf-8") as f:
    columns_to_keep = [line.strip() for line in f if line.strip()]

csv_filename = "/jbod2/def-gimap5/20250923_BQC19/data/clinical/redcap_clinical_data_label_2024-04-11.csv"

# --- FIX: Use utf-8-sig to automatically strip the hidden \ufeff BOM character ---
with open(csv_filename, "r", encoding="utf-8-sig") as f:
    reader = csv.reader(f)
    header = next(reader)
# ---------------------------------------------------------------------------------


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
id_col = "BQC ID" if "BQC ID" in header else ("BQCID" if "BQCID" in header else header[0])
if id_col not in valid_columns:
    valid_columns.insert(0, id_col)

# This will now succeed cleanly because valid_columns only contains strings
valid_columns = list(dict.fromkeys(valid_columns))
print(f"Successfully matched {len(valid_columns)} columns cleanly.")


# 4. Stream and filter using chunks into a temporary file
chunk_size = 20000
temp_filename = "filtered_output_raw.tmp"
final_filename = "filtered_output.csv"

for f in [temp_filename, final_filename]:
    if os.path.exists(f): os.remove(f)

print("Streaming and extracting matched columns natively supporting internal commas...")
first_chunk = True

# Explicitly using engine='c' with standard quoting parameters isolates comma fields
for chunk in pd.read_csv(csv_filename, usecols=valid_columns, chunksize=chunk_size, low_memory=False, engine='c', quoting=csv.QUOTE_MINIMAL):
    chunk = chunk[valid_columns]
    if id_col in chunk.columns:
        chunk = chunk.rename(columns={id_col: "BQC ID"})
    chunk.to_csv(temp_filename, mode='a', index=False, header=first_chunk)
    first_chunk = False

# 5. Load extracted data and re-evaluate conflicts safely
print("Analyzing dataset for longitudinal value conflicts...")
raw_extracted_df = pd.read_csv(temp_filename, low_memory=False)

unique_counts = raw_extracted_df.groupby("BQC ID").nunique(dropna=True)
conflicting_cols = unique_counts.columns[(unique_counts > 1).any()].tolist()

print(f"\n--- COLUMNS WITH MULTIPLE DISTINCT VALUES FOUND ({len(conflicting_cols)}) ---")
if conflicting_cols:
    for col in conflicting_cols:
        try:
            raw_index = header.index(col)
        except ValueError:
            raw_index = "Extracted"
        num_patients_conflicted = (unique_counts[col] > 1).sum()
        print(f"  - [Raw Col Index: {raw_index}] {col} ({num_patients_conflicted} patients conflicted)")
else:
    print("  None! Every column holds perfectly consistent values across rows for all patients.")
print("----------------------------------------------------------------\n")

# 6. Flatten longitudinal records with advanced clinical logic and save
print("Flattening repeated longitudinal records to one row per patient using targeted logic...")
raw_extracted_df = raw_extracted_df.copy()

max_value_columns = [
    "Did or does the patient receive ventilatory support?",
    "Venous lactate:",
    "D-Dimer:",
    "IL-6:"
]

vent_col = "Did or does the patient receive ventilatory support?"
if vent_col in raw_extracted_df.columns:
    mapping_dict = {"yes": 1, "no": 0, "1": 1, "0": 0, "true": 1, "false": 0, 1: 1, 0: 0}
    raw_extracted_df[vent_col] = raw_extracted_df[vent_col].astype(str).str.lower().str.strip().map(mapping_dict)

for col in ["Venous lactate:", "D-Dimer:", "IL-6:"]:
    if col in raw_extracted_df.columns:
        raw_extracted_df[col] = pd.to_numeric(raw_extracted_df[col], errors='coerce')

base_flattened = raw_extracted_df.groupby("BQC ID", as_index=False).first()
metrics_flattened = raw_extracted_df.groupby("BQC ID", as_index=False)[["BQC ID"] + [c for c in max_value_columns if c in raw_extracted_df.columns]].max()

flattened_df = pd.merge(
    base_flattened.drop(columns=[c for c in max_value_columns if c in base_flattened.columns]),
    metrics_flattened,
    on="BQC ID"
)

flattened_df.to_csv(final_filename, index=False)
if os.path.exists(temp_filename): os.remove(temp_filename)
print(f"Done! Cleaned file saved as {final_filename}")
