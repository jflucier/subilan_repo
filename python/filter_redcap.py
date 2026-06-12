import pandas as pd

# 1. Load the target clinical features
with open("redcap_tokeep.txt", "r", encoding="utf-8") as f:
    columns_to_keep = [line.strip() for line in f if line.strip()]

# SWITCHED TO THE LABEL CSV WHERE THE COLUMNS MATCH YOUR TEXT FILE BETTER
csv_filename = "/jbod2/def-gimap5/20250923_BQC19/data/clinical/redcap_clinical_data_label_2024-04-11.csv"
header = pd.read_csv(csv_filename, nrows=0).columns.tolist()


# 2. Build a normalized mapping dictionary
def normalize(text):
    if not isinstance(text, str): return ""
    text = text.lower().strip().replace(":", "").replace("?", "").replace(" ", "").replace("_", "")
    # Standardize broken encoding variants for Age and Fever
    text = text.replace("√çge", "age").replace("âg", "age").replace("âge", "age").replace("age", "age")
    text = text.replace("‚â•38.0", "380").replace("≥38.0", "380").replace(">=38.0", "380")
    return text


raw_mapping = {normalize(col): col for col in header}

# 3. Match columns dynamically
# 3. Match columns dynamically
valid_columns = []
missed_columns = []

for col in columns_to_keep:
    norm_col = normalize(col)

    # Check 1: Direct match
    if col in header:
        valid_columns.append(col)
    # Check 2: Fuzzy normalized match
    elif norm_col in raw_mapping:
        valid_columns.append(raw_mapping[norm_col])
    # Check 3: Targeted suffix/contains match for the Age column
    elif "au recrutement" in norm_col:
        # Pull the exact string header from the file that contains "recrutement"
        age_match = [h for h in header if "recrutement" in h.lower()]
        if age_match:
            valid_columns.append(age_match[0])  # Append the string, not a list
        else:
            missed_columns.append(col)
    # Check 4: Explicit overrides for the Fever column variations
    elif "fever" in norm_col and "38.0" in norm_col:
        matches = [h for h in header if "Fever" in h and "38.0" in h]
        for m in matches:
            if m not in valid_columns:
                valid_columns.append(m)
    else:
        missed_columns.append(col)

# Ensure BQC ID is at the top
id_col = "BQC ID" if "BQC ID" in header else ("BQCID" if "BQCID" in header else header[0])
if id_col not in valid_columns:
    valid_columns.insert(0, id_col)

valid_columns = list(dict.fromkeys(valid_columns))

print(f"Successfully matched {len(valid_columns)} columns.")

print(f"\n--- MISSED COLUMNS ({len(missed_columns)}) ---")
if missed_columns:
    for col in missed_columns:
        print(f"  - {col}")
else:
    print("  None! Perfect 100% column match achieved.")
print("---------------------------------\n")

# 4. Stream and filter using chunks
chunk_size = 10000
output_filename = "filtered_output.csv"

first_chunk = True
for chunk in pd.read_csv(csv_filename, usecols=valid_columns, chunksize=chunk_size, low_memory=False):
    chunk = chunk[valid_columns]

    if id_col in chunk.columns:
        chunk = chunk.rename(columns={id_col: "BQC ID"})

    chunk.to_csv(output_filename, mode='a', index=False, header=first_chunk)
    first_chunk = False

print(f"Done! Filtered file saved as {output_filename}")
