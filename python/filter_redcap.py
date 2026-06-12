import pandas as pd

# 1. Load the target clinical features and normalize them
with open("redcap_tokeep.txt", "r", encoding="utf-8") as f:
    # Strip spaces/newlines, lower case, and drop trailing punctuation like ':' or '?'
    columns_to_keep = [line.strip() for line in f if line.strip()]

csv_filename = "/jbod2/def-gimap5/20250923_BQC19/data/clinical/redcap_clinical_data_label_2024-04-11.csv"
header = pd.read_csv(csv_filename, nrows=0).columns.tolist()


# 2. Build a fuzzy mapping dictionary (maps normalized strings to exact raw columns)
def normalize(text):
    if not isinstance(text, str): return ""
    # Lowercase, strip spaces, remove common punctuation/accents that break mappings
    text = text.lower().strip().replace(":", "").replace("?", "").replace(" ", "").replace("_", "")
    text = text.replace("√çge", "age").replace("âg", "age").replace("é", "e").replace("à", "a")
    return text


raw_mapping = {normalize(col): col for col in header}

# 3. Match columns dynamically
valid_columns = []
for col in columns_to_keep:
    norm_col = normalize(col)

    # Check 1: Direct exact match
    if col in header:
        valid_columns.append(col)
    # Check 2: Fuzzy normalized match (handles Sexe: -> female/sexe variations)
    elif norm_col in raw_mapping:
        valid_columns.append(raw_mapping[norm_col])
    # Check 3: Common structural synonyms
    elif "sexe" in norm_col and "female" in header:
        valid_columns.append("female")
    elif "age" in norm_col and "age" in header:
        valid_columns.append("age")

# Ensure our variant of BQCID is definitely added at the top
id_col = "BQCID" if "BQCID" in header else ("BQC ID" if "BQC ID" in header else header[0])
if id_col not in valid_columns:
    valid_columns.insert(0, id_col)

# Deduplicate valid columns list while keeping order
valid_columns = list(dict.fromkeys(valid_columns))

print(f"Successfully matched {len(valid_columns)} columns out of {len(columns_to_keep)} requests.")

# 4. Stream and filter using chunks
chunk_size = 10000
output_filename = "filtered_output.csv"

first_chunk = True
for chunk in pd.read_csv(csv_filename, usecols=valid_columns, chunksize=chunk_size, low_memory=False):
    chunk = chunk[valid_columns]

    # Standardize the ID header name to match your PCA script expectations
    if id_col in chunk.columns:
        chunk = chunk.rename(columns={id_col: "BQC ID"})

    chunk.to_csv(output_filename, mode='a', index=False, header=first_chunk)
    first_chunk = False

print(f"Done! Filtered file saved as {output_filename}")
