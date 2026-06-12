import pandas as pd

# Load the list of columns to keep and remove trailing whitespaces/newlines
with open("redcap_tokeep.txt", "r", encoding="utf-8") as f:
    columns_to_keep = [line.strip() for line in f if line.strip()]

# Read only the header of your big CSV to verify which columns actually exist
# Replace 'your_huge_file.csv' with your actual filename
csv_filename = "/jbod2/def-gimap5/20250923_BQC19/data/clinical/redcap_clinical_data_raw_2024-04-11.csv"
header = pd.read_csv(csv_filename, nrows=0).columns.tolist()

# Match columns (handles cases where text file has columns missing in the CSV)
valid_columns = [col for col in columns_to_keep if col in header]

# Read and filter the CSV using chunks to save memory
chunk_size = 10000
output_filename = "filtered_output.csv"

first_chunk = True
for chunk in pd.read_csv(csv_filename, usecols=valid_columns, chunksize=chunk_size):
    # Ensure correct column ordering as specified in your text file
    chunk = chunk[valid_columns]

    # Write to new CSV
    chunk.to_csv(output_filename, mode='a', index=False, header=first_chunk)
    first_chunk = False

print(f"Done! Filtered file saved as {output_filename}")
