import pandas as pd
import numpy as np
import sys
import os


def reverse_log2_tpm(log2_value):
    """
    Reverses the transformation log2(TPM + 0.001) to get the linear TPM value.
    TPM = 2^(log2_value) - 0.001
    Ensures resulting TPM values are non-negative.
    """
    pseudo_count = 0.001
    linear_tpm = np.power(2, log2_value) - pseudo_count
    linear_tpm[linear_tpm < 0] = 0
    return linear_tpm


def process_tsv_file(input_filename, output_filename):
    print(f"Reading data from: {input_filename}")

    # Adding low_memory=False to suppress the DtypeWarning you saw earlier
    df = pd.read_csv(input_filename, sep='\t', index_col=0, low_memory=False)

    print(f"Original data shape: {df.shape}")

    if 'ref_group' in df.index:
        print("Removing 'ref_group' metadata row...")
        df = df.drop('ref_group')

    print(f"Data shape after cleaning: {df.shape}")

    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    df = df.dropna(how='all')

    if df.select_dtypes(include=[np.number]).empty:
        print("Error: No numeric columns found to transform after cleaning.")
        return

    print(f"Applying inverse log2 transformation to data...")
    df = df.apply(reverse_log2_tpm, axis=0)

    print(f"Saving linear TPM data to: {output_filename}")
    df.to_csv(output_filename, sep='\t', index=True)
    print("Transformation complete.")


if __name__ == "__main__":
    # Check if exactly two extra arguments (input and output file paths) were provided
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <input_file_path> <output_file_path>")
        print("Example: python rsem_to_tpm.py input.tsv output.tsv")
        sys.exit(1)  # Exit the script if arguments are incorrect

    # The first argument is the script name itself (sys.argv[0])
    # The second argument is the input file path (sys.argv[1])
    # The third argument is the output file path (sys.argv[2])
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    if os.path.exists(input_file):
        process_tsv_file(input_file, output_file)
    else:
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
