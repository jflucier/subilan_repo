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


def save_suppa_file(df, output_path):
    """Helper function to save a DataFrame in the specific SUPPA format."""
    print(f"Saving SUPPA-compatible data to: {output_path}")
    header = df.columns.tolist()
    # Header line with ONLY sample names
    header_line = '\t'.join(header) + '\n'

    with open(output_path, 'w') as f:
        f.write(header_line)
        # Data rows with index but no header
        df.to_csv(f, sep='\t', header=False, index=True)


def process_tsv_file(input_filename, output_basepath):
    print(f"Reading data from: {input_filename}")

    # Read the TSV file. low_memory=False to suppress DtypeWarning.
    df_full = pd.read_csv(input_filename, sep='\t', index_col=0, low_memory=False)

    print(f"Original data shape: {df_full.shape}")

    if 'ref_group' in df_full.index:
        ref_groups = df_full.loc['ref_group']
        df_expression = df_full.drop('ref_group')
        print("Extracted 'ref_group' metadata and removed row.")
    else:
        print("Error: 'ref_group' row not found in input file.")
        sys.exit(1)

    for col in df_expression.columns:
        df_expression[col] = pd.to_numeric(df_expression[col], errors='coerce')

    df_expression = df_expression.dropna(how='all')

    if df_expression.select_dtypes(include=[np.number]).empty:
        print("Error: No numeric columns found to transform after cleaning.")
        sys.exit(1)

    print(f"Applying inverse log2 transformation to data...")
    df_transformed = df_expression.apply(reverse_log2_tpm, axis=0)

    print("Removing version numbers (.X) from isoform IDs in the index...")
    # --- CORRECTED LINE ---
    df_transformed.index = df_transformed.index.map(lambda x: x.split('.')[0])
    # ----------------------

    # --- Split data based on ref_groups ---
    samples_low = ref_groups[ref_groups == 'LOW'].index.tolist()
    samples_high = ref_groups[ref_groups == 'HIGH'].index.tolist()

    df_low = df_transformed[samples_low]
    df_high = df_transformed[samples_high]

    # Define output paths
    output_low_path = output_basepath + '_LOW.tsv'
    output_high_path = output_basepath + '_HIGH.tsv'

    # Save the files using the helper function
    save_suppa_file(df_low, output_low_path)
    save_suppa_file(df_high, output_high_path)

    print("Transformation and splitting complete.")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <input_file_path> <output_base_filename>")
        print("Example: python rsem_to_tpm.py input.tsv output/results/tcga_ACC")
        sys.exit(1)

    # Corrected variable assignment from sys.argv list to string paths
    input_file = sys.argv[1]
    output_base = sys.argv[2]

    if os.path.exists(input_file):
        # Ensure the output directory exists before writing files
        output_dir = os.path.dirname(output_base)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            print(f"Created directory: {output_dir}")

        process_tsv_file(input_file, output_base)
    else:
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
