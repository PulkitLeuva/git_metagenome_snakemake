import pandas as pd
import sys

def parse_mpa4(input_file, output_file, output_file_df1, output_file_df2):
    # Read the file, skipping the initial metadata lines
    with open(input_file) as f:
        lines = f.readlines()
    
    # Extract the SampleID from the 3rd line
    sample_id = lines[3].strip().split('\t')[1]
    print("sample_id ===" + sample_id)
    
    # Load the rest of the file into a DataFrame, skipping the first 4 lines
    df = pd.read_csv(input_file, sep="\t", skiprows=4, header=None)
    
    # Set the header correctly
    df.columns = ['clade_name', 'NCBI_tax_id', 'relative_abundance', 'additional_species']
    
    # Split the 'clade_name' column into multiple columns based on '|'
    clade_split = df['clade_name'].str.split('|', expand=True)
    
    # Rename columns dynamically based on how many levels there are
    clade_split.columns = [f'level_{i+1}' for i in range(clade_split.shape[1])]
    
    # Concatenate the original dataframe with the new clade levels
    df = pd.concat([clade_split, df[['NCBI_tax_id', 'relative_abundance', 'additional_species']]], axis=1)
    
    # Check if 'level_2' exists in df and create df1
    if 'level_2' in df.columns:
        df1 = df[['level_2', 'relative_abundance']]
        df1.dropna(inplace=True)  # Remove rows with None or NaN values
        df1.drop_duplicates(inplace=True)  # Remove duplicate rows
        # Rename columns
        df1.rename(columns={'level_2': 'Phylum', 'relative_abundance': sample_id}, inplace=True)
        df1.to_excel(output_file_df1, index=False)
    else:
        print("Warning: 'level_2' column not found in the DataFrame.")
    
    # Check if 'level_7' exists in df and create df2
    if 'level_7' in df.columns:
        df2 = df[['level_7', 'relative_abundance']]
        df2.dropna(inplace=True)  # Remove rows with None or NaN values
        df2.drop_duplicates(inplace=True)  # Remove duplicate rows
        # Rename columns
        df2.rename(columns={'level_7': 'Species', 'relative_abundance': sample_id}, inplace=True)
        df2.to_excel(output_file_df2, index=False)
    else:
        print("Warning: 'level_7' column not found in the DataFrame.")
    
    # Optionally save the full DataFrame to the primary output file
    df.to_excel(output_file, index=False)

if __name__ == "__main__":
    # Extract command line arguments from snakemake
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    output_file_df1 = snakemake.output[1]
    output_file_df2 = snakemake.output[2]
    
    # Call the function
    parse_mpa4(input_file, output_file, output_file_df1, output_file_df2)
