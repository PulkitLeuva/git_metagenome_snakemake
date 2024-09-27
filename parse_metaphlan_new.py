import pandas as pd
import glob  # For reading multiple files

def parse_mpa4(input_file):
    # Read the file, skipping the initial metadata lines
    with open(input_file) as f:
        lines = f.readlines()

    # Extract the SampleID from the 3rd line
    sample_id = lines[3].strip().split('\t')[1]
    print("sample_id ===" + sample_id)

    # Load the rest of the file into a DataFrame, skipping the first 5 lines
    df = pd.read_csv(input_file, sep="\t", skiprows=5, header=None)

    # Set the header correctly
    df.columns = ['clade_name', 'NCBI_tax_id', 'relative_abundance', 'additional_species']

    # Split the 'clade_name' column into multiple columns based on '|'
    clade_split = df['clade_name'].str.split('|', expand=True)

    # Rename columns dynamically based on how many levels there are
    clade_split.columns = [f'level_{i+1}' for i in range(clade_split.shape[1])]

    # Concatenate the original dataframe with the new clade levels
    df = pd.concat([clade_split, df[['NCBI_tax_id', 'relative_abundance', 'additional_species']]], axis=1)

    corrupted_level2 = False
    corrupted_level7 = False

    # Check if 'level_2' exists in df and create df1 (Phylum-level data)
    if 'level_2' in df.columns:
        df1 = df[['level_2', 'relative_abundance']]
        df1 = df1.dropna()  # Remove rows with None or NaN values
        df1 = df1.groupby('level_2', as_index=False).agg('first')  # Group by 'level_2' and keep the first non-null value
        df1.rename(columns={'level_2': 'Phylum', 'relative_abundance': sample_id}, inplace=True)
        print(len(df1))
    else:
        print(f"Warning: 'level_2' column not found in the DataFrame for {sample_id}.")
        corrupted_level2 = True
        df1 = pd.DataFrame()  # Return an empty DataFrame if level_2 is not present

    # Check if 'level_7' exists in df and create df2 (Species-level data)
    if 'level_7' in df.columns:
        df2 = df[['level_7', 'relative_abundance']]
        df2 = df2.dropna()  # Remove rows with None or NaN values
        df2 = df2.groupby('level_7', as_index=False).agg('first')  # Group by 'level_7' and keep the first non-null value
        df2.rename(columns={'level_7': 'Species', 'relative_abundance': sample_id}, inplace=True)
        print(len(df2))
    else:
        print(f"Warning: 'level_7' column not found in the DataFrame for {sample_id}.")
        corrupted_level7 = True
        df2 = pd.DataFrame()  # Return an empty DataFrame if level_7 is not present

    return df1, df2, corrupted_level2, corrupted_level7, sample_id


# Initialize empty DataFrames to hold the merged results
main_df1 = pd.DataFrame()
main_df2 = pd.DataFrame()

# Lists to track corrupted files
corrupted_level2_files = []
corrupted_level7_files = []

# List of input files
input_files = glob.glob('/lustre/rohan.b/Gut_Mcrobiome_sequencing_runs/S4_run13_1/git_metagenome_snakemake/metaphlan_out_13/*.txt')  # Update the path and file format

# Loop through each file and process it
for input_file in input_files:
    df1, df2, corrupted_level2, corrupted_level7, sample_id = parse_mpa4(input_file)  # Run your script logic here
    
    # Append to main DataFrames
    main_df1 = pd.concat([main_df1, df1], ignore_index=True)
    main_df2 = pd.concat([main_df2, df2], ignore_index=True)
    
    # Track corrupted files by Sample ID
    if corrupted_level2:
        corrupted_level2_files.append(sample_id)
    if corrupted_level7:
        corrupted_level7_files.append(sample_id)

# Group by 'Phylum' and 'Species' to remove duplicates in the final DataFrames
main_df1 = main_df1.groupby('Phylum', as_index=False).agg('first')
main_df2 = main_df2.groupby('Species', as_index=False).agg('first')

# Save the final merged DataFrames
main_df1.to_csv('/lustre/rohan.b/Gut_Mcrobiome_sequencing_runs/parse_metaphlan/main_output_df1.csv', index=False)
main_df2.to_csv('/lustre/rohan.b/Gut_Mcrobiome_sequencing_runs/parse_metaphlan/main_output_df2.csv', index=False)

# Save the lists of corrupted files to separate CSVs
pd.DataFrame(corrupted_level2_files, columns=['Sample_ID']).to_csv('/lustre/rohan.b/Gut_Mcrobiome_sequencing_runs/parse_metaphlan/corrupted_level2_files.csv', index=False)
pd.DataFrame(corrupted_level7_files, columns=['Sample_ID']).to_csv('/lustre/rohan.b/Gut_Mcrobiome_sequencing_runs/parse_metaphlan/corrupted_level7_files.csv', index=False)
