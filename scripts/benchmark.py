import pandas as pd
from pathlib import Path

# Define the main dir containing the benchmark files
main_dir = Path("/home/valentin/Documents/rhomax/benchmark/")

# Initialize an empty list to store individual dataframes
all_data = []

# Iterate through each subdir inside the main dir
for subdir in main_dir.iterdir():
    # Find all tsv files witing the subdir
    for tsv in subdir.glob("*.tsv"):
        # Reat tsv file
        df_temp = pd.read_csv(tsv, sep="\t")
        # Add a new column for source
        df_temp["source"] = tsv.name
        # Append the processed df to the list created previously
        all_data.append(df_temp)

# Create the full df by concatenating all the individual dfs
df = pd.concat(all_data, ignore_index=True)

# Remove irrelevant columns
df = df.drop(columns=["max_uss", "max_pss", "io_in", "io_out"])

# Divide the df into 2 df's for comparison by the source name
df_fm = df[df["source"].str.contains("fm")]
df_tm = df[df["source"].str.contains("tm")]

# Print statistical summary to analyze by
print(df_fm.describe())
print(df_tm.describe())