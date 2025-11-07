import pandas as pd
import os

# FOR REMOVING THE GENE-GENE AND LNC-LNC INTERACTIONS
#directory = r'C:\Users\Sandeep Kushwaha\Desktop\Major_project_data\Major_project_data\New folder'
#if not os.path.exists(directory):
  #  os.makedirs(directory)

# Load the TSV file into a DataFrame
df = pd.read_csv(r'C:\Users\Sandeep Kushwaha\Desktop\pca_latest\project29\wgcna\New folder\CytoscapeInput5-edges-blue.txt', sep='\t', header=None)
# Check if both column 1 and column 2 start with "TCON" and filter those rows out
df = df[~((df[0].str.contains('TCON')) & (df[1].str.contains('TCON')))]
df = df[~((df[0].str.contains('MSTRG')) & (df[1].str.contains('MSTRG')))]
df = df[~((df[0].str.contains('GENE')) & (df[1].str.contains('GENE')))]
df = df[~((df[0].str.contains('STRG')) & (df[1].str.contains('STRG')))]
df = df[~((df[0].str.contains('rna')) & (df[1].str.contains('rna')))]
df = df[~((df[0].str.contains('GENE')) & (df[1].str.contains('MSTRG')))]
df = df[~((df[0].str.contains('gene')) & (df[1].str.contains('gene')))]

print("Number of rows before filtering:", df.shape[0])
print("Number of rows after filtering:", df.shape[0])
print("Data before filtering:")
print(df.head())  # Shows the first few rows before filtering

# Apply your filters

print("Data after filtering:")
print(df.head())  # Shows the first few rows after filtering
print("Number of rows before filtering:", df.shape[0])


# Continue for other filters


# Save the modified DataFrame back to a new TSV file
df.to_csv(r'C:\Users\Sandeep Kushwaha\Desktop\pca_latest\project29\wgcna\New folder\blue_interaction_file.csv', sep='\t', index=False, header=True)
print("DONE")