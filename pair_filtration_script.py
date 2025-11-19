import pandas as pd

import pandas as pd

# Load file1 (correlation data)
correlation_data = pd.read_csv(r'C:\Users\Sandeep Kushwaha\Desktop\Major_project_data\Major_project_data\New folder\thres3_filtered_edgelist.tsv', sep='\t')

# Load file2 (gene list)
gene_data = pd.read_csv(r'C:\Users\Sandeep Kushwaha\Desktop\Major_project_data\Major_project_data\mok23\ctrl_t24\down_both_last_column.csv', header=None)

# Extract gene names from the second file (assuming they are in the first column)
gene_list = gene_data[0].str.split('|').str[0].unique()  # Split by '|' if necessary and take the unique genes

# Filter the correlation data where both gene1 and gene2 are in the gene_list
filtered_data = correlation_data[correlation_data['gene1'].isin(gene_list) & correlation_data['gene2'].isin(gene_list)]

# Save the filtered data to a new file
filtered_data.to_csv('filtered_d_correlated_genes.csv', index=False)
