import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import argparse

def translate_genes(entrez_file, genes_file, output_file):
    """
    Translates gene names in a DataFrame using a mapping provided in another CSV file,
    and saves the result to a new CSV file.

    Parameters:
    entrez_file (str): Path to the CSV file containing the initial and converted gene aliases.
                       This file should have columns 'initial_alias' and 'converted_alias'.
    genes_file (str): Path to the CSV file containing the gene data to be translated.
                      This file should have gene names as its index or in a column.
    output_file (str): Path where the translated DataFrame will be saved as a CSV file.

    Returns:
    None
    """
    # Load the initial and converted gene aliases
    entrez = pd.read_csv(entrez_file)
    
    # Load the genes data, assuming the index column is the first column
    gf = pd.read_csv(genes_file, index_col=0)
    
    # Reset the index and rename the column to 'genes'
    gf.reset_index(inplace=True)
    gf.rename(columns={'index': 'genes'}, inplace=True)
    
    # Create a mapping from initial_alias to converted_alias
    mapping = dict(zip(entrez['initial_alias'], entrez['converted_alias']))
    
    # Translate gene names using the mapping
    gf['translated'] = gf['genes'].map(mapping)
    
    # Drop rows where the 'translated' column is NaN
    gf.dropna(subset=['translated'], inplace=True)
    
    # Set the 'translated' column as the new index
    gf.set_index('translated', inplace=True)
    
    # Drop the original 'genes' column
    gf.drop(columns=['genes'], inplace=True)
    
    # Save the translated DataFrame to a CSV file
    gf.to_csv(output_file)

def perform_gsea(input_file, output_dir, gene_set='KEGG_2019_Human', processes=4, permutation_num=1000):
    """
    Performs GSEA analysis on the given gene expression data.

    Parameters:
    input_file (str): Path to the CSV file containing the translated gene expression data.
    output_dir (str): Directory where the GSEA results will be saved.
    gene_set (str): The gene set to be used for GSEA. Default is 'KEGG_2019_Human'.
    processes (int): Number of CPU processes to use. Default is 4.
    permutation_num (int): Number of permutations for GSEA. Default is 100.

    Returns:
    DataFrame: GSEA results.
    """
    # Load the gene expression data
    df = pd.read_csv(input_file)

    # Set the 'translated' column as the index
    df.set_index('translated', inplace=True)

    # Print basic info about the DataFrame (optional)
    print(df.info())

    # Calculate the rank based on the sum of expression values
    df['rank'] = df.sum(axis=1)

    # Prepare the ranking DataFrame
    rank_df = df[['rank']].reset_index()
    rank_df.columns = ['gene', 'rank']

    # Sort by rank
    rank_df = rank_df.sort_values('rank', ascending=False)

    # Perform GSEA using the specified gene set
    pre_res = gp.prerank(rnk=rank_df, gene_sets=gene_set, outdir=output_dir, processes=processes, permutation_num=permutation_num)

    # Print and return the GSEA results
    print(pre_res.res2d)
    return pre_res.res2d

def plot_gsea_results(input_file, output_file, pval_threshold=0.05, figsize=(12, 16)):
    """
    Creates a bar plot for GSEA results showing gene sets with p-value <= threshold and saves it as a JPEG file.

    Parameters:
    input_file (str): Path to the CSV file containing GSEA results.
    output_file (str): Path where the plot will be saved as a JPEG file.
    pval_threshold (float): p-value threshold for filtering gene sets. Default is 0.05.
    figsize (tuple): Size of the figure for the plot. Default is (12, 16).

    Returns:
    None
    """
    # Load GSEA results data
    df = pd.read_csv(input_file)

    # Filter gene sets based on the p-value threshold
    filtered_df = df[df['fdr'] <= pval_threshold]

    # Sort the results by NES in descending order
    sorted_df = filtered_df.sort_values(by='nes', ascending=False)

    # Create a bar plot
    plt.figure(figsize=figsize)
    plt.barh(sorted_df['Term'], sorted_df['nes'], color='skyblue')
    plt.xlabel('NES (Normalized Enrichment Score)')
    plt.ylabel('Gene Set')
    plt.title(f'GSEA Results for Gene Sets with p-value <= {pval_threshold}')
    plt.gca().invert_yaxis()  # Invert y-axis to have the highest NES at the top
    plt.axvline(x=0, color='gray', linestyle='--')

    # Save the plot as a JPEG file
    plt.savefig(output_file, format='jpg')

    # Display the plot (optional)
    # plt.show()

def main():
    parser = argparse.ArgumentParser(description="Run GSEA analysis and plot results.")
    parser.add_argument('entrez_file', type=str, help='Path to the entrezgene CSV file.')
    parser.add_argument('genes_file', type=str, help='Path to the genes filtered CSV file.')
    parser.add_argument('--gene_set', type=str, default='KEGG_2019_Human', help='The gene set to be used for GSEA. Default is KEGG_2019_Human.')
    parser.add_argument('--translated_file', type=str, default='result_translated_file.csv', help='Name of the output translated CSV file. Default is result_translated_file.csv.')
    parser.add_argument('--gsea_output_dir', type=str, default='gsea_results', help='Directory where the GSEA results will be saved. Default is gsea_results.')
    parser.add_argument('--plot_file', type=str, default='gsea_results_plot.jpg', help='Path where the GSEA results plot will be saved as a JPEG file. Default is gsea_results_plot.jpg.')
    args = parser.parse_args()

    # Translate genes
    translate_genes(args.entrez_file, args.genes_file, args.translated_file)

    # Perform GSEA analysis
    gsea_results = perform_gsea(args.translated_file, args.gsea_output_dir, gene_set=args.gene_set, processes=14, permutation_num=100)

    # Plot GSEA results
    gsea_results_file = f"{args.gsea_output_dir}/gseapy.prerank.gene_sets.report.csv"
    plot_gsea_results(gsea_results_file, args.plot_file, pval_threshold=0.25, figsize=(12, 16))

if __name__ == "__main__":
    main()
