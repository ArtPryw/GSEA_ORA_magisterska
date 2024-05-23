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

def perform_ora(input_file, output_dir, gene_set='KEGG_2019_Human'):
    """
    Performs ORA analysis on the given gene expression data.

    Parameters:
    input_file (str): Path to the CSV file containing the translated gene expression data.
    output_dir (str): Directory where the ORA results will be saved.
    gene_set (str): The gene set to be used for ORA. Default is 'KEGG_2019_Human'.

    Returns:
    DataFrame: ORA results.
    """
    # Load the gene expression data
    df = pd.read_csv(input_file)

    # Get the list of translated gene names
    gene_list = df.index.to_list()

    # Perform ORA using the specified gene set
    ora_res = gp.enrichr(gene_list=gene_list, gene_sets=gene_set, outdir=output_dir)

    # Print and return the ORA results
    print(ora_res.res2d)
    return ora_res.res2d

def plot_ora_results(input_file, output_file, pval_threshold=0.05, figsize=(12, 16)):
    """
    Creates a bar plot for ORA results showing gene sets with p-value <= threshold and saves it as a JPEG file.

    Parameters:
    input_file (str): Path to the CSV file containing ORA results.
    output_file (str): Path where the plot will be saved as a JPEG file.
    pval_threshold (float): p-value threshold for filtering gene sets. Default is 0.05.
    figsize (tuple): Size of the figure for the plot. Default is (12, 16).

    Returns:
    None
    """
    # Load ORA results data
    df = pd.read_csv(input_file)

    # Filter gene sets based on the p-value threshold
    filtered_df = df[df['Adjusted P-value'] <= pval_threshold]

    # Sort the results by Combined Score in descending order
    sorted_df = filtered_df.sort_values(by='Combined Score', ascending=False)

    # Create a bar plot
    plt.figure(figsize=figsize)
    plt.barh(sorted_df['Term'], sorted_df['Combined Score'], color='skyblue')
    plt.xlabel('Combined Score')
    plt.ylabel('Gene Set')
    plt.title(f'ORA Results for Gene Sets with p-value <= {pval_threshold}')
    plt.gca().invert_yaxis()  # Invert y-axis to have the highest Combined Score at the top
    plt.axvline(x=0, color='gray', linestyle='--')

    # Save the plot as a JPEG file
    plt.savefig(output_file, format='jpg')

    # Display the plot (optional)
    # plt.show()

def main():
    parser = argparse.ArgumentParser(description="Run ORA analysis and plot results.")
    parser.add_argument('entrez_file', type=str, help='Path to the entrezgene CSV file.')
    parser.add_argument('genes_file', type=str, help='Path to the genes filtered CSV file.')
    parser.add_argument('--gene_set', type=str, default='KEGG_2019_Human', help='The gene set to be used for ORA. Default is KEGG_2019_Human.')
    parser.add_argument('--translated_file', type=str, default='result_translated_file.csv', help='Name of the output translated CSV file. Default is result_translated_file.csv.')
    parser.add_argument('--ora_output_dir', type=str, default='ora_results', help='Directory where the ORA results will be saved. Default is ora_results.')
    parser.add_argument('--plot_file', type=str, default='ora_results_plot.jpg', help='Path where the ORA results plot will be saved as a JPEG file. Default is ora_results_plot.jpg.')
    args = parser.parse_args()

    # Translate genes
    translate_genes(args.entrez_file, args.genes_file, args.translated_file)

    # Perform ORA analysis
    ora_results = perform_ora(args.translated_file, args.ora_output_dir, gene_set=args.gene_set)

    # Plot ORA results
    ora_results_file = f"{args.ora_output_dir}/enrichr.{args.gene_set}.report.csv"
    plot_ora_results(ora_results_file, args.plot_file, pval_threshold=0.05, figsize=(12, 16))

if __name__ == "__main__":
    main()
