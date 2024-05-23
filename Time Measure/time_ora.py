import pandas as pd
from ora_script import perform_ora, plot_ora_results
from utils import measure_time, translate_genes

def measure_translate_genes(entrez_file, genes_file, output_file):
    result, execution_time = measure_time(translate_genes, entrez_file, genes_file, output_file)
    return result, execution_time

def measure_perform_ora(translated_file, output_dir, gene_set):
    result, execution_time = measure_time(perform_ora, translated_file, output_dir, gene_set)
    return result, execution_time

def measure_plot_ora_results(input_file, output_file, pval_threshold, figsize):
    result, execution_time = measure_time(plot_ora_results, input_file, output_file, pval_threshold, figsize)
    return result, execution_time

if __name__ == "__main__":
    entrez_file = 'entrez.csv'
    genes_file = 'genes.csv'
    translated_file = 'translated.csv'
    ora_output_dir = 'ora_results'
    plot_file = 'ora_results_plot.jpg'
    
    measure_translate_genes(entrez_file, genes_file, translated_file)
    
    ora_results, _ = measure_perform_ora(translated_file, ora_output_dir, gene_set='KEGG_2019_Human')
    
    ora_results_file =
