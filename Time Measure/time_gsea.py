import pandas as pd
from gsea_script import perform_gsea, plot_gsea_results
from utils import measure_time, translate_genes

def measure_translate_genes(entrez_file, genes_file, output_file):
    result, execution_time = measure_time(translate_genes, entrez_file, genes_file, output_file)
    return result, execution_time

def measure_perform_gsea(translated_file, output_dir, gene_set, processes, permutation_num):
    result, execution_time = measure_time(perform_gsea, translated_file, output_dir, gene_set, processes, permutation_num)
    return result, execution_time

def measure_plot_gsea_results(input_file, output_file, pval_threshold, figsize):
    result, execution_time = measure_time(plot_gsea_results, input_file, output_file, pval_threshold, figsize)
    return result, execution_time

if __name__ == "__main__":
    entrez_file = 'entrez.csv'
    genes_file = 'genes.csv'
    translated_file = 'translated.csv'
    gsea_output_dir = 'gsea_results'
    plot_file = 'gsea_results_plot.jpg'
    
    measure_translate_genes(entrez_file, genes_file, translated_file)
    
    gsea_results, _ = measure_perform_gsea(translated_file, gsea_output_dir, gene_set='KEGG_2019_Human', processes=4, permutation_num=100)
    
    gsea_results_file = f"{gsea_output_dir}/gseapy.prerank.gene_sets.report.csv"
    measure_plot_gsea_results(gsea_results_file, plot_file, pval_threshold=0.05, figsize=(12, 16))
