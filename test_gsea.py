import pytest
import pandas as pd
from gsea_script import translate_genes, perform_gsea, plot_gsea_results

@pytest.fixture
def setup_files(tmp_path):
    # Create dummy entrez file
    entrez_file = tmp_path / "entrez.csv"
    entrez_data = pd.DataFrame({
        "initial_alias": ["gene1", "gene2", "gene3"],
        "converted_alias": ["ENTREZ1", "ENTREZ2", "ENTREZ3"]
    })
    entrez_data.to_csv(entrez_file, index=False)

    # Create dummy genes file
    genes_file = tmp_path / "genes.csv"
    genes_data = pd.DataFrame({
        "gene1": [1, 2, 3],
        "gene2": [4, 5, 6],
        "gene3": [7, 8, 9]
    }).T
    genes_data.to_csv(genes_file, header=False)

    return entrez_file, genes_file, tmp_path

def test_translate_genes(setup_files):
    entrez_file, genes_file, tmp_path = setup_files
    output_file = tmp_path / "translated.csv"
    
    translate_genes(entrez_file, genes_file, output_file)
    result = pd.read_csv(output_file, index_col=0)

    assert "ENTREZ1" in result.index
    assert "ENTREZ2" in result.index
    assert "ENTREZ3" in result.index

def test_perform_gsea(setup_files):
    entrez_file, genes_file, tmp_path = setup_files
    translated_file = tmp_path / "translated.csv"
    output_dir = tmp_path / "gsea_results"
    output_dir.mkdir()
    
    translate_genes(entrez_file, genes_file, translated_file)
    
    gsea_results = perform_gsea(translated_file, output_dir, gene_set='KEGG_2019_Human', processes=1, permutation_num=10)
    
    assert not gsea_results.empty

def test_plot_gsea_results(setup_files):
    entrez_file, genes_file, tmp_path = setup_files
    translated_file = tmp_path / "translated.csv"
    output_dir = tmp_path / "gsea_results"
    output_dir.mkdir()
    plot_file = tmp_path / "gsea_plot.jpg"
    
    translate_genes(entrez_file, genes_file, translated_file)
    perform_gsea(translated_file, output_dir, gene_set='KEGG_2019_Human', processes=1, permutation_num=10)
    
    gsea_results_file = output_dir / "gseapy.prerank.gene_sets.report.csv"
    plot_gsea_results(gsea_results_file, plot_file, pval_threshold=0.05)
    
    assert plot_file.exists()
