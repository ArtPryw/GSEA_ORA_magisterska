import pytest
import pandas as pd
from ora_script import translate_genes, perform_ora, plot_ora_results

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

def test_perform_ora(setup_files):
    entrez_file, genes_file, tmp_path = setup_files
    translated_file = tmp_path / "translated.csv"
    output_dir = tmp_path / "ora_results"
    output_dir.mkdir()
    
    translate_genes(entrez_file, genes_file, translated_file)
    
    ora_results = perform_ora(translated_file, output_dir, gene_set='KEGG_2019_Human')
    
    assert not ora_results.empty

def test_plot_ora_results(setup_files):
    entrez_file, genes_file, tmp_path = setup_files
    translated_file = tmp_path / "translated.csv"
    output_dir = tmp_path / "ora_results"
    output_dir.mkdir()
    plot_file = tmp_path / "ora_plot.jpg"
    
    translate_genes(entrez_file, genes_file, translated_file)
    perform_ora(translated_file, output_dir, gene_set='KEGG_2019_Human')
    
    ora_results_file = output_dir / "enrichr.KEGG_2019_Human.report.csv"
    plot_ora_results(ora_results_file, plot_file, pval_threshold=0.05)
    
    assert plot_file.exists()
