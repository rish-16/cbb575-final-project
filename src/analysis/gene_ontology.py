import gseapy as gp
import pandas as pd
import matplotlib.pyplot as plt

go_df = pd.read_csv("../../data/GSE119834.top.table.tsv", delimiter="\t")

genes_of_interest = go_df['Symbol'].tolist()

enr = gp.enrichr(
    gene_list=genes_of_interest,
    gene_sets='GO_Cellular_Component_2018',  # You can use other GO terms like 'GO_Molecular_Function_2018' or 'GO_Cellular_Component_2018'
    organism='Human',  # Change this depending on your organism
    outdir='enrichment_results',
    cutoff=0.05,  # Adjust this cutoff for p-values (e.g., 0.05 for significance)
    no_plot=True
)

# enr.results.head(10).plot.barh(x='Term', y='Adjusted P-value', color='skyblue', figsize=(10, 6))
enr.plot(kind='bubble', title='Top GO Enrichment Terms', color_by='Adjusted P-value')
plt.title('Top GO Enriched Terms')
plt.xlabel('Adjusted P-value')
plt.ylabel('GO Term')
plt.show()