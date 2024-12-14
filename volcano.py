import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

volcano_data = pd.read_csv('GSE119834.top.table.tsv', sep='\t')

fc_threshold = 1
p_value_threshold = 0.05

volcano_data['neg_log10_padj'] = -np.log10(volcano_data['padj'])

volcano_data['significant'] = (volcano_data['log2FoldChange'].abs() > fc_threshold) & (volcano_data['padj'] < p_value_threshold)

specific_genes = ['2214', '1536', '5724', '72', '58189', '3855']
specific_gene_data = volcano_data[volcano_data['GeneID'].astype(str).isin(specific_genes)]

plt.figure(figsize=(10, 8))
plt.scatter(volcano_data['log2FoldChange'], volcano_data['neg_log10_padj'], alpha=0.5, label="Non-significant")
plt.scatter(volcano_data.loc[volcano_data['significant'], 'log2FoldChange'],
            volcano_data.loc[volcano_data['significant'], 'neg_log10_padj'],
            color='red', alpha=0.7, label="Significant")


for index, row in specific_gene_data.iterrows():
    plt.scatter(row['log2FoldChange'], row['neg_log10_padj'], color='black', marker='o', s=50, label=f"Gene {row['GeneID']}")
    plt.text(row['log2FoldChange'], row['neg_log10_padj'] + 0.2, f"{row['GeneID']}", fontsize=10, color='black', ha='center')

plt.axhline(-np.log10(p_value_threshold), color='blue', linestyle='--', label="p = 0.05")
plt.axvline(-fc_threshold, color='green', linestyle='--', label="Fold Change = -1")
plt.axvline(fc_threshold, color='green', linestyle='--', label="Fold Change = +1")

plt.title("Volcano Plot of Differential Gene Expression")
plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log10 Adjusted P-value")
plt.legend()
plt.tight_layout()
plt.show()