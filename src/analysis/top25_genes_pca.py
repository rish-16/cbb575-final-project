import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

test_samples = ['GSM3384758', 'GSM3384759', 'GSM3384760', 'GSM3384761', 'GSM3384762', 'GSM3384763', 'GSM3384764', 'GSM3384765', 'GSM3384766', 'GSM3384767', 'GSM3384768', 'GSM3384769', 'GSM3384770', 'GSM3384771', 'GSM3384772', 'GSM3384773', 'GSM3384774', 'GSM3384775', 'GSM3384776', 'GSM3384777', 'GSM3384778', 'GSM3384779', 'GSM3384780', 'GSM3384781', 'GSM3384782', 'GSM3384783', 'GSM3384784', 'GSM3384785', 'GSM3384786', 'GSM3384787', 'GSM3384788', 'GSM3384789', 'GSM3384790', 'GSM3384791', 'GSM3384792', 'GSM3384793', 'GSM3384794', 'GSM3384795', 'GSM3384796', 'GSM3384797', 'GSM3384798', 'GSM3384799', 'GSM3384800', 'GSM3384801', 'GSM3384802']
control_samples = ['GSM3384847', 'GSM3384849', 'GSM3384850', 'GSM3384851', 'GSM3384852', 'GSM3384853', 'GSM3384854', 'GSM3384855']

print (len(test_samples), len(control_samples))

N = 50
topN_gene_ids = pd.read_csv("../../data/GSE119834.top.table.tsv", delimiter="\t")['GeneID'].tolist()[:N]
norm_count_df = pd.read_csv("../../data/GSE119834_raw_counts_GRCh38.p13_NCBI.tsv", delimiter="\t")
topN_norm_count_df = norm_count_df[norm_count_df['GeneID'].isin(topN_gene_ids)]

all_cols = norm_count_df.columns[1:]
print (all_cols)
stem_cell_samples = list(set(all_cols).difference(set(test_samples + control_samples)))

topN_norm_count_df_ctrl = topN_norm_count_df[control_samples]
topN_norm_count_df_test = topN_norm_count_df[test_samples]
topN_norm_count_df_stem = topN_norm_count_df[stem_cell_samples]

tsne_ctrl = TSNE(n_components=2, random_state=42)
tsne_result_ctrl = tsne_ctrl.fit_transform(topN_norm_count_df_ctrl)

tsne_test = TSNE(n_components=2, random_state=42)
tsne_result_test = tsne_test.fit_transform(topN_norm_count_df_test)

tsne_stem = TSNE(n_components=2, random_state=42)
tsne_result_stem = tsne_stem.fit_transform(topN_norm_count_df_stem)

fig = plt.figure(figsize=(8, 6))
plt.scatter(tsne_result_test[:, 0], tsne_result_test[:, 1], label="Test", color="#10ac84")
plt.scatter(tsne_result_ctrl[:, 0], tsne_result_ctrl[:, 1], label="Control", color="#5f27cd")
plt.scatter(tsne_result_stem[:, 0], tsne_result_stem[:, 1], label="Stem Cells", color="#54a0ff")
plt.grid()
plt.legend(fontsize=16)
plt.xlabel('t-SNE1', fontsize=16)
plt.ylabel('t-SNE2', fontsize=16)
plt.title('t-SNE of Top-50 genes (Test vs Control vs Stem Cells)')
plt.show()
fig.savefig("top50_genes_test_v_ctrl_v_stem_tsne.pdf")