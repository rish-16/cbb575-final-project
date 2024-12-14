import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

rna_seq_data = pd.read_csv('GSE119834_norm_counts_TPM_GRCh38.p13_NCBI.tsv', sep='\t')

test_samples = ['GSM3384758', 'GSM3384759', 'GSM3384760', 'GSM3384761', 'GSM3384762', 'GSM3384763', 'GSM3384764', 'GSM3384765', 'GSM3384766', 'GSM3384767', 'GSM3384768', 'GSM3384769', 'GSM3384770', 'GSM3384771', 'GSM3384772', 'GSM3384773', 'GSM3384774', 'GSM3384775', 'GSM3384776', 'GSM3384777', 'GSM3384778', 'GSM3384779', 'GSM3384780', 'GSM3384781', 'GSM3384782', 'GSM3384783', 'GSM3384784', 'GSM3384785', 'GSM3384786', 'GSM3384787', 'GSM3384788', 'GSM3384789', 'GSM3384790', 'GSM3384791', 'GSM3384792', 'GSM3384793', 'GSM3384794', 'GSM3384795', 'GSM3384796', 'GSM3384797', 'GSM3384798', 'GSM3384799', 'GSM3384800', 'GSM3384801', 'GSM3384802']
control_samples = ['GSM3384847', 'GSM3384849', 'GSM3384850', 'GSM3384851', 'GSM3384852', 'GSM3384853', 'GSM3384854', 'GSM3384855']
selected_genes = [7805, 11309, 81704, 1536, 2321, 3689, 221395, 2212, 5787, 8857, 713, 4481, 9938, 2533, 2214, 117289, 714, 5996, 84868, 5724, 8832, 963, 712, 3043, 9332]

subset_data = rna_seq_data[rna_seq_data['GeneID'].isin(selected_genes)]
subset_data = subset_data.set_index('GeneID')  # Set GeneID as the index

test_data = subset_data[test_samples]
control_data = subset_data[control_samples]

heatmap_data = pd.concat([test_data, control_data], axis=1)

heatmap_data = heatmap_data.apply(lambda x: np.log2(x + 1))

print(heatmap_data.shape)
print(heatmap_data.head())

# plt.figure(figsize=(12, 8))
# sns.heatmap(heatmap_data, cmap='viridis', xticklabels=True, yticklabels=True)
# plt.title('Gene Expression Heatmap')
# plt.xlabel('Samples')
# plt.ylabel('Genes')
# plt.tight_layout()
# plt.show()

