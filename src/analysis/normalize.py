import pandas as pd
import numpy as np

count_df = pd.read_csv("../../data/GSE119834_raw_counts_GRCh38.p13_NCBI.tsv", delimiter="\t")

size_df = pd.read_csv("../../data/Human.GRCh38.p13.annot.tsv", delimiter="\t")
size_df = size_df['Length']

assert len(count_df) == len(size_df), "Gene ID and length mismatch"

merged_data = count_df.merge(size_df, left_index=True, right_index=True)

rpk_count_df = merged_data.iloc[:, 1:-1].div(merged_data['Length'] / 1000, axis=0)
rpk_sum = rpk_count_df.sum(axis=0)
tpm = rpk_count_df.div(rpk_sum, axis=1) * 1e6 # final TPM-normalized

merged = pd.concat([df['GeneID'], tpm], axis=1)
print (merged)
merged.to_csv("../../data/GSE119834_tpm_normalized.csv")