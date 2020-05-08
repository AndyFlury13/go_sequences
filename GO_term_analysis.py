import pandas as pd
import numpy as np

counts = pd.read_csv('/Users/anna/htseq_counts_stranded_reverse.txt', sep='\t')
counts = counts.drop(['Gene_2', 'Gene_3', 'Gene_4'], axis=1)
counts = counts.loc[~((counts['AF16_Control'] == 0) & (counts['AF16_Treatment'] == 0))]
counts = counts.loc[~((counts['HK104_Control'] == 0) & (counts['HK104_Treatment'] == 0))]
counts = counts.replace(0,1)
counts['log2foldchange_control'] = counts['AF16_Control']/counts['HK104_Control']
counts['log2foldchange_control'] = np.log2(counts.log2foldchange_control)
counts['log2foldchange_treatment'] = counts['AF16_Treatment']/counts['HK104_Treatment']
counts['log2foldchange_treatment'] = np.log2(counts.log2foldchange_treatment)
# df2 = df.applymap(np.log2)
print(counts)

# control = control[(control.start >= 6068461) & (control.end <= 6070855)]
# control = control[(control['chrom'] == "AF16_II") | (control['chrom'] == "HK104_II")]
# print(control)

# control.to_csv(r'/Users/anna/hsp4_not_stranded_locus_only_for_htseq_control.txt', sep='\t')


# treatment = treatment[(treatment.start >= 6068461) & (treatment.end <= 6070855)]
# treatment = treatment[(treatment['chrom'] == "AF16_II") | (treatment['chrom'] == "HK104_II")]
# print(treatment)

# treatment.to_csv(r'/Users/anna/hsp4_not_stranded_locus_only_for_htseq_treatment.txt', sep='\t')