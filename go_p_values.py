#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
from scipy.stats import chisquare


# First we'll save the 'counts' dataframe as its own csv, replacing all the 0s with ones, for some reason. This python code was taken from GO_term_analysis.py

# In[ ]:





# In[2]:


#Save counts as its own csv, create chi-square values
counts = pd.read_csv('htseq_counts_stranded_reverse.txt', sep='\t')
counts = counts.replace(to_replace=0, value=1)


# Now we'll do a chi square test (this one: https://www.khanacademy.org/math/ap-statistics/chi-square-tests/chi-square-tests-two-way-tables/v/chi-square-test-association-independence) on the values. The test_value column specifies whether or not the log2(chi-square value) is greater than 0.

# In[3]:


counts['AF_sum'] = counts['AF16_Control'] + counts['AF16_Treatment']
counts['HK_sum'] = counts['HK104_Control'] + counts['HK104_Treatment']
counts['control_sum'] = counts['AF16_Control'] + counts['HK104_Control']
counts['treatment_sum'] = counts['AF16_Treatment'] + counts['HK104_Treatment']
counts['AFHK_sum'] = counts['AF_sum'] + counts['HK_sum']
counts['AF_control_exp'] = (counts['control_sum'] / counts['AFHK_sum'])*(counts['AF_sum'] / counts['AFHK_sum']) * counts['AFHK_sum']
counts['AF_treatment_exp'] = (counts['treatment_sum'] / counts['AFHK_sum'])*(counts['AF_sum'] / counts['AFHK_sum']) * counts['AFHK_sum']
counts['HK_control_exp'] = (counts['control_sum'] / counts['AFHK_sum'])*(counts['HK_sum'] / counts['AFHK_sum']) * counts['AFHK_sum']
counts['HK_treatment_exp'] = (counts['treatment_sum'] / counts['AFHK_sum'])*(counts['HK_sum'] / counts['AFHK_sum']) * counts['AFHK_sum']
counts['test_value'] = (counts['AF16_Control'] - counts["AF_control_exp"]) ** 2 / counts['AF_control_exp']
counts['test_value'] += (counts['AF16_Treatment'] - counts["AF_treatment_exp"]) ** 2 / counts['AF_treatment_exp']
counts['test_value'] += (counts['HK104_Treatment'] - counts["HK_treatment_exp"]) ** 2 / counts['HK_treatment_exp']
counts['test_value'] += (counts['HK104_Control'] - counts["HK_control_exp"]) ** 2 / counts['HK_control_exp']

counts['test_value'] = np.log2(counts["test_value"])
counts['test_value'] = counts['test_value'] > 0
counts['test_value'] = counts['test_value'].astype(int)
counts


# Now we can load all the relevant data into pandas dataframes. In the below cell we'll load what was the 'counts' dataframe into experiment, and clean it up. We'll use only the gene name, control, and treatment values, and drop all the other columns. I divided control by treatment, took the absolute value, and stuck that in the test_value column.

# In[4]:


experiment = pd.read_csv('experimental_data.csv')

experiment = experiment[['Gene', 'log2foldchange_control', 'log2foldchange_treatment']].copy()
experiment = experiment.rename({'log2foldchange_control': 'control', 'log2foldchange_treatment':'treatment'}, axis=1)
experiment['test_value'] = np.abs(experiment['control'] / experiment['treatment'])
experiment


# Now we'll load 'GO_terms_no_dups.txt' into the df go_genes. I parsed it so that in the 'genes' column I have a list of all the genes corresponding to that particular go squence.

# In[5]:


go_genes = pd.read_csv('briggsae_GO_terms_no_dups.txt',sep='\t')
go_genes = go_genes.rename({go_genes.columns[0]:'go_sequence',
                            go_genes.columns[1]:'genes'}, axis=1)
go_genes


# Now we're going to remove all genes in experiment that don't have a go_sequence.

# In[ ]:


def find_my_GO(string):
    t = go_genes.loc[go_genes['genes'].str.contains(string)]
    return t['go_sequence'].str.cat(sep=',')

experiment['go_sequences'] = np.vectorize(find_my_GO)(experiment['Gene'])


# In[ ]:


experiment = experiment.loc[experiment['go_sequences'] != '']
experiment = experiment.loc[experiment['test_value'] == experiment['test_value']]
experiment


# In[ ]:


def count_genes(go_seq):
    t = experiment.loc[experiment['go_sequences'].str.contains(go_seq)]
    
    return t.shape[0]
    
go_genes['count'] = np.vectorize(count_genes)(go_genes['go_sequence'])
go_genes = go_genes.iloc[:2881, :].copy()

go_genes = go_genes.loc[go_genes['count'] != 0]


# Now we add the p values.

# In[ ]:


def test_sum(go_seq):
    t = experiment.loc[experiment['go_sequences'].str.contains(go_seq)]
    return np.mean(t['test_value'])

go_genes['test_sum'] = np.vectorize(test_sum)(go_genes['go_sequence'])
go_genes

t = experiment.loc[experiment['go_sequences'].str.contains('GO:0000003')]


# In[ ]:


def p_value(go_seq):
    t = experiment.loc[~experiment['go_sequences'].str.contains(go_seq)]
    sample_number = go_genes.loc[go_genes['go_sequence']==go_seq]['count'].values[0]
    test_sum = go_genes.loc[go_genes['go_sequence']==go_seq]['test_sum'].values[0]
    sums = []
    for i in range(1000):
        samples = np.random.choice(a=np.array(t['test_value']), size=sample_number)
        sample_sum = np.mean(samples)
        sums += [sample_sum]
    sums = np.array(sums)
    return np.sum(test_sum >= sums) / 1000



go_genes['p_value'] = np.vectorize(p_value)(go_genes['go_sequence'])


# In[ ]:


go_genes[['go_sequence', 'p_value']].to_csv('pvalues.csv')


# In[ ]:


go_genes


# In[ ]:




