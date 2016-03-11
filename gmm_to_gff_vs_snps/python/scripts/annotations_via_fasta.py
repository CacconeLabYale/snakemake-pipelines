"""Use fasta headers in cuffmerge produced transcript fasta to determine existing gene annotations."""

import numpy as np
import pandas as pd


from pyfasta import Fasta


# ## File paths:

# In[3]:

# define paths to files

tx_fasta = snakemake.input.tx_fasta
ortholog_table = snakemake.input.ortholog_table
tx_orthos_out = snakemake.output.annotations_xls


# In[14]:

tx_names = [x.split(' gene=') for x in Fasta(tx_fasta, key_fn=lambda h: h.lstrip('>')).keys()]


# In[16]:

tx_df = pd.DataFrame(tx_names, columns=['tx_tcons','tx_gene'])


# In[17]:

orthos = pd.read_csv(ortholog_table)


# In[18]:

tx_orthos = pd.merge(left=tx_df, right=orthos,
                     how='inner', on=None,
                     left_on='tx_gene', right_on='Gene stable ID',
                     left_index=False, right_index=False,
                     sort=False, suffixes=('_x', '_y'),
                     copy=True, indicator=False)


# In[19]:

tx_orthos.shape


# In[20]:

# unique GENES
tx_orthos.tx_gene.unique().shape


# In[21]:

# unique TX
tx_orthos.tx_tcons.unique().shape


# In[10]:

tx_orthos.drop("tx_gene", axis=1).to_excel(tx_orthos_out,index=False)
