"""Filter hits in PSL then record data and dump output files."""

import pandas as pd
import numpy as np

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")

from python.functions import *

FILTERED_TX_DATAFRAME = snakemake.config["FILTERED_TX_DATAFRAME"]

orthos = pd.read_excel(snakemake.input.orthos)
psl = psl_to_dataframe(snakemake.input.psl)

add_coverage_to_psl(psl)

filtered_psl = filter_psl(psl,
                          coverage=snakemake.params["bed_from_psl_coverage"],
                          qsize=snakemake.params["bed_from_psl_qsize"])

# record some interesting data
hits_per_tx = filtered_psl.groupby('qname').match.agg(len)
qname_qsize = filtered_psl[["qname","qsize"]].drop_duplicates().set_index('qname')
qname_qsize_hits = pd.merge(left=qname_qsize, right=pd.DataFrame(hits_per_tx),
                                how='inner',
                                on=None, left_on=None, right_on=None,
                                left_index=True, right_index=True,
                                sort=False, suffixes=('_x', '_y'), copy=True, indicator=False)

# Ortho Tx length vs number of hits after filtering
sns.jointplot('qsize', 'match', data=qname_qsize_hits)
plt.suptitle("Ortho Tx length vs number of hits after filtering", fontsize=12, y=1)
plt.savefig(filename=snakemake.output.tx_length_vs_hits, dpi=None,
            facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format=None,
            transparent=False, bbox_inches=None, pad_inches=0.1,
            frameon=None)

# filtered hits with annotated orthologs
hits_per_tx_ortho = pd.merge(left=pd.DataFrame(hits_per_tx),
                             right=orthos[["tx_tcons","Gene stable ID","{main_sp} gene stable ID".format(main_sp=snakemake.config["COMMON"]["MAIN_SPECIES"]),"Homology type"]],
                             how='inner',
                             on=None, left_on=None, right_on="tx_tcons",
                             left_index=True, right_index=False)

FILTERED_TX_DATAFRAME["max_hits"] = [hits_per_tx.max()]
FILTERED_TX_DATAFRAME["mean_hits"] = [hits_per_tx.mean()]
FILTERED_TX_DATAFRAME["median_hits"] = [hits_per_tx.median()]
FILTERED_TX_DATAFRAME["filtered_tx"] = [len(hits_per_tx)]
FILTERED_TX_DATAFRAME["filtered_tx_matching_ortho"] = [len(hits_per_tx_ortho.tx_tcons.unique())]
FILTERED_TX_DATAFRAME["filtered_tx_not_matching_ortho"] = FILTERED_TX_DATAFRAME["filtered_tx"] - FILTERED_TX_DATAFRAME["filtered_tx_matching_ortho"]


FILTERED_TX_DATAFRAME.T.to_csv(path_or_buf=snakemake.output.filtered_tx_data, sep='\t',)


bed = psl_to_bed(filtered_psl)
bed.to_csv(path_or_buf=snakemake.output.bed_from_psl, sep='\t', header=False, index=False)
