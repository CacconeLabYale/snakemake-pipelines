"""Functions used in this pipeline."""

import pandas as pd
import numpy as np


def psl_to_dataframe(psl_path):
    """Read PSL into dataframe."""
    head = pd.read_csv(filepath_or_buffer=psl_path, sep='\t',
                    header=None,
                    skiprows=2,
                    nrows=3,)

    head = head.applymap(lambda i: i if pd.notnull(i) else '')

    h1 = head.iloc[0,:]
    h2 = head.iloc[1,:]

    col_names = (h1 + h2).apply(lambda i: i.replace(' ','').replace('-','').replace("'",'').replace(".",'')).str.lower()


    df = pd.read_csv(filepath_or_buffer=psl_path, sep='\t',
                    header=None,
                    names=col_names,
                    skiprows=5,
                    nrows=None,)

    return df



def add_coverage_to_psl(psl):
    """Modify psl to have new column named 'coverage' containing coverage fraction.

    Return: None
    """
    psl['coverage'] = psl.match / psl.qsize

def filter_psl(psl, coverage=0.90, qsize=60):
    """Return new psl dataframe after filtering."""
    psl_ = psl.query(''' coverage >= {coverage} '''.format(coverage=coverage)).query(''' qsize >= {qsize} '''.format(qsize=qsize))

    return psl_

def psl_to_bed(psl):
    """Return new bed dataframe converted from psl dataframe.

    | BED         | PSL           |
    | --          | --            |
    | chrom       | tname         |
    | chromStart  | tstart        |
    | chromEnd    | tend          |
    | name        | qname         |
    | score       | 0             |
    | strand      | strand        |
    | thickStart  | tname         |
    | thickEnd    | tstart        |
    | itemRgb     | "0,0,0" |
    | blockCount  | blockcount    |
    | blockSizes  | blocksizes    |
    | blockStarts | tstarts       |
    """
    bed = pd.DataFrame()

    bed["chrom"] = psl["tname"].copy()
    bed["chromStart"] = psl["tstart"].copy()
    bed["chromEnd"] = psl["tend"].copy()
    bed["name"] = psl["qname"].copy()
    bed["score"] = 0
    bed["strand"] = psl["strand"].copy()
    bed["thickStart"] = psl["tstart"].copy()
    bed["thickEnd"] = psl["tend"].copy()
    bed["itemRgb"] = "0,0,0"
    bed["blockCount"] = psl["blockcount"].copy()
    bed["blockSizes"] = psl["blocksizes"].copy()
    bed["blockStarts"] = psl["tstarts"].copy()

    return bed
