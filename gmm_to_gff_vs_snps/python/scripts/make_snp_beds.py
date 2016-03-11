"""Convert SNPs with Yao's deprecated scaffold names to the current official names if needed and write BED files."""

import pandas as pd



do_cleaning = snakemake.params.do_cleaning
scaffold_name_map = snakemake.params.scaffold_name_map
p_thresh = snakemake.params.p_thresh

snp_files = snakemake.input
snp_beds = snakemake.output


def make_name_map(df):
    """Return a kk-to-scaff name dictionary."""
    return {key:value for key,value in zip(df.kk_name.values,df.scaf_name.values)}

def add_chrom_loc_with_map(df, name_map):
    """Add CHROM and LOC columns when not using a name_map.."""
    df["CHROM"] = df.SNP.apply(lambda x: name_map[x.split(':')[1]])
    df["LOC"] = df.SNP.apply(lambda x: int(x.split(':')[2]))

def add_chrom_loc_no_map(df, name_map=None):
    """Add CHROM and LOC columns when not using a name_map."""
    df["CHROM"] = df.SNP.apply(lambda x: x.split(':')[1])
    df["LOC"] = df.SNP.apply(lambda x: int(x.split(':')[2]))




if do_cleaning:
    # create name map
    name_df = pd.read_csv(scaffold_name_map)
    name_map = make_name_map(name_df)

    add_chrom_loc = add_chrom_loc_with_map
else:
    name_map = None

    add_chrom_loc = add_chrom_loc_no_map


for snp_file, snp_bed in zip(snp_files, snp_beds):

    snp_list = pd.read_table(snp_file)
    add_chrom_loc(df=snp_list, name_map=name_map)

    snp_list_passed = snp_list.query(""" P <= {thresh} """.format(thresh=p_thresh))

    bed_table = snp_list_passed[["CHROM","LOC"]].copy()
    bed_table = bed_table.rename(columns={"LOC":"END"}).copy()
    bed_table['START'] = bed_table.END - 1
    bed_table = bed_table[["CHROM","START","END"]].sort_values(by=["CHROM","START"]).copy()


    bed_table.to_csv(snp_bed, header=False, index=False, sep='\t')
