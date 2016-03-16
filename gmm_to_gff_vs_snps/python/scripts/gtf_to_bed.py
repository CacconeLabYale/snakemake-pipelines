"""Use gffutils to convert the GTF file into a blocked BED12 file."""

import gffutils
from time import time

gtf = snakemake.input.gtf

bed = snakemake.output.bed
gtf_db = snakemake.output.gtf_db

def seconds_to_hms(seconds):
    """Return seconds formatted to hours:minutes:seconds."""
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)


def first_n_features(data, n=5000):
    """Only use the first `n` features of source data."""
    for i, feature in enumerate(gffutils.iterators.DataIterator(data)):
        if i > n:
            break
        yield feature


# Note: this function is optional; if you don't want these IDs then comment out
# the lines at [1] below
def subfeature_handler(f):
    """Given a gffutils.Feature object (which does not yet have its ID assigned),figure out what its ID should be.

    This is intended to be used for CDS, UTR, start_codon, and stop_codon
    features in the Ensembl release 81 GTF files.  I figured a reasonable
    unique ID would consist of the parent transcript and the feature type,
    followed by an autoincrementing number.

    See https://pythonhosted.org/gffutils/database-ids.html#id-spec for
    details and other options.
    """
    return ''.join(
        ['autoincrement:',
         f.attributes['transcript_id'][0],
         '_',
         f.featuretype])




def write_genes_bed(db, bed_path, verbose=True):
    """Convert to bed12 format and write to file."""
    genes = db.features_of_type('gene')

    i = 0

    with open(bed_path, 'w') as out:
        for g in genes:
            i += 1
            if verbose:
                if i == 1:
                    print("Writing records:\n1 ..")
                elif i%2500 == 0:
                    print(i)
            out.write(db.bed12(g.id, name_field='gene_id')+'\n')




# gffutils can spend a lot of time trying to decide on a unique ID for each
# feature. So we have to give it hints of where to look in the attributes.
#
# We also tell it to use our subfeature_handler function for featuretypes with
# no unique IDs.
id_spec = {
    'gene': 'gene_id',
    'transcript': 'transcript_id',

    # [1] These aren't needed for speed, but they do give nicer IDs.
    'CDS': [subfeature_handler],
    'stop_codon': [subfeature_handler],
    'start_codon': [subfeature_handler],
    'exon': [subfeature_handler]
}


try:
    db = gffutils.FeatureDB(gtf_db)
    print("Using existing gffutils database.")
except ValueError as exc:
    if "does not exist" in exc.args[0]:
        print("Building gffutils sqlite database.  This may take a LONG time.")
        start = time()
        db = gffutils.create_db(gtf, dbfn=gtf_db,
                                merge_strategy='create_unique',
                                disable_infer_genes=True, disable_infer_transcripts=True,
                                verbose=True,
                                id_spec=id_spec)
        end = time()
        t = seconds_to_hms(end-start)
        print("gffutils sqlite database build completed in: {elapsed}.".format(elapsed=t))

    else:
        raise exc

write_genes_bed(db=db, bed_path=bed, verbose=True)
