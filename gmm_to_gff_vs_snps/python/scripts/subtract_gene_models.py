"""Generate a bed file that represents only the BLAT mappings that do not overlap annotated features in the main species."""

import pybedtools as pbt


# input:
#     gene_models=GENE_MODELS,
#     bed_from_psl=rules.filter_psl_to_bed.output.filter_psl_to_bed,
# output:
#     gene_model_subtracted=GENE_MODEL_SUBTRACTED,



def biotype_filter(feature, btype='CDS'):
    """Return True if feature's biotype equals btype."""
    return feature.fields[2] == btype




gene_models = pbt.BedTool(snakemake.input.gene_models)

genes_only = gene_models.filter(biotype_filter, btype='gene').sort()

bed_from_psl = pbt.BedTool(snakemake.input.bed_from_psl)

subtracted = bed_from_psl.subtract(genes_only).sort()



subtracted.saveas(snakemake.output.gene_model_subtracted)
genes_only.saveas(snakemake.output.genes_only_sorted)
