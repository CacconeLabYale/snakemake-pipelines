"""Generate a bed file that represents only the BLAT mappings that do not overlap annotated features in the main species."""

import pybedtools as pbt
import gffutils



gene_models_bed = snakemake.input.gene_models_bed
bed_from_psl = snakemake.input.bed_from_psl

gene_model_subtracted = snakemake.output.gene_model_subtracted



# load BedTools objects
genes_bed = pbt.BedTool(gene_models_bed)
bed_from_psl = pbt.BedTool(bed_from_psl)

# Do subtraction and save to file
subtracted = bed_from_psl.subtract(genes_bed)
subtracted.saveas(gene_model_subtracted)
