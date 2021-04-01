# pyCistrome

pyCistrome is a python module to perform motif enrichment analysis in sets of regions with different tools and identify high confidence TF cistromes.

Under development

## cistarget analysis

```python
from pyCistrome.motif_enrichment_cistarget import *
from pyCistrome.utils import region_sets_to_signature
import pyranges as pr
import time
start_time = time.time()

#define the cistarget database and motif annotation files
cistarget_database_fname = '/staging/leuven/stg_00002/icistarget-data/make_rankings/v9/CTX_hg38/CTX_hg38_SCREEN3_1kb_bg_with_mask/CTX_hg38_SCREEN3_1kb_bg_with_mask.regions_vs_motifs.rankings.feather'
motif_annotations_fname  = '/staging/leuven/res_00001/databases/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl'

# we will do motif enrichment analysis in regions of topic 5
topic_5 = '/staging/leuven/stg_00002/lcb/sdewin/PhD/De_Winter_hNTorg/ATAC_DOWNSTREAM/pycisTopic/7.TOPIC_BINARIZATION/Topic5.bed'

# first we read the bed file into a pyranges object
pr_topic_5 = pr.read_bed(topic_5)

# next we initialize the cistarget database
ctx_db = cisTargetDatabase(ctx_db_fname = cistarget_database_fname, index_name = 'motifs')

# before we can start the analysis we first have to overlap/query the regions in topic_5 with the reigons in the cistarget database.
# we keep only regions with a minimal fractional overlap of 0.4, i.e. a region in the cistarget database is only kept when at least
# 40 % of that region overlaps with a region in topic 5.
ctx_topic_overlap = overlap_regions_w_cistarget_db(ctx_db, pr_topic_5, min_fraction_overlap = 0.4)

#98.5 % of the regions in topic 5 overlap with regions in the cistarget database

# here we convert the topic 5 PyRanges object to a gene signature object, which will be used to find enriched motifs
input_signatures = region_sets_to_signature(ctx_topic_overlap)

# Finaly we will run the motif enrichment analysis.
# This analysis will also for each motif return sets of regions which are enriched for that motif (i.e. cistromes)

topic_5_cistarget_result = find_enriched_features(  ctx_db = ctx_db,
                                                    region_set_signature = input_signatures[0],
                                                    motif_annotations_fname = motif_annotations_fname,
                                                    get_cistrome = True)
print("--- %s seconds ---" % (time.time() - start_time))

```
> --- 334.348491191864 seconds ---