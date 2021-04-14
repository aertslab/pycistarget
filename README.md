# pyCistrome

pyCistrome is a python module to perform motif enrichment analysis in sets of regions with different tools and identify high confidence TF cistromes.

Under development

## cistarget analysis

```python
import pyCistrome.motif_enrichment_cistarget as motenr
import pyranges as pr
import os

#load region sets
path_to_region_sets = '/staging/leuven/stg_00002/lcb/cbravo/Liver/Multiome/pycistopic/GEMSTAT/ChIP/All_summits'
region_sets_files = ['Cebpa_ERR235722_summits_order_by_score_extended_250bp_top5K.bed', 'Foxa1_ERR235786_summits_order_by_score_extended_250bp_top5K.bed', 'Hnf4a_ERR235763_summits_order_by_score_extended_250bp_top5K.bed', 'Onecut1_ERR235752_summits_order_by_score_extended_250bp_top5K.bed']
region_sets = {x.replace('.bed', ''):pr.read_bed(os.path.join(path_to_region_sets, x)) for x in region_sets_files}

#initialize cistarget database
cistarget_database_fname = '/staging/leuven/stg_00002/icistarget-data/make_rankings/v9/CTX_mm10/CTX_mm10_SCREEN3_1kb_bg_with_mask/CTX_mm10_SCREEN3_1kb_bg_with_mask.regions_vs_motifs.rankings.feather'

ctx_db = motenr.cisTargetDatabase(cistarget_database_fname, 'motifs')

#run motif enrichment analysis
cistarget_results = motenr.cistarget_motif_enrichment(ctx_db = ctx_db,
                                                      pr_regions_dict = region_sets,
                                                      species = 'mus_musculus',
                                                      n_cpu = 4)
```