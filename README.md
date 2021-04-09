# pyCistrome

pyCistrome is a python module to perform motif enrichment analysis in sets of regions with different tools and identify high confidence TF cistromes.

Under development

## cistarget analysis

```python
from pyCistrome.motif_enrichment_cistarget import *
from pyCistrome.utils import region_sets_to_signature, add_motif_url
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

Add motif images urls and export html

```python
enrichment_table = topic_5_cistarget_result['Enrichment'].copy()

# Add urls to motif images
enrichment_table = add_motif_url(enrichment_table, motif_names = enrichment_table.index.to_list(), base_url = "http://motifcollections.aertslab.org/v9/logos/", key_to_add ='Motif_url')

# Generate html table

def path_to_image_html(path, image_width = 120):
        return '<img src="'+ path + '" width="'+str(image_width)+'" >'

# Remove columns with a lot of text, sort on NES and reorder columns
enrichment_table_small = enrichment_table.drop(['TargetGenes', 'Annotation'], axis = 1)
enrichment_table_small = enrichment_table_small.sort_values('NES', ascending = False)
enrichment_table_small = enrichment_table_small[['Motif_url', 'NES', 'AUC', 'gene_name','RankAtMax']]

#write html file
with open(os.path.join('/staging/leuven/stg_00002/lcb/sdewin/PhD/De_Winter_hNTorg/ATAC_DOWNSTREAM/pycisTopic/cistarget_test', 'topic_5.html'), 'w') as f:
    f.write(enrichment_table_small.to_html(escape = False, formatters = dict(Motif_url=path_to_image_html)))
```

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Motif_url</th>
      <th>NES</th>
      <th>AUC</th>
      <th>gene_name</th>
      <th>RankAtMax</th>
    </tr>
    <tr>
      <th>MotifID</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>hocomoco__SOX3_HUMAN.H11MO.0.B</th>
      <td><img src="http://motifcollections.aertslab.org/v9/logos/hocomoco__SOX3_HUMAN.H11MO.0.B" width="120" ></td>
      <td>8.214834</td>
      <td>0.044119</td>
      <td>NANOG;SMAD1;SOX2;SOX3;SOX4;SOX6;SOX10</td>
      <td>19814</td>
    </tr>
    <tr>
      <th>hocomoco__SOX3_MOUSE.H11MO.0.C</th>
      <td><img src="http://motifcollections.aertslab.org/v9/logos/hocomoco__SOX3_MOUSE.H11MO.0.C" width="120" ></td>
      <td>8.069373</td>
      <td>0.043453</td>
      <td>SMAD1;SOX2;SOX3;SOX4;SOX6;SOX9;SOX10</td>
      <td>19707</td>
    </tr>
    <tr>
      <th>hocomoco__SOX2_HUMAN.H11MO.0.A</th>
      <td><img src="http://motifcollections.aertslab.org/v9/logos/hocomoco__SOX2_HUMAN.H11MO.0.A" width="120" ></td>
      <td>6.877392</td>
      <td>0.037998</td>
      <td>NANOG;NANOGP8;SMAD1;SOX2;SOX3;SOX4;SOX6;SOX9;SOX10;SOX14;SOX17;SOX21;TCF7</td>
      <td>19775</td>
    </tr>
    <tr>
      <th>homer__CCWTTGTY_Sox3</th>
      <td><img src="http://motifcollections.aertslab.org/v9/logos/homer__CCWTTGTY_Sox3" width="120" ></td>
      <td>6.847717</td>
      <td>0.037862</td>
      <td>NANOG;NANOGP8;POU5F1;SMAD1;SOX2;SOX3;SOX4;SOX6;SOX14;SOX17;SOX21</td>
      <td>19899</td>
    </tr>
    <tr>
      <th>homer__CCWTTGTYYB_Sox10</th>
      <td><img src="http://motifcollections.aertslab.org/v9/logos/homer__CCWTTGTYYB_Sox10" width="120" ></td>
      <td>6.642238</td>
      <td>0.036921</td>
      <td>SOX2;SOX3;SOX4;SOX6</td>
      <td>19948</td>
    </tr>
    <tr>
      <th>hocomoco__SP1_HUMAN.H11MO.1.A</th>
      <td><img src="http://motifcollections.aertslab.org/v9/logos/hocomoco__SP1_HUMAN.H11MO.1.A" width="120" ></td>
      <td>6.513225</td>
      <td>0.036331</td>
      <td>BRF1;BRF2;CTCF;EGR1;ESRRA;EZH2;FOS;IRF1;KLF3;KLF5;KLF6;KLF9;KLF12;KLF14;KLF15;KLF16;KLF17;KLF18;MAZ;NELFE;NFYA;NFYB;PATZ1;PBX3;POLR2A;SMARCA4;SP1;SP2;SP3;SP4;SP6;SP7;SP8;SP9;STAT1</td>
      <td>19872</td>
    </tr>
    <tr>
      <th>hocomoco__SOX10_MOUSE.H11MO.1.A</th>
      <td><img src="http://motifcollections.aertslab.org/v9/logos/hocomoco__SOX10_MOUSE.H11MO.1.A" width="120" ></td>
      <td>6.505768</td>
      <td>0.036297</td>
      <td>NANOG;NANOGP8;SMAD1;SOX2;SOX3;SOX4;SOX6;SOX9;SOX10;SOX11</td>
      <td>19468</td>
    </tr>
    <tr>
      <th>transfac_public__M00196</th>
      <td><img src="http://motifcollections.aertslab.org/v9/logos/transfac_public__M00196" width="120" ></td>
      <td>6.312826</td>
      <td>0.035414</td>
      <td>CTCF;EZH2;FOS;IRF1;KLF12;KLF13;KLF14;KLF15;KLF16;KLF17;KLF18;KLF2;KLF3;KLF5;KLF9;MAZ;NFYA;NFYB;PBX3;SP1;SP2;SP3;SP4;SP6;SP7;SP8;SP9;STAT1</td>
      <td>19853</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
  </tbody>
</table>