**************
SCENIC+ motif collection
**************

The SCENIC+ motif collection includes more than 49,504 motifs from 29 motif collections, with curated TF motif annotations based on 
direct evidence and orthology between species for human, mouse and fly. In order to account for motif redundancy (i.e. the same or a 
very similar version of the same motif can be found in more than one of these collections), we have implemented a new approach to create
non-redundant (or clustered) motif collections using a two-step clustering strategy. 

First, identical PWMs across collection (after rescaling) were merged, resulting in 34,524 motifs. A matrix with motif-to-motif similarity
values was computed using Tomtom (MEME v5.4.1), and motifs with equal length and q-value < 10-40 were merged, resulting in 32,766 motifs 
(unclustered motif collection). For clustering, we used motifs that are similar to at least another motif with q-value < 10-5 (11,526), 
while the remaining were kept as unique motifs, or singlets (9,685). Dimer motifs (1,265) were excluded from the clustering, as well as 
motifs from factorbook and desso, as they do not have direct annotations since they are derived from AI models, as well as motifs with 
an Information Content below 5. With the -log10(Tomtom similiarity)+10-45 matrix as input, we used Seurat (v4.0.3) to normalize,
scale and perform PCA. Using 100 PCs, we performed Leiden clustering with resolution 25, resulting in 199 clusters. Clusters were 
refined by running STAMP (v1.3; using the -cc -sd –chp option) 119 resulting in 1,985 clusters. For each cluster, we used STAMP’s 
consensus motif to generate the motif logo. The TF annotation of a cluster is inferred by merging the TF annotations (direct and orthology) of 
all its members. Overall, the clustered motif collection contains 9,685 singlets, 1,265 dimers and 1,985 clusters (with a mean of 5.8 
motifs per cluster).

The SCENIC+ motif collection contains 8,384, 8,045, 958 annotated clusters for 1,553, 1,357 and 467 TFs (with an average of 5, 6, 2 
motifs per TF) for human, mouse, and fly; respectively. Importantly, motifs are not only annotated based on direct TF-motif experimental
evidence; but also based on orthology (by downloading gene orthologs for the selected species from Ensembl Biomart 105), which permits 
the incorporation of annotations based on experiments in different species. In fact, 433 mouse TFs are only found via orthology, 
augmenting TF-annotations by 47%, as more experiments have been performed in human systems than in mouse. 

The non-redundant motif collection (logos and annotations) can be explored in the table below.

.. raw:: html

    <iframe src="../../source/Motif_collection_collapsed.html" height="500px" width="100%" frameborder="0" "transform": scale(0.25)></iframe>



    


