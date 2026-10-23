# Changelog

## 1.2.0

### Packaging
- Requires Python >= 3.10; supports numpy 2, pandas 2 and 3.
- Dependencies are declared completely (`numba`, `tables` were missing) and locked with uv (`uv.lock`, replaces `poetry.lock`).
- `ctxcore` is no longer a dependency (it failed to import with recent setuptools and caused very slow DEM runs).
- `pyarrow` is capped below 25: pyarrow 25 and 26 use ~7x more memory when reading a subset of regions from a database.
- `ipython` is optional (`pip install "pycistarget[notebook]"`), only needed for `show_result`.
- Added tests (`uv run pytest`), ruff, mypy and a CI workflow.
- The build backend is maturin: the motif scoring for `create_db` is a Rust extension (`rust/`). Wheels are abi3 (one per platform); installing from source needs a Rust toolchain.

### New features
- `pycistarget create_db`: create cisTarget motif databases (scores and rankings), replacing `create_cistarget_motif_databases.py` from create_cisTarget_databases. Motifs are scored with a reimplementation of Cluster-Buster that gives the same scores as `cbust -f 4 -c 0.0 -r 10000 -b <padding> -t 1 [-l]` (it does the same floating point operations and libm calls, so it matches Cluster-Buster run on the same platform), and the rankings use the same seeded tie-breaking, so the databases are identical.
- `pycistarget cistarget --bed_fname a.bed b.bed ...`: analyse several region sets in one call (one output file per bed file).
- `pycistarget dem --one_vs_all_beds a.bed b.bed ...`: DEM of each bed file against all other bed files, reading the database only once. The result for `a.bed` is the same as `--foreground_beds a.bed --background_beds <all others> --name a_vs_all`.
- `DEM.run` accepts scores loaded once with `DEMDatabase.load_scores`, to run several analyses on the same database regions.
- `--force` (`force=True` for `cisTarget` / `DEM`): run without motif-to-TF annotations when they can not be loaded (e.g. no connection to the cisTarget webservers, or no annotations for the species), instead of failing. The results have empty annotation columns and no cistromes.
- `pycistarget plot barplot` and `pycistarget plot heatmap`: plot the logos of enriched motifs (drawn with logomaker from the Cluster-Buster motif files) with their NES or Log2FC, for one region set (bars) or across region sets (heatmap). The same plots are available in Python (`pycistarget.plotting.plot_logo_barplot` / `plot_logo_heatmap`, returning the figure and axes for customisation). Needs the new `plot` extra (`pip install "pycistarget[plot]"`).

### Performance
- DEM: loading the background regions was quadratic in the number of regions (hours for a 1-vs-all contrast with version 1.1).
- DEM: the Wilcoxon rank-sum test counts the (few thousand) distinct scores of each motif instead of sorting all scores, and the motif hit thresholds are computed without scikit-learn; results are bit-identical.
- DEM: scores are read in chunks into a single matrix (about half the memory).
- DEM 1-vs-all contrast on a 488k region database (379k background regions): ~1 minute and 13 GB.
  `--one_vs_all_beds` for 75 region sets: 14 minutes (one process).
- cisTarget: ~1.9x faster (75 region sets, 8 in parallel: 3.2 instead of 6.1 minutes), lower peak memory.
- numba functions are cached on disk, so they are compiled only once.
- `create_db`: motif scoring needs ~2.2x less CPU time than Cluster-Buster (488k regions x 3,730 motifs: 51 CPU-hours instead of an estimated 115, 3.2 hours with 16 threads); ranking and writing the databases takes ~2 minutes.

### Bug fixes
- `pycistarget dem` failed with a `NameError` when `--genome_annotation` was not given.
- `pycistarget dem` failed when the output folder did not exist; `.gz` is stripped from bed file names.
- `max_bg_regions` was ignored unless `balance_number_of_promoters` was set. **Changes DEM results when `--max_bg_regions` is used.**
- DEM: the motif hit threshold and number of motif hits were computed but not saved; the motif enrichment table now has `Motif_hit_thr` and `Motif_hits` columns.
- DEM: `motifs_to_use` was ignored.
- Cistromes: TFs were matched as a suffix of the annotation, e.g. the SP1 cistrome also contained regions of motifs annotated to MESP1 (also LTF/HLTF, EBF1/SREBF1, ESX1/HESX1, ...). TFs are now matched by name. **Changes the cistromes of these TFs.**
- The order of TFs in annotation columns (e.g. "SIX1, SIX2") was random between runs; it is now sorted.
- Failing to load motif-to-TF annotations raised no error (results silently lacked annotations and cistromes); it now raises an error.
- SSL certificate verification was disabled globally on import.
- `read_hdf5`: species was returned as bytes and `annotation_to_use` as a one-shot iterator.
- `read_hdf5` failed on results without enriched motifs (pandas does not write empty tables).
- Importing pycistarget changed the global pandas option `display.max_colwidth`, and logging was configured inside library functions.

Apart from the changes marked above, results are identical to version 1.1.
