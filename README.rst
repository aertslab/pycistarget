.. image:: https://zenodo.org/badge/351019583.svg
   :target: https://zenodo.org/badge/latestdoi/351019583

pycistarget
============

pycistarget is a python module to perform motif enrichment analysis in sets of regions with different tools and identify high confidence TF cistromes.

Installation
************

pycistarget requires Python 3.10 or newer. To install pycistarget::

	$ pip install git+https://github.com/aertslab/pycistarget

or, from a clone of this repository, using `uv <https://docs.astral.sh/uv/>`_::

	$ git clone https://github.com/aertslab/pycistarget
	$ cd pycistarget
	$ uv sync

Installing from source compiles the motif scoring code (``rust/``), which needs a `Rust toolchain <https://rustup.rs>`_.

To show results in a notebook (``show_result``), install the ``notebook`` extra: ``pip install "pycistarget[notebook]"``.

Usage
*****

pycistarget needs a cisTarget database, created with ``pycistarget create_db`` (see below) or `create_cisTarget_databases <https://github.com/aertslab/create_cisTarget_databases>`_.

Motif enrichment with the cisTarget algorithm (uses the ``*.regions_vs_motifs.rankings.feather`` database)::

	$ pycistarget cistarget \
		--cistarget_db_fname db.regions_vs_motifs.rankings.feather \
		--bed_fname regions.bed \
		--output_folder out \
		--species homo_sapiens

Differentially enriched motifs between foreground and background regions (uses the ``*.regions_vs_motifs.scores.feather`` database)::

	$ pycistarget dem \
		--dem_db_fname db.regions_vs_motifs.scores.feather \
		--foreground_beds foreground.bed \
		--background_beds background_1.bed background_2.bed \
		--output_folder out \
		--species homo_sapiens

To analyse several region sets, pass several bed files: ``--bed_fname`` analyses each bed file separately (one output file per bed file) and ``--one_vs_all_beds`` runs DEM for each bed file against all other bed files, reading the database only once::

	$ pycistarget cistarget --bed_fname topics/*.bed ...
	$ pycistarget dem --one_vs_all_beds topics/*.bed ...

Motif-to-TF annotations are downloaded from the cisTarget webservers, unless a local file is provided with ``--path_to_motif_annotations``.
See ``pycistarget cistarget --help`` and ``pycistarget dem --help`` for all options.

Results written to hdf5 can be read with ``pycistarget.input_output.read_hdf5``.

Creating a database
*******************

``pycistarget create_db`` creates the cisTarget motif databases (scores and rankings) for a set of regions.
It gives the same databases as ``create_cistarget_motif_databases.py`` from `create_cisTarget_databases <https://github.com/aertslab/create_cisTarget_databases>`_, but faster:
motifs are scored with a reimplementation of `Cluster-Buster <https://github.com/weng-lab/cluster-buster>`_ which gives the same scores (on the same platform), and ties in the rankings are broken with the same seeded random permutations.

First create a FASTA file of the regions with background padding (``create_fasta_with_padded_bg_from_bed.sh`` from create_cisTarget_databases), then::

	$ pycistarget create_db \
		--fasta regions.with_1kb_bg_padding.fa \
		--motifs_dir v10nr_clust_public/singletons \
		--motifs motifs.txt \
		--output_prefix db \
		--bg_padding 1000 \
		--mask \
		--threads 32

This writes ``db.motifs_vs_regions.scores.feather``, ``db.regions_vs_motifs.scores.feather`` and ``db.regions_vs_motifs.rankings.feather``.
All scores are kept in memory: the peak memory use is about 3 x 4 bytes x regions x motifs (e.g. 19 GB for 488,000 regions and 3,730 motifs).

Development
***********

::

	$ uv sync
	$ uv run pytest
	$ uv run ruff check src tests && uv run ruff format src tests
	$ uv run mypy src
	$ cd rust && cargo fmt --check && cargo clippy --all-targets && cargo test

Check version
**********************

To check your pycistarget version::

	$ pycistarget --version

or::

	import pycistarget
	pycistarget.__version__


Questions?
**********************

* If you have **technical questions or problems**, such as bug reports or ideas for new features, please open an issue under the issues tab.
* If you have **questions about the interpretation of results or your analysis**, please start a Discussion under the Discussions tab.



Reference
**********************

`Bravo Gonzalez-Blas, C. & De Winter, S. et al. (2022). SCENIC+: single-cell multiomic inference of enhancers and gene regulatory networks <https://www.biorxiv.org/content/10.1101/2022.08.19.504505v1>`_
