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

To show results in a notebook (``show_result``), install the ``notebook`` extra: ``pip install "pycistarget[notebook]"``.

Usage
*****

pycistarget needs a cisTarget database, created with `create_cisTarget_databases <https://github.com/aertslab/create_cisTarget_databases>`_.

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

Development
***********

::

	$ uv sync
	$ uv run pytest
	$ uv run ruff check src tests && uv run ruff format src tests
	$ uv run mypy src

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
