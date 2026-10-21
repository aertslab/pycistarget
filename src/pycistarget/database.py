"""Reading cisTarget databases (Feather files created with create_cisTarget_databases)."""

from __future__ import annotations

import logging
from collections.abc import Iterable
from os import PathLike
from pathlib import Path

import numba
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.feather as pf
import pyarrow.ipc as ipc

log = logging.getLogger(__name__)

_INDEX_COLUMNS = ("motifs", "tracks")
_CHUNK_BYTES = 1 << 30


def _read_column_names(path: Path) -> list[str]:
    try:
        # Feather v2 (Arrow IPC): only the schema is read.
        with ipc.open_file(path) as reader:
            return reader.schema.names
    except pa.ArrowInvalid:
        # Feather v1: fall back to a memory mapped read.
        return pf.read_table(path, memory_map=True).column_names


@numba.njit(parallel=True, cache=True)
def _transpose_into(block, out, start):
    """``out[:, start : start + block.shape[0]] = block.T``, in cache friendly tiles."""
    tile = 64
    n_regions, n_features = block.shape
    for t in numba.prange((n_features + tile - 1) // tile):
        f0 = t * tile
        f1 = min(f0 + tile, n_features)
        for r0 in range(0, n_regions, tile):
            r1 = min(r0 + tile, n_regions)
            for f in range(f0, f1):
                for r in range(r0, r1):
                    out[f, start + r] = block[r, f]


class FeatherCisTargetDatabase:
    """
    A cisTarget scores or rankings database with regions as columns and motifs (or tracks) as rows.

    Only the columns that are requested are read from disk.

    Parameters
    ----------
    fname
        Path to a ``*.regions_vs_motifs.{rankings,scores}.feather`` file.
    """

    def __init__(self, fname: str | PathLike):
        self.fname = Path(fname)
        if not self.fname.is_file():
            raise FileNotFoundError(f'Database "{self.fname}" does not exist.')
        column_names = _read_column_names(self.fname)
        index_columns = [c for c in column_names if c in _INDEX_COLUMNS]
        if len(index_columns) != 1:
            raise ValueError(
                f'"{self.fname}" is not a regions_vs_motifs (or regions_vs_tracks) cisTarget database: '
                f"expected exactly one of the columns {_INDEX_COLUMNS}."
            )
        self.index_column: str = index_columns[0]
        self.region_ids: tuple[str, ...] = tuple(c for c in column_names if c != self.index_column)
        self._region_position = {region: i for i, region in enumerate(self.region_ids)}

    @property
    def total_regions(self) -> int:
        return len(self.region_ids)

    def _in_database(self, region_ids: Iterable[str]) -> list[str]:
        """The regions which are in the database, in the order of the database."""
        return sorted(
            self._region_position.keys() & set(region_ids), key=self._region_position.__getitem__
        )

    def load(self, region_ids: Iterable[str]) -> pd.DataFrame:
        """
        Load the database for a subset of regions.

        Regions which are not in the database are ignored.
        Columns and rows are returned in the order of the database.
        """
        return self._read(self._in_database(region_ids))

    def load_full(self) -> pd.DataFrame:
        return self._read(list(self.region_ids))

    def load_array(
        self, region_ids: Iterable[str], chunk_size: int | None = None
    ) -> tuple[np.ndarray, pd.Index, list[str]]:
        """
        Load the database for a subset of regions as a C-contiguous (motifs x regions) array.

        Regions which are not in the database are ignored. Columns are read in chunks of
        ``chunk_size`` regions (by default about 1 GiB), so memory use is about the size of
        the returned array plus two chunks.

        Returns the array, the motifs (rows) and the regions (columns), both in the order
        of the database.
        """
        found = self._in_database(region_ids)
        features = pd.Index(
            pf.read_table(self.fname, columns=[self.index_column]).column(0).to_pylist(),
            dtype=object,
            name=self.index_column,
        )
        if chunk_size is None:
            # Every read has a fixed cost (parsing the schema of all regions), so read large chunks.
            chunk_size = max(1, _CHUNK_BYTES // (4 * max(1, len(features))))
        values: np.ndarray | None = None
        for start in range(0, len(found), chunk_size):
            columns = found[start : start + chunk_size]
            table = pf.read_table(self.fname, columns=columns, memory_map=False, use_threads=True)
            # Regions x motifs, transposed into the motifs x regions array.
            block = np.stack([table.column(column).to_numpy() for column in columns])
            del table
            if values is None:
                values = np.empty((len(features), len(found)), dtype=block.dtype)
            _transpose_into(block, values, start)
        if values is None:
            values = np.empty((len(features), 0), dtype=np.float32)
        return values, features, found

    def _read(self, columns: list[str]) -> pd.DataFrame:
        table = pf.read_table(
            self.fname, columns=[*columns, self.index_column], memory_map=False, use_threads=True
        )
        df = table.to_pandas().set_index(self.index_column)
        df.columns.name = "regions"
        return df
