"""Cluster-Buster compatible motif scoring (Rust extension, see rust/)."""

from os import PathLike

import numpy as np
import numpy.typing as npt

class Sequences:
    """Sequences of a FASTA file, encoded for scoring."""

    def __init__(self, path: str | PathLike[str], mask: bool = False) -> None: ...
    @property
    def names(self) -> list[str]: ...
    @property
    def lengths(self) -> list[int]: ...
    def __len__(self) -> int: ...

def score_motifs(
    sequences: Sequences,
    motif_files: list[str],
    bg_padding: int = 0,
    bg_range: int = 10000,
    n_threads: int = 0,
) -> npt.NDArray[np.float32]:
    """
    Score sequences with Cluster-Buster motif files.

    Returns a float32 array of shape (motifs, sequences) with the top cluster score of each
    sequence, as ``cbust -f 4 -c 0.0 -r bg_range -b bg_padding -t 1`` reports it (0 when there is
    no cluster).
    """
