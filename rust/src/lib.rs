//! Cluster-Buster compatible motif scoring for creating cisTarget databases.
//!
//! Scores are identical to those of `cbust -f 4 -c 0.0 -r <bg_range> -b <bg_padding> -t 1 [-l]`
//! as read by create_cisTarget_databases (the printed 3 significant digits as float32).

pub mod fasta;
pub mod motif;
pub mod scan;

use std::path::PathBuf;
use std::sync::Arc;

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray2};
use pyo3::exceptions::{PyOSError, PyValueError};
use pyo3::prelude::*;
use rayon::prelude::*;

const BLOCK_SIZE: usize = 256;

/// Round to 3 significant digits like Cluster-Buster's output (`%.3g`, which rounds the exact
/// value to nearest, ties to even, as Rust's formatting does) and read back as float32 like
/// create_cisTarget_databases.
pub fn as_printed(score: f64) -> f32 {
    let text = format!("{score:.2e}");
    text.parse::<f64>().expect("a formatted float is a number") as f32
}

/// Top cluster scores (as printed) of all `sequences` (columns) for all `motifs` (rows).
pub fn score(
    sequences: &fasta::Sequences,
    motifs: &[motif::Motif],
    bg_padding: usize,
    bg_range: usize,
) -> Array2<f32> {
    let n_seqs = sequences.len();
    let blocks: Vec<Vec<f32>> = (0..n_seqs.div_ceil(BLOCK_SIZE))
        .into_par_iter()
        .map_init(scan::Workspace::default, |ws, block| {
            let first = block * BLOCK_SIZE;
            let seqs: Vec<&[u8]> = (first..(first + BLOCK_SIZE).min(n_seqs))
                .map(|i| sequences.seq(i))
                .collect();
            let bgs: Vec<scan::Background> = seqs
                .iter()
                .map(|seq| scan::Background::new(seq, bg_range))
                .collect();
            let bg_refs: Vec<&scan::Background> = bgs.iter().collect();
            let mut out = vec![0.0f32; motifs.len() * seqs.len()];
            let mut top = [None; scan::GROUP];
            // Score sequences of the same length in groups, others one by one.
            let mut i = 0;
            while i < seqs.len() {
                let size = if seqs.len() - i >= scan::GROUP
                    && seqs[i..i + scan::GROUP]
                        .iter()
                        .all(|s| s.len() == seqs[i].len())
                {
                    scan::GROUP
                } else {
                    1
                };
                let group = i..i + size;
                for (j, motif) in motifs.iter().enumerate() {
                    let top = &mut top[..size];
                    scan::top_cluster_scores(
                        &seqs[group.clone()],
                        &bg_refs[group.clone()],
                        motif,
                        bg_padding,
                        ws,
                        top,
                    );
                    for (k, score) in top.iter().enumerate() {
                        if let Some(score) = score {
                            out[j * seqs.len() + i + k] = as_printed(*score);
                        }
                    }
                }
                i += size;
            }
            out
        })
        .collect();

    let mut scores = Array2::<f32>::zeros((motifs.len(), n_seqs));
    for (block, out) in blocks.iter().enumerate() {
        let start = block * BLOCK_SIZE;
        let block_len = out.len() / motifs.len().max(1);
        for (j, mut row) in scores.rows_mut().into_iter().enumerate() {
            row.as_slice_mut().unwrap()[start..start + block_len]
                .copy_from_slice(&out[j * block_len..(j + 1) * block_len]);
        }
    }
    scores
}

/// Sequences of a FASTA file, encoded for scoring.
#[pyclass(frozen, module = "pycistarget._cbust")]
struct Sequences {
    inner: Arc<fasta::Sequences>,
}

#[pymethods]
impl Sequences {
    /// Read a FASTA file. With `mask`, lowercase bases are not scored (like Cluster-Buster `-l`).
    #[new]
    #[pyo3(signature = (path, mask = false))]
    fn new(py: Python<'_>, path: PathBuf, mask: bool) -> PyResult<Self> {
        let bytes = std::fs::read(&path)
            .map_err(|e| PyOSError::new_err(format!("{}: {e}", path.display())))?;
        let inner = py.detach(|| fasta::parse(&bytes, mask));
        Ok(Sequences {
            inner: Arc::new(inner),
        })
    }

    #[getter]
    fn names(&self) -> Vec<String> {
        self.inner.names.clone()
    }

    #[getter]
    fn lengths(&self) -> Vec<usize> {
        self.inner.offsets.windows(2).map(|w| w[1] - w[0]).collect()
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }
}

/// Read Cluster-Buster motif files.
fn read_motifs(paths: &[PathBuf]) -> PyResult<Vec<motif::Motif>> {
    paths
        .iter()
        .map(|path| {
            let bytes = std::fs::read(path)
                .map_err(|e| PyOSError::new_err(format!("{}: {e}", path.display())))?;
            motif::parse(&bytes).map_err(|e| {
                PyValueError::new_err(format!("Could not read motif file {}: {e}", path.display()))
            })
        })
        .collect()
}

/// Score sequences with Cluster-Buster motif files.
///
/// Returns a float32 array of shape (motifs, sequences) with the top cluster score of each sequence,
/// as `cbust -f 4 -c 0.0 -r bg_range -b bg_padding -t 1` reports it (0 when there is no cluster).
#[pyfunction]
#[pyo3(signature = (sequences, motif_files, bg_padding = 0, bg_range = 10000, n_threads = 0))]
fn score_motifs<'py>(
    py: Python<'py>,
    sequences: &Sequences,
    motif_files: Vec<PathBuf>,
    bg_padding: usize,
    bg_range: usize,
    n_threads: usize,
) -> PyResult<Bound<'py, PyArray2<f32>>> {
    let motifs = read_motifs(&motif_files)?;
    let seqs = &sequences.inner;
    for (name, length) in seqs
        .names
        .iter()
        .zip(seqs.offsets.windows(2).map(|w| w[1] - w[0]))
    {
        if length < 2 * bg_padding + 1 {
            return Err(PyValueError::new_err(format!(
                "Sequence {name} should be at least {} bp long (bg_padding = {bg_padding}).",
                2 * bg_padding + 1
            )));
        }
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let scores = py.detach(|| pool.install(|| score(seqs, &motifs, bg_padding, bg_range)));
    Ok(scores.into_pyarray(py))
}

#[cfg(test)]
mod tests {
    use super::as_printed;

    #[test]
    fn rounds_like_printf() {
        // Values printed by cbust (`%.3g`), including exact ties (rounded to even).
        for (score, printed) in [
            (15.4049, "15.4"),
            (1.125, "1.12"),
            (12.25, "12.2"),
            (0.00625, "0.00625"),
            (1234.5, "1.23e+03"),
            (9.996, "10"),
        ] {
            assert_eq!(
                as_printed(score),
                printed.parse::<f32>().unwrap(),
                "{score}"
            );
        }
    }
}

#[pymodule]
fn _cbust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Sequences>()?;
    m.add_function(wrap_pyfunction!(score_motifs, m)?)?;
    Ok(())
}
