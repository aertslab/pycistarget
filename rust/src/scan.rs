//! The top cluster score of a sequence, as `cbust -f 4 -t 1 -c 0` reports it.
//!
//! This follows Cluster-Buster's `scan_seq`, `forward` and `backward` operation by operation
//! (same floating point operations in the same order, same libm functions), but:
//!   - skips the per-matrix contribution scores of the clusters (not printed with `-f 4`),
//!   - skips terms `log1p(exp(d))` which are exactly 0 or too small to change the running score,
//!   - precomputes `matrix value + background` for each base when the background is the same
//!     along the sequence (these are the operands Cluster-Buster adds),
//!   - computes the matrix window sums of up to 8 consecutive positions side by side. A window
//!     sum ending at position n starts from the forward score at n - width, so the sums of
//!     `width` consecutive positions are independent (each still sums its terms in order).

use crate::fasta::N;
use crate::motif::Motif;

/// The local background of a sequence (`get_bg`): log(total / count) of the base at each position.
pub enum Background {
    /// The counting window covers the whole sequence: one value per base (A, C, G, T, N).
    PerBase([f64; 8]),
    PerPosition(Vec<f64>),
}

impl Background {
    /// Base frequencies counted within `bg_range` bases with Laplace's rule. N bases are not counted.
    pub fn new(seq: &[u8], bg_range: usize) -> Self {
        let mut counts = [1u32; 5];
        let mut total = 4u32;
        for &b in &seq[..seq.len().min(bg_range)] {
            if b < N {
                counts[b as usize] += 1;
                total += 1;
            }
        }
        if seq.len() <= bg_range {
            let mut values = [0.0; 8];
            for (value, count) in values.iter_mut().zip(counts) {
                *value = (total as f64 / count as f64).ln();
            }
            return Background::PerBase(values);
        }
        let mut bg = Vec::with_capacity(seq.len());
        for n in 0..seq.len() {
            if seq.len() - n > bg_range {
                let b = seq[n + bg_range];
                if b < N {
                    counts[b as usize] += 1;
                    total += 1;
                }
            }
            bg.push((total as f64 / counts[seq[n] as usize] as f64).ln());
            if n >= bg_range {
                let b = seq[n - bg_range];
                if b < N {
                    counts[b as usize] -= 1;
                    total -= 1;
                }
            }
        }
        Background::PerPosition(bg)
    }
}

/// `matrix value + background` of a base at a position: the operands of Cluster-Buster's window
/// sums.
trait Term {
    /// Row `k` of the motif's matrices.
    fn row(&self, k: usize) -> &[f64; 8];
    fn term(&self, row: &[f64; 8], base: u8, pos: usize) -> f64;
}

/// Precomputed `matrix value + background` per row and base.
struct Table<'a>(&'a [[f64; 8]]);

struct PerPosition<'a> {
    rows: &'a [[f64; 8]],
    bg: &'a [f64],
}

impl Term for Table<'_> {
    #[inline(always)]
    fn row(&self, k: usize) -> &[f64; 8] {
        &self.0[k]
    }

    #[inline(always)]
    fn term(&self, row: &[f64; 8], base: u8, _pos: usize) -> f64 {
        row[base as usize & 7]
    }
}

impl Term for PerPosition<'_> {
    #[inline(always)]
    fn row(&self, k: usize) -> &[f64; 8] {
        &self.rows[k]
    }

    #[inline(always)]
    fn term(&self, row: &[f64; 8], base: u8, pos: usize) -> f64 {
        row[base as usize & 7] + self.bg[pos]
    }
}

/// `score += log1p(exp(s - score))`, skipping the libm calls when the result is known exactly.
#[inline(always)]
fn add_log1p_exp(score: f64, s: f64) -> f64 {
    let d = s - score;
    // exp(d) is 0 for d < -745.2 (e.g. for windows with an N base, which score about -1e25),
    // and log1p(0) = 0.
    if d < -746.0 {
        return score;
    }
    // log1p(exp(d)) <= exp(-40) < 2^-57: less than half a unit in the last place of |score| >= 2^-3.
    if d < -40.0 && score.abs() >= 0.125 {
        return score;
    }
    score + d.exp().ln_1p()
}

/// Marks a window sum which is not used (Cluster-Buster's `continue`): its term is 0.
const SKIP: f64 = f64::NEG_INFINITY;

/// Number of sequences (of the same length) whose forward algorithms are interleaved, so the
/// processor can overlap their (otherwise strictly sequential) libm calls.
pub const GROUP: usize = 4;

/// Reusable buffers.
#[derive(Default)]
pub struct Workspace {
    tables: Vec<Vec<[f64; 8]>>,
    scores: Vec<Vec<f64>>,
    segments: Vec<(usize, usize)>,
    /// Window sums of each matrix and sequence for the current block of positions.
    sums: Vec<[[f64; 8]; GROUP]>,
}

/// Top cluster scores of sequences of the same length (`None` when there is no cluster with a
/// score > 0), written to `out`. `seqs`, `bgs` and `out` have length 1 or `GROUP`.
pub fn top_cluster_scores(
    seqs: &[&[u8]],
    bgs: &[&Background],
    motif: &Motif,
    bg_padding: usize,
    ws: &mut Workspace,
    out: &mut [Option<f64>],
) {
    let len = seqs[0].len();
    assert!(
        seqs.iter().all(|s| s.len() == len) && bgs.len() == seqs.len() && out.len() == seqs.len()
    );
    ws.scores.resize_with(seqs.len(), Vec::new);
    for scores in ws.scores.iter_mut() {
        scores.clear();
        scores.resize(len, 0.0);
    }
    ws.sums.resize(motif.mats.len(), [[0.0; 8]; GROUP]);
    ws.tables.resize_with(seqs.len(), Vec::new);
    // The background type only depends on the length.
    if let Background::PerBase(_) = bgs[0] {
        for (table, bg) in ws.tables.iter_mut().zip(bgs) {
            let Background::PerBase(values) = bg else {
                unreachable!()
            };
            table.clear();
            table.extend(
                motif
                    .rows
                    .iter()
                    .map(|row| std::array::from_fn(|b| row[b] + values[b])),
            );
        }
        let terms: Vec<Table> = ws.tables.iter().map(|t| Table(t)).collect();
        scan(
            seqs,
            &terms,
            motif,
            bg_padding,
            &mut ws.scores,
            &mut ws.segments,
            &mut ws.sums,
            out,
        );
    } else {
        let terms: Vec<PerPosition> = bgs
            .iter()
            .map(|bg| match bg {
                Background::PerPosition(values) => PerPosition {
                    rows: &motif.rows,
                    bg: values,
                },
                Background::PerBase(_) => unreachable!(),
            })
            .collect();
        scan(
            seqs,
            &terms,
            motif,
            bg_padding,
            &mut ws.scores,
            &mut ws.segments,
            &mut ws.sums,
            out,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn scan<T: Term>(
    seqs: &[&[u8]],
    terms: &[T],
    motif: &Motif,
    bg_padding: usize,
    scores: &mut [Vec<f64>],
    segs: &mut Vec<(usize, usize)>,
    sums: &mut [[[f64; 8]; GROUP]],
    out: &mut [Option<f64>],
) {
    let (start, end) = (bg_padding, seqs[0].len() - bg_padding);
    macro_rules! forward {
        ($q:expr) => {
            match motif.min_width {
                8.. => forward::<8, $q, T>(seqs, terms, motif, start, end, scores, sums),
                4.. => forward::<4, $q, T>(seqs, terms, motif, start, end, scores, sums),
                _ => forward::<1, $q, T>(seqs, terms, motif, start, end, scores, sums),
            }
        };
    }
    match seqs.len() {
        GROUP => forward!(GROUP),
        1 => forward!(1),
        n => panic!("scan takes 1 or {GROUP} sequences, not {n}"),
    }

    for q in 0..seqs.len() {
        let scores = &mut scores[q];
        segments(scores, start, end, segs);
        let mut best: Option<f64> = None;
        for &(first, second) in segs.iter() {
            let start = if first + 2 > bg_padding + motif.max_width {
                first + 2 - motif.max_width
            } else {
                bg_padding
            };
            let hi = backward(seqs[q], &terms[q], motif, start, second, scores);
            let score = scores[hi];
            if score > 0.0 && best.is_none_or(|b| score > b) {
                best = Some(score);
            }
        }
        out[q] = best;
    }
}

/// HMM forward algorithm over `start..end` for `Q` sequences of the same length, in blocks of `B`
/// (at most the smallest matrix width) positions.
fn forward<const B: usize, const Q: usize, T: Term>(
    seqs: &[&[u8]],
    terms: &[T],
    motif: &Motif,
    start: usize,
    end: usize,
    scores: &mut [Vec<f64>],
    sums: &mut [[[f64; 8]; GROUP]],
) {
    let seqs: &[&[u8]; Q] = seqs.try_into().unwrap();
    let terms: &[T; Q] = terms.try_into().unwrap();
    let scores: &mut [Vec<f64>; Q] = scores.try_into().unwrap();
    let len = seqs[0].len();
    let mut n0 = start;
    while n0 < end {
        let b = B.min(end - n0);
        // Window sums of the matrices ending at n0..n0 + b. They only use scores before n0.
        let all_fit = n0 + 1 >= motif.max_width && n0 + B <= len;
        for (m, sums) in motif.mats.iter().zip(sums.iter_mut()) {
            let w = m.width;
            for q in 0..Q {
                let (seq, term, scores) = (seqs[q], &terms[q], &scores[q]);
                let sums: &mut [f64; B] = (&mut sums[q][..B]).try_into().unwrap();
                if all_fit {
                    // Lanes beyond b are computed but not used.
                    for (j, sum) in sums.iter_mut().enumerate() {
                        let n = n0 + j;
                        *sum = if n >= w { scores[n - w] } else { 0.0 };
                    }
                    // The matrix is applied backwards, ending at position n.
                    for step in 0..w {
                        let row = term.row(m.offset + w - 1 - step);
                        let pos = n0 - step;
                        let bases: &[u8; B] = seq[pos..pos + B].try_into().unwrap();
                        for (j, sum) in sums.iter_mut().enumerate() {
                            *sum += term.term(row, bases[j], pos + j);
                        }
                    }
                } else {
                    for (j, sum) in sums.iter_mut().enumerate().take(b) {
                        let n = n0 + j;
                        *sum = if n >= w {
                            scores[n - w]
                        } else if n + 1 == w {
                            0.0
                        } else {
                            SKIP
                        };
                        if *sum != SKIP {
                            for step in 0..w {
                                let row = term.row(m.offset + w - 1 - step);
                                *sum += term.term(row, seq[n - step], n - step);
                            }
                        }
                    }
                }
            }
        }

        for j in 0..b {
            let n = n0 + j;
            let mut score = [motif.gap_score; Q];
            if n > 0 {
                for q in 0..Q {
                    score[q] += scores[q][n - 1];
                }
            }
            for s in sums.iter() {
                for q in 0..Q {
                    score[q] = add_log1p_exp(score[q], s[q][j]);
                }
            }
            for q in 0..Q {
                scores[q][n] = score[q];
            }
        }
        n0 += b;
    }
}

/// Segments with maximal score increases, from the forward scores (the second part of
/// Cluster-Buster's `forward`).
fn segments(scores: &[f64], start: usize, end: usize, segs: &mut Vec<(usize, usize)>) {
    segs.clear();
    let mut lo = 0;
    for n in start..end {
        let score = scores[n];
        if score <= scores[lo] {
            lo = n;
        } else if score > scores[n - 1] {
            let mut seg_start = n - 1;
            for i in (0..segs.len()).rev() {
                let (first, second) = segs[i];
                if scores[second] >= score {
                    break;
                }
                if scores[first] < scores[seg_start] {
                    seg_start = first;
                    segs.truncate(i);
                    if seg_start == lo {
                        break;
                    }
                }
            }
            segs.push((seg_start, n));
        }
    }
}

/// HMM backward algorithm from `end` down to `start`; returns the position of the maximum score.
fn backward<T: Term>(
    seq: &[u8],
    term: &T,
    motif: &Motif,
    start: usize,
    end: usize,
    scores: &mut [f64],
) -> usize {
    let mut hi = end;
    for n in (start..=end).rev() {
        let mut score = motif.gap_score;
        if n != end {
            score += scores[n + 1];
        }
        for m in &motif.mats {
            let w = m.width;
            let mut s = if end - n >= w {
                scores[n + w]
            } else if end - n + 1 == w {
                0.0
            } else {
                continue;
            };
            for k in 0..w {
                s += term.term(term.row(m.offset + k), seq[n + k], n + k);
            }
            score = add_log1p_exp(score, s);
        }
        scores[n] = score;
        if score > scores[hi] {
            hi = n;
        }
    }
    hi
}
