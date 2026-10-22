//! Cluster-Buster motif files and their preprocessing (`get_cbust_pssm` and `get_matrices`).
//!
//! The arithmetic follows Cluster-Buster exactly (including which steps are done in single
//! precision), so that scores are bit-identical.

/// Cluster-Buster's `minus_infinity` (a float constant), used for log(0) and for N bases.
pub const MINUS_INFINITY: f32 = -1e25;
const PSEUDOCOUNT: f32 = 0.375;
const DEFAULT_EXPECTED_GAP: f64 = 35.0;

/// Matrices of one motif file, ready for scanning.
pub struct Motif {
    /// Rows of all matrices: for each matrix its forward and reverse complement version, with
    /// 5 columns (A, C, G, T, N, padded to 8) and the HMM transition probability added to the
    /// first row.
    pub rows: Vec<[f64; 8]>,
    /// The matrices, in Cluster-Buster's order.
    pub mats: Vec<MatrixRows>,
    pub gap_score: f64,
    pub max_width: usize,
    pub min_width: usize,
}

/// The rows `offset..offset + width` of `Motif::rows`.
#[derive(Clone, Copy)]
pub struct MatrixRows {
    pub offset: usize,
    pub width: usize,
}

struct Matrix {
    rows: Vec<[f64; 8]>,
}

impl Matrix {
    fn width(&self) -> usize {
        self.rows.len()
    }
}

struct RawMatrix {
    rows: Vec<[f32; 4]>,
    weight: f32,
    gap: f32,
}

/// A minimal `std::istream` for the `>>` operations Cluster-Buster uses on motif files.
struct Stream<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl Stream<'_> {
    /// `strm >> c`: the next non-whitespace character.
    fn next_char(&mut self) -> Option<u8> {
        self.skip_whitespace();
        let c = *self.bytes.get(self.pos)?;
        self.pos += 1;
        Some(c)
    }

    fn skip_whitespace(&mut self) {
        while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
    }

    /// `std::getline`: the rest of the current line.
    fn line(&mut self) -> &[u8] {
        let start = self.pos;
        while self.pos < self.bytes.len() && self.bytes[self.pos] != b'\n' {
            self.pos += 1;
        }
        let line = &self.bytes[start..self.pos];
        self.pos = (self.pos + 1).min(self.bytes.len());
        line
    }

    /// `strm >> float`.
    fn float(&mut self) -> Option<f32> {
        self.skip_whitespace();
        let start = self.pos;
        while self.pos < self.bytes.len()
            && matches!(
                self.bytes[self.pos],
                b'0'..=b'9' | b'+' | b'-' | b'.' | b'e' | b'E'
            )
        {
            self.pos += 1;
        }
        std::str::from_utf8(&self.bytes[start..self.pos])
            .ok()?
            .parse()
            .ok()
    }
}

/// `strm >> float` on a string, which gives 0 when it fails.
fn parse_float_or_zero(words: &mut std::str::SplitAsciiWhitespace) -> f32 {
    words
        .next()
        .and_then(|w| {
            Stream {
                bytes: w.as_bytes(),
                pos: 0,
            }
            .float()
        })
        .unwrap_or(0.0)
}

/// Read the next matrix (`get_cbust_pssm`). `Ok(None)` at the end of the file.
fn read_matrix(strm: &mut Stream) -> Result<Option<RawMatrix>, String> {
    let mut weight = 1.0;
    let mut gap = -1.0;
    let mut rows = Vec::new();
    let mut has_title = false;

    while let Some(c) = strm.next_char() {
        if c == b'>' {
            if has_title || !rows.is_empty() {
                strm.pos -= 1;
                break;
            }
            strm.line();
            has_title = true;
            weight = 1.0;
        } else if c == b'#' {
            let line = String::from_utf8_lossy(strm.line()).into_owned();
            let mut words = line.split_ascii_whitespace();
            match words.next() {
                Some("WEIGHT") => weight = parse_float_or_zero(&mut words),
                Some("GAP") => gap = parse_float_or_zero(&mut words),
                _ => {}
            }
        } else {
            strm.pos -= 1;
            let mut row = [0.0; 4];
            for value in row.iter_mut() {
                *value = strm.float().ok_or("could not read matrix row")?;
            }
            rows.push(row);
        }
    }
    if !has_title && rows.is_empty() {
        return Ok(None);
    }
    if rows.is_empty() {
        return Err("empty matrix".to_string());
    }
    Ok(Some(RawMatrix { rows, weight, gap }))
}

/// Pseudocounts, normalisation and log (`normalize_pssm` and `log_pssm`), in single precision.
fn log_probabilities(rows: &[[f32; 4]]) -> Vec<[f32; 4]> {
    rows.iter()
        .map(|row| {
            let mut row = row.map(|count| count + PSEUDOCOUNT);
            let total = row.iter().fold(0.0f32, |total, &p| total + p);
            let inverse = 1.0 / total;
            for p in row.iter_mut() {
                *p *= inverse;
                // Cluster-Buster calls the double precision log here.
                *p = if *p == 0.0 {
                    MINUS_INFINITY
                } else {
                    (*p as f64).ln() as f32
                };
            }
            row
        })
        .collect()
}

fn with_n_column(rows: &[[f32; 4]]) -> Matrix {
    Matrix {
        rows: rows
            .iter()
            .map(|r| {
                let mut row = [MINUS_INFINITY as f64; 8];
                for (value, &p) in row.iter_mut().zip(r) {
                    *value = p as f64;
                }
                row
            })
            .collect(),
    }
}

/// Parse a Cluster-Buster motif file (one or more matrices) and preprocess it (`get_matrices`).
pub fn parse(bytes: &[u8]) -> Result<Motif, String> {
    let mut strm = Stream { bytes, pos: 0 };
    let mut mats = Vec::new();
    let mut weights = Vec::new();
    let mut total_weight = 0.0f32;
    let mut expected_gap = DEFAULT_EXPECTED_GAP;

    while let Some(raw) = read_matrix(&mut strm)? {
        let forward = log_probabilities(&raw.rows);
        // Reverse complement: rotate the 4-column matrix by 180 degrees.
        let reverse: Vec<[f32; 4]> = forward
            .iter()
            .rev()
            .map(|r| [r[3], r[2], r[1], r[0]])
            .collect();
        mats.push(with_n_column(&forward));
        mats.push(with_n_column(&reverse));
        weights.extend([raw.weight, raw.weight]);
        total_weight += raw.weight * 2.0;
        if raw.gap >= 0.0 {
            expected_gap = raw.gap as f64;
        }
    }
    if mats.is_empty() {
        return Err("no matrices".to_string());
    }
    if total_weight == 0.0 {
        return Err("all matrix weights are 0".to_string());
    }

    // tau = 0, so log(1 - tau) = 0.
    let continue_logprob = 1.0f64.ln();
    let cis_start_logprob = continue_logprob - (expected_gap + 1.0).ln();
    for (mat, &weight) in mats.iter_mut().zip(&weights) {
        let w = (weight / total_weight) as f64;
        let w = if w == 0.0 {
            MINUS_INFINITY as f64
        } else {
            w.ln()
        };
        let term = cis_start_logprob + w - mat.width() as f64 * continue_logprob;
        for value in mat.rows[0][..5].iter_mut() {
            *value += term;
        }
    }

    let mut offset = 0;
    let refs = mats
        .iter()
        .map(|mat| {
            let rows = MatrixRows {
                offset,
                width: mat.width(),
            };
            offset += mat.width();
            rows
        })
        .collect();
    Ok(Motif {
        gap_score: (expected_gap / (expected_gap + 1.0)).ln(),
        max_width: mats.iter().map(Matrix::width).max().unwrap(),
        min_width: mats.iter().map(Matrix::width).min().unwrap(),
        rows: mats.into_iter().flat_map(|mat| mat.rows).collect(),
        mats: refs,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_matrices_weights_and_gap() {
        let motif = parse(b">m1\n# WEIGHT 2\n1 0 0 3\n0 4 0 0\n>m2\n#GAP 10\n1 1 1 1\n").unwrap();
        assert_eq!(motif.mats.len(), 4);
        assert_eq!(motif.mats[0].width, 2);
        assert_eq!(motif.mats[2].width, 1);
        assert_eq!(motif.max_width, 2);
        assert_eq!(motif.gap_score, (10.0f64 / 11.0).ln());
        // Row 1 of the reverse complement of m1 is row 0 of m1 complemented: G <- C.
        let p_c = 0.375f32 * (1.0 / 5.5f32);
        assert_eq!(
            motif.rows[motif.mats[1].offset + 1][2],
            (p_c as f64).ln() as f32 as f64
        );
        assert!(parse(b"").is_err());
        assert!(parse(b">empty\n").is_err());
    }
}
