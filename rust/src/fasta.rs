//! FASTA reading, following Cluster-Buster's `get_fasta` and `-l` (mask lowercase) option.

/// Code for bases which are not scored (N, other letters and, with masking, lowercase bases).
pub const N: u8 = 4;

/// All sequences of a FASTA file, encoded as 0 (A), 1 (C), 2 (G), 3 (T) or `N`.
pub struct Sequences {
    pub names: Vec<String>,
    /// Sequence `i` is `data[offsets[i]..offsets[i + 1]]`.
    pub offsets: Vec<usize>,
    pub data: Vec<u8>,
}

impl Sequences {
    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    pub fn seq(&self, i: usize) -> &[u8] {
        &self.data[self.offsets[i]..self.offsets[i + 1]]
    }
}

fn translation_table(mask: bool) -> [u8; 256] {
    let mut table = [N; 256];
    for (i, base) in b"ACGT".iter().copied().enumerate() {
        table[base as usize] = i as u8;
        if !mask {
            table[base.to_ascii_lowercase() as usize] = i as u8;
        }
    }
    table
}

/// The first word of a FASTA title line (like `istringstream >> string`).
fn first_word(title: &[u8]) -> String {
    let title = String::from_utf8_lossy(title);
    match title.split_ascii_whitespace().next() {
        Some(word) => word.to_string(),
        None => title.into_owned(),
    }
}

/// Parse FASTA data like Cluster-Buster does.
///
/// Whitespace is skipped, letters are bases, `#` starts a comment until the end of the line
/// and all other characters are ignored. Sequence names are the first word of the title line.
pub fn parse(bytes: &[u8], mask: bool) -> Sequences {
    let table = translation_table(mask);
    let mut names = Vec::new();
    let mut offsets = vec![0];
    let mut data = Vec::with_capacity(bytes.len());
    let mut title: Option<&[u8]> = None;
    let mut i = 0;

    let line_end = |i: usize| {
        bytes[i..]
            .iter()
            .position(|&c| c == b'\n')
            .map_or(bytes.len(), |p| i + p)
    };

    while i < bytes.len() {
        let c = bytes[i];
        if c.is_ascii_alphabetic() {
            data.push(table[c as usize]);
            i += 1;
        } else if c == b'>' {
            if title.is_some() || data.len() != *offsets.last().unwrap() {
                names.push(first_word(title.unwrap_or(b"")));
                offsets.push(data.len());
            }
            let end = line_end(i + 1);
            title = Some(&bytes[i + 1..end]);
            i = end;
        } else if c == b'#' {
            i = line_end(i + 1);
        } else {
            i += 1;
        }
    }
    if title.is_some() || data.len() != *offsets.last().unwrap() {
        names.push(first_word(title.unwrap_or(b"")));
        offsets.push(data.len());
    }
    data.shrink_to_fit();
    Sequences {
        names,
        offsets,
        data,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_like_cluster_buster() {
        let s = parse(b">a x\nACgtN\n# comment ACGT\nA-C\n>b\n\n>c\nT", true);
        assert_eq!(s.names, ["a", "b", "c"]);
        assert_eq!(s.seq(0), [0, 1, N, N, N, 0, 1]);
        assert_eq!(s.seq(1), [] as [u8; 0]);
        assert_eq!(s.seq(2), [3]);
        let s = parse(b">a\nacgtn", false);
        assert_eq!(s.seq(0), [0, 1, 2, 3, N]);
    }
}
