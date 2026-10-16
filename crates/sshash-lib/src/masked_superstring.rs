//! Masked superstrings in the KmerCamel format
//!
//! Supports BiGreedy masked superstrings as output by `kmercamel compute`:
//! one FASTA record with a single sequence over `ACGTacgt`, where an
//! uppercase base marks the start of a represented k-mer (mask bit 1) and a
//! lowercase base an unrepresented start (mask bit 0). The represented
//! k-mers are exactly the k-mers that start at uppercase positions; with the
//! default min-one mask, each represented k-mer has exactly one uppercase
//! occurrence, as required for the dictionary to assign unique k-mer ids.
//!
//! The maximal runs of represented k-mer starts become the dictionary's
//! strings: possibly overlapping intervals of the one superstring text (see
//! [`crate::spectrum_preserving_string_set::SpectrumPreservingStringSet`]).

/// An interval `[begin, end)` of the superstring whose k-mers are represented.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MaskedSuperstringInterval {
    /// First base of the interval (an uppercase position).
    pub begin: u64,
    /// One past the last base of the interval.
    pub end: u64,
}

/// True if the character marks the start of a represented k-mer.
#[inline]
pub fn is_on_mask_char(c: u8) -> bool {
    c.is_ascii_uppercase()
}

/// Validate that a masked superstring contains only `ACGTacgt`.
pub fn validate_masked_superstring(seq: &[u8]) -> Result<(), String> {
    for (i, &c) in seq.iter().enumerate() {
        if !matches!(c, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') {
            return Err(format!(
                "masked-superstring input contains an invalid character {:?} at position {}",
                c as char, i
            ));
        }
    }
    Ok(())
}

/// Return the maximal runs of represented k-mer starts, as intervals
/// `[begin, end)` over the superstring: a run of uppercase positions
/// `[first, last]` yields the interval `[first, last + k)`, whose k-mers are
/// exactly the k-mers starting at the run's positions. Intervals are sorted
/// by begin, have length at least k, and consecutive intervals overlap by at
/// most k - 2 characters (their uppercase runs are separated by at least one
/// lowercase position), so each k-mer occurrence of the superstring lies in
/// at most one interval.
pub fn masked_superstring_intervals(
    seq: &[u8],
    k: usize,
) -> Result<Vec<MaskedSuperstringInterval>, String> {
    let n = seq.len();
    if n < k {
        return Err("masked-superstring input is shorter than k".to_string());
    }

    if seq[(n - k + 1)..].iter().any(|&c| is_on_mask_char(c)) {
        return Err(
            "masked-superstring input has an uppercase character in its final k-1 positions"
                .to_string(),
        );
    }

    let mut intervals = Vec::new();
    let mut pos = 0usize;
    while pos <= n - k {
        if !is_on_mask_char(seq[pos]) {
            pos += 1;
            continue;
        }
        let begin = pos;
        loop {
            pos += 1;
            if pos > n - k || !is_on_mask_char(seq[pos]) {
                break;
            }
        }
        intervals.push(MaskedSuperstringInterval {
            begin: begin as u64,
            end: (pos + k - 1) as u64,
        });
    }
    if intervals.is_empty() {
        return Err("masked-superstring input has no represented kmer starts".to_string());
    }
    Ok(intervals)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_run() {
        // One run of 3 uppercase starts with k=5: interval [0, 7).
        let seq = b"ACGtacgtac";
        let seq = {
            let mut s = seq.to_vec();
            s[0] = b'A';
            s[1] = b'C';
            s[2] = b'G';
            s
        };
        let intervals = masked_superstring_intervals(&seq, 5).unwrap();
        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0], MaskedSuperstringInterval { begin: 0, end: 7 });
    }

    #[test]
    fn test_two_overlapping_runs() {
        // k=5; runs at positions {0,1} and {4,5}: intervals [0,6) and [4,10),
        // overlapping by 2 < k.
        let seq = b"ACgtACgtac";
        let intervals = masked_superstring_intervals(seq, 5).unwrap();
        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0], MaskedSuperstringInterval { begin: 0, end: 6 });
        assert_eq!(intervals[1], MaskedSuperstringInterval { begin: 4, end: 10 });
    }

    #[test]
    fn test_degenerate_adjacent_runs() {
        // Runs separated by exactly one lowercase position touch without
        // overlapping k-mers: {0}, {3} with k=3 gives [0,3) and [3,6).
        let seq = b"AcgTcgaa";
        let intervals = masked_superstring_intervals(seq, 3).unwrap();
        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0], MaskedSuperstringInterval { begin: 0, end: 3 });
        assert_eq!(intervals[1], MaskedSuperstringInterval { begin: 3, end: 6 });
    }

    #[test]
    fn test_errors() {
        // Shorter than k.
        assert!(masked_superstring_intervals(b"ACg", 5).is_err());
        // Uppercase in the final k-1 positions.
        assert!(masked_superstring_intervals(b"acgtacgTa", 3).is_err());
        // No represented k-mer start.
        assert!(masked_superstring_intervals(b"acgtacgt", 3).is_err());
        // Invalid character.
        assert!(validate_masked_superstring(b"ACGNacgt").is_err());
        assert!(validate_masked_superstring(b"ACGTacgt").is_ok());
    }
}
