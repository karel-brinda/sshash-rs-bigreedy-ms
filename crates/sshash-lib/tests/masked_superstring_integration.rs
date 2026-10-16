//! End-to-end tests for BiGreedy masked-superstring indexing.
//!
//! A masked superstring (KmerCamel format) is one sequence over `ACGTacgt`
//! where an uppercase base marks the start of a represented k-mer. With a
//! min-one mask, k-mer ids must be the ranks of the uppercase positions in
//! text order, the dictionary's strings are the maximal uppercase runs
//! (overlapping intervals of the text), and membership is content-based.

use sshash_lib::builder::config::BuildConfiguration;
use sshash_lib::builder::dictionary_builder::DictionaryBuilder;
use sshash_lib::kmer::{Kmer, KmerBits};
use sshash_lib::masked_superstring_intervals;
use std::collections::HashMap;

fn lcg_text(seed: u64, len: usize) -> Vec<u8> {
    let mut state = seed | 1;
    let mut next = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (state >> 33) as usize
    };
    (0..len).map(|_| b"acgt"[next() % 4]).collect()
}

/// Apply uppercase runs (inclusive ranges of represented k-mer starts).
fn apply_mask(text: &mut [u8], runs: &[(usize, usize)]) {
    for &(first, last) in runs {
        for c in text.iter_mut().take(last + 1).skip(first) {
            *c = c.to_ascii_uppercase();
        }
    }
}

/// Canonical key of a k-mer given as ASCII: min(seq, revcomp(seq)),
/// case-insensitive. Any canonical representative works as a map key.
fn canonical_key(seq: &[u8]) -> Vec<u8> {
    let fwd: Vec<u8> = seq.iter().map(|c| c.to_ascii_uppercase()).collect();
    let rc: Vec<u8> = fwd
        .iter()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => panic!("invalid base"),
        })
        .collect();
    fwd.min(rc)
}

/// The ground truth: intervals, and the map canonical k-mer -> expected
/// (kmer_id, string_id, kmer_offset). Asserts the mask is min-one.
struct GroundTruth {
    intervals: Vec<(u64, u64)>,
    by_canonical: HashMap<Vec<u8>, (u64, u64, u64)>,
    num_kmers: u64,
}

fn ground_truth(ms: &str, k: usize) -> GroundTruth {
    let intervals: Vec<(u64, u64)> = masked_superstring_intervals(ms.as_bytes(), k)
        .unwrap()
        .iter()
        .map(|iv| (iv.begin, iv.end))
        .collect();

    let mut by_canonical = HashMap::new();
    let mut kmer_id = 0u64;
    for (string_id, &(begin, end)) in intervals.iter().enumerate() {
        for pos in begin..=(end - k as u64) {
            let key = canonical_key(&ms.as_bytes()[pos as usize..pos as usize + k]);
            let prev = by_canonical.insert(key, (kmer_id, string_id as u64, pos));
            assert!(prev.is_none(), "test mask is not min-one: duplicate k-mer at position {pos}");
            kmer_id += 1;
        }
    }
    GroundTruth { intervals, by_canonical, num_kmers: kmer_id }
}

fn build_dict(ms: &str, k: usize, m: usize, force_external: bool) -> sshash_lib::Dictionary {
    let mut config = BuildConfiguration::new(k, m).unwrap();
    config.verbose = false;
    config.num_threads = 1;
    config.force_external_sort = force_external;
    if force_external {
        let dir = std::env::temp_dir().join(format!("sshash_ms_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        config.tmp_dirname = dir;
    }
    DictionaryBuilder::new(config)
        .unwrap()
        .build_from_masked_superstring(ms)
        .unwrap()
}

fn check_all<const K: usize>(dict: &sshash_lib::Dictionary, ms: &str, truth: &GroundTruth)
where
    Kmer<K>: KmerBits,
{
    assert_eq!(dict.num_kmers(), truth.num_kmers);
    assert_eq!(dict.num_strings(), truth.intervals.len() as u64);

    // Every represented k-mer: id = rank of its uppercase position in text
    // order, with the correct string (interval) metadata — on both strands.
    for (string_id, &(begin, end)) in truth.intervals.iter().enumerate() {
        for pos in begin..=(end - K as u64) {
            let bytes = &ms.as_bytes()[pos as usize..pos as usize + K];
            let kmer = Kmer::<K>::from_ascii_unchecked(bytes);
            let expected = truth.by_canonical[&canonical_key(bytes)];

            let fwd = dict.query::<K>(&kmer);
            assert!(fwd.is_found(), "k-mer at position {pos} not found");
            assert_eq!(fwd.kmer_id, expected.0, "wrong kmer_id at position {pos}");
            assert_eq!(fwd.string_id, string_id as u64);
            assert_eq!(fwd.kmer_offset, pos);
            assert_eq!(fwd.kmer_id_in_string, pos - begin);
            assert_eq!(fwd.string_begin, begin);
            assert_eq!(fwd.string_end, end);
            assert_eq!(fwd.kmer_orientation, 1);

            let rev = dict.query::<K>(&kmer.reverse_complement());
            assert_eq!(rev.kmer_id, expected.0);
            assert_eq!(rev.kmer_offset, pos);
            assert_eq!(rev.kmer_orientation, -1);
        }
    }

    // Membership is content-based: an occurrence outside every interval
    // (an unrepresented, lowercase-start occurrence) answers positive iff
    // the same canonical k-mer is represented somewhere else.
    let n = ms.len() as u64;
    for pos in 0..=(n - K as u64) {
        let contained = truth
            .intervals
            .iter()
            .any(|&(b, e)| b <= pos && pos + K as u64 <= e);
        if contained {
            continue;
        }
        let bytes = &ms.as_bytes()[pos as usize..pos as usize + K];
        let kmer = Kmer::<K>::from_ascii_unchecked(bytes);
        let res = dict.query::<K>(&kmer);
        match truth.by_canonical.get(&canonical_key(bytes)) {
            Some(&(id, _, off)) => {
                assert!(res.is_found(), "represented k-mer (other occurrence) at {pos} not found");
                assert_eq!(res.kmer_id, id);
                assert_eq!(res.kmer_offset, off);
            }
            None => {
                assert!(!res.is_found(), "unrepresented k-mer at position {pos} answered positive");
            }
        }
    }

    // Random absent k-mers answer negative.
    let mut state = 12345u64;
    for _ in 0..1000 {
        let text: Vec<u8> = (0..K)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[((state >> 33) & 3) as usize]
            })
            .collect();
        if truth.by_canonical.contains_key(&canonical_key(&text)) {
            continue;
        }
        let kmer = Kmer::<K>::from_ascii_unchecked(&text);
        assert!(!dict.query::<K>(&kmer).is_found());
    }
}

/// The main small test case: a 300-base text with six uppercase runs,
/// including overlapping intervals (overlap < k) and a singleton run.
fn small_masked_superstring() -> String {
    let mut text = lcg_text(42, 300);
    apply_mask(
        &mut text,
        &[
            (0, 19),    // interval [0, 32)
            (25, 39),   // interval [25, 52), overlaps the previous by 7 < k
            (45, 45),   // interval [45, 58), overlaps the previous by 7 < k
            (70, 99),   // interval [70, 112)
            (110, 180), // interval [110, 193), overlaps the previous by 2
            (200, 260), // interval [200, 273)
        ],
    );
    String::from_utf8(text).unwrap()
}

#[test]
fn masked_superstring_lookup_ids_and_membership() {
    let ms = small_masked_superstring();
    let truth = ground_truth(&ms, 13);
    let dict = build_dict(&ms, 13, 5, false);
    check_all::<13>(&dict, &ms, &truth);
}

#[test]
fn masked_superstring_external_sort_path_matches() {
    let ms = small_masked_superstring();
    let truth = ground_truth(&ms, 13);
    let dict = build_dict(&ms, 13, 5, true);
    check_all::<13>(&dict, &ms, &truth);
}

/// With m = 3 there are only 32 canonical minimizers, so a few-thousand-kmer
/// input drives bucket sizes past the skew threshold (64): this exercises
/// the heavy-bucket lookup path and the skew-index build, which must decode
/// the string-id-encoded super-kmer offsets.
#[test]
fn masked_superstring_heavy_buckets() {
    const K: usize = 21;
    let mut text = lcg_text(7, 4200);
    apply_mask(
        &mut text,
        &[
            (0, 3000),    // interval [0, 3021)
            (3005, 4100), // interval [3005, 4121), overlaps the previous by 16 < k
        ],
    );
    let ms = String::from_utf8(text).unwrap();
    let truth = ground_truth(&ms, K);
    assert!(truth.num_kmers > 4000);
    let dict = build_dict(&ms, K, 3, false);
    check_all::<K>(&dict, &ms, &truth);
}

#[test]
fn masked_superstring_streaming_matches_point_lookups() {
    const K: usize = 13;
    let ms = small_masked_superstring();
    let truth = ground_truth(&ms, K);
    let dict = build_dict(&ms, K, 5, false);

    // Stream the whole superstring (uppercased, as a query sequence): every
    // streaming result must equal the point lookup, and positivity must match
    // content-based membership. (Debug builds additionally cross-check every
    // streaming answer inside lookup_with_dict.)
    let query = ms.to_ascii_uppercase();
    let bytes = query.as_bytes();
    let mut engine = dict.create_streaming_query::<K>();
    for i in 0..=(bytes.len() - K) {
        let res = engine.lookup(&bytes[i..i + K]);
        let kmer = Kmer::<K>::from_ascii_unchecked(&bytes[i..i + K]);
        let point = dict.query::<K>(&kmer);
        assert_eq!(res.kmer_id, point.kmer_id, "streaming mismatch at {i}");
        assert_eq!(
            res.is_found(),
            truth.by_canonical.contains_key(&canonical_key(&bytes[i..i + K])),
            "membership mismatch at {i}"
        );
    }
}

#[test]
fn masked_superstring_save_load_roundtrip() {
    const K: usize = 13;
    let ms = small_masked_superstring();
    let truth = ground_truth(&ms, K);
    let dict = build_dict(&ms, K, 5, false);

    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("ms_index");
    dict.save(&prefix).unwrap();
    let loaded = sshash_lib::Dictionary::load(&prefix).unwrap();

    assert_eq!(loaded.k(), dict.k());
    assert_eq!(loaded.num_strings(), dict.num_strings());
    assert_eq!(loaded.num_kmers(), dict.num_kmers());
    check_all::<K>(&loaded, &ms, &truth);
}

#[test]
fn masked_superstring_rejects_invalid_inputs() {
    const K: usize = 13;
    let config = {
        let mut c = BuildConfiguration::new(K, 5).unwrap();
        c.verbose = false;
        c.num_threads = 1;
        c
    };
    let builder = DictionaryBuilder::new(config).unwrap();
    // No represented k-mer start.
    assert!(builder.build_from_masked_superstring(&"acgt".repeat(10)).is_err());
    // Invalid character.
    assert!(builder.build_from_masked_superstring("ACGTNacgtacgtacgt").is_err());
    // Uppercase within the final k-1 positions.
    let mut bad = "acgt".repeat(10);
    bad.replace_range(35..36, "A");
    assert!(builder.build_from_masked_superstring(&bad).is_err());
    // Shorter than k.
    assert!(builder.build_from_masked_superstring("ACGTA").is_err());
}

#[test]
fn masked_superstring_heavy_buckets_use_skew_index() {
    const K: usize = 21;
    let mut text = lcg_text(7, 4200);
    apply_mask(&mut text, &[(0, 3000), (3005, 4100)]);
    let ms = String::from_utf8(text).unwrap();
    let dict = build_dict(&ms, K, 3, false);
    assert!(
        dict.index_ref().skew_index.num_partitions() > 0,
        "test input does not drive any bucket past the skew threshold"
    );
}
