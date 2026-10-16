//! Spectrum-Preserving String Set (SPSS)
//!
//! This is the core data structure that stores the strings in a compressed form.
//! Each string is stored as a sequence of 2-bit encoded DNA bases.
//!
//! Strings are embedded, possibly overlapping, intervals of one single text,
//! as for a BiGreedy masked superstring, where string j is the j-th maximal
//! run of represented k-mer starts: string j spells `text[begin_j..end_j)`,
//! intervals are sorted by begin, have length at least k, and consecutive
//! intervals overlap by fewer than k characters, so every k-mer occurrence of
//! the text lies in at most one interval. Only the k-mers of the intervals
//! are indexed, and global k-mer ids are assigned in interval order. A plain
//! SPSS (unitigs, simplitigs, eulertigs) is the degenerate case
//! `begin_{j+1} == end_j`.
//!
//! Super-kmer (minimizer) offsets are encoded as `| string-id | relative
//! offset |` so that the same text position remains distinguishable when it
//! anchors k-mers of two adjacent overlapping intervals.

use crate::constants::ceil_log2;
use crate::kmer::{Kmer, KmerBits};
use crate::offsets::{EliasFanoOffsets, OffsetsVector};
use std::io::{self, Read, Write};
use sux::bits::bit_field_vec::BitFieldVec;
use value_traits::slices::{SliceByValue, SliceByValueMut};

/// A decoded super-kmer offset: the absolute text position of a minimizer
/// occurrence plus the string (interval) whose super-kmer recorded it.
#[derive(Clone, Copy, Debug)]
pub struct DecodedSkOffset {
    /// Absolute base position of the minimizer occurrence in the text.
    pub absolute: u64,
    /// The string whose super-kmer recorded this occurrence.
    pub string_id: u64,
    /// Begin of that string (in bases).
    pub string_begin: u64,
    /// End of that string (in bases, exclusive).
    pub string_end: u64,
}

/// The spectrum-preserving string set
///
/// Stores the shared text in a bit-packed format (2 bits per base) plus the
/// string intervals over it:
/// - `begins`: interval begins, strictly increasing, O(1) access;
/// - `ends`: a leading zero followed by the interval ends, Elias-Fano encoded
///   for successor queries (locate the interval fully containing a k-mer);
/// - `cum`: a leading zero followed by the cumulative k-mer counts, giving
///   each string's first global k-mer id.
pub struct SpectrumPreservingStringSet {
    /// Bit-packed text (2 bits per base), shared by all strings
    strings: Vec<u8>,
    /// Interval begins (in bases), strictly increasing
    begins: BitFieldVec,
    /// Leading zero + interval ends (in bases), strictly increasing
    ends: EliasFanoOffsets,
    /// Leading zero + cumulative k-mer counts per interval
    cum: EliasFanoOffsets,
    /// Bits of the relative-offset field in an encoded super-kmer offset
    rel_offset_bits: u32,
    /// K-mer size of the strings
    k: usize,
    /// M (minimizer) size
    m: usize,
}

impl SpectrumPreservingStringSet {
    /// Create a new empty SPSS
    pub fn new(k: usize, m: usize) -> Self {
        Self {
            strings: vec![0u8; 8],
            begins: BitFieldVec::new(1, 0),
            ends: EliasFanoOffsets::from_vec(&[0]),
            cum: EliasFanoOffsets::from_vec(&[0]),
            rel_offset_bits: 1,
            k,
            m,
        }
    }

    /// Create a new SPSS from existing strings and endpoint offsets (the
    /// plain-SPSS path: string j spells `text[offsets[j]..offsets[j+1])`).
    ///
    /// # Arguments
    /// * `strings` - Encoded string data (2-bit packed)
    /// * `offsets` - Offset vector for string boundaries (first entry 0)
    /// * `k` - K-mer size
    /// * `m` - Minimizer size
    pub fn from_parts(strings: Vec<u8>, offsets: OffsetsVector, k: usize, m: usize) -> Self {
        let endpoints = offsets.into_inner();
        assert!(!endpoints.is_empty() && endpoints[0] == 0);
        let begins = endpoints[..endpoints.len() - 1].to_vec();
        let ends = endpoints[1..].to_vec();
        Self::from_interval_vecs(strings, begins, ends, k, m)
    }

    /// Create a new SPSS whose strings are the given intervals of the text.
    ///
    /// Intervals must be sorted with strictly increasing begins and ends,
    /// have length at least k, and consecutive intervals must overlap by
    /// fewer than k characters.
    pub fn from_interval_vecs(
        mut strings: Vec<u8>,
        begins_v: Vec<u64>,
        ends_v: Vec<u64>,
        k: usize,
        m: usize,
    ) -> Self {
        assert_eq!(begins_v.len(), ends_v.len());
        // Pad strings with 8 zero bytes to allow safe unaligned u64 reads
        // near the end without out-of-bounds access.
        strings.resize(strings.len() + 8, 0);

        let p = begins_v.len();
        if p == 0 {
            let mut empty = Self::new(k, m);
            empty.strings = strings;
            return empty;
        }

        let mut max_len = 0u64;
        for j in 0..p {
            let len = ends_v[j] - begins_v[j];
            assert!(len >= k as u64, "string {} is shorter than k", j);
            if len > max_len {
                max_len = len;
            }
            if j > 0 {
                assert!(begins_v[j] > begins_v[j - 1] && ends_v[j] > ends_v[j - 1]);
                debug_assert!(
                    begins_v[j] + k as u64 > ends_v[j - 1],
                    "consecutive intervals must overlap by fewer than k characters"
                );
            }
        }

        let universe = *ends_v.last().unwrap();
        let begins_width = ceil_log2(universe + 1).max(1);
        let mut begins = BitFieldVec::new(begins_width, p);
        for (j, &b) in begins_v.iter().enumerate() {
            begins.set_value(j, b as usize);
        }

        let mut ends_arr = Vec::with_capacity(p + 1);
        ends_arr.push(0u64);
        ends_arr.extend_from_slice(&ends_v);
        let ends = EliasFanoOffsets::from_vec(&ends_arr);

        let mut cum_arr = Vec::with_capacity(p + 1);
        cum_arr.push(0u64);
        for j in 0..p {
            let num_kmers = ends_v[j] - begins_v[j] - k as u64 + 1;
            cum_arr.push(cum_arr[j] + num_kmers);
        }
        let cum = EliasFanoOffsets::from_vec(&cum_arr);

        // Relative minimizer offsets within a string span [0, len - m].
        debug_assert!(m < k);
        let rel_offset_bits = ceil_log2(max_len - m as u64 + 1).max(1) as u32;
        let id_bits = ceil_log2(p as u64) as u32;
        assert!(
            id_bits + rel_offset_bits <= 64,
            "encoded super-kmer offset does not fit within 64 bits"
        );

        Self { strings, begins, ends, cum, rel_offset_bits, k, m }
    }

    /// Get the string offsets (begin, end) for a string ID
    #[inline]
    pub fn string_offsets(&self, string_id: u32) -> (u64, u64) {
        let id = string_id as usize;
        let begin = self.begins.index_value(id) as u64;
        let end = self.ends.access(id + 1);
        (begin, end)
    }

    /// Get the number of strings stored
    #[inline]
    pub fn num_strings(&self) -> u64 {
        self.begins.len() as u64
    }

    /// Get the starting offset (in bases) of a string
    #[inline]
    pub fn string_offset(&self, string_id: u64) -> u64 {
        self.begins.index_value(string_id as usize) as u64
    }

    /// First global k-mer id of a string: the number of k-mers of the
    /// strings preceding it.
    #[inline]
    pub fn first_kmer_id(&self, string_id: u64) -> u64 {
        self.cum.access(string_id as usize)
    }

    /// Total number of indexed k-mers (over all strings).
    #[inline]
    pub fn num_kmers(&self) -> u64 {
        self.cum.access(self.begins.len())
    }

    /// Get the k-mer size
    pub fn k(&self) -> usize {
        self.k
    }

    /// Get the minimizer size
    pub fn m(&self) -> usize {
        self.m
    }

    /// Bits of the relative-offset field in an encoded super-kmer offset.
    #[inline]
    pub fn rel_offset_bits(&self) -> u32 {
        self.rel_offset_bits
    }

    /// Width of an encoded super-kmer offset: `| string-id | relative offset |`.
    #[inline]
    pub fn num_bits_per_offset(&self) -> usize {
        ceil_log2(self.num_strings()) + self.rel_offset_bits as usize
    }

    /// Encode a minimizer occurrence at absolute position `absolute`,
    /// recorded by the super-kmer of string `string_id` (which begins at
    /// `string_begin`), as `| string-id | relative offset |`.
    #[inline]
    pub fn encode_offset(&self, absolute: u64, string_begin: u64, string_id: u64) -> u64 {
        debug_assert!(absolute >= string_begin);
        let relative = absolute - string_begin;
        debug_assert!(self.rel_offset_bits == 64 || relative < (1u64 << self.rel_offset_bits));
        (string_id << self.rel_offset_bits) | relative
    }

    /// Decode an encoded super-kmer offset.
    #[inline]
    pub fn decode_offset(&self, encoded: u64) -> DecodedSkOffset {
        let relative = encoded & ((1u64 << self.rel_offset_bits) - 1);
        let string_id = encoded >> self.rel_offset_bits;
        debug_assert!((string_id as usize) < self.begins.len());
        let string_begin = self.begins.index_value(string_id as usize) as u64;
        let string_end = self.ends.access(string_id as usize + 1);
        DecodedSkOffset {
            absolute: string_begin + relative,
            string_id,
            string_begin,
            string_end,
        }
    }

    /// Get raw access to the packed string data (for buffered reading)
    #[inline]
    pub fn strings_data(&self) -> &[u8] {
        &self.strings
    }

    /// Get the total number of bases covered by the strings (the end of the
    /// last interval; any trailing unindexed bases of the text are excluded).
    pub fn total_bases(&self) -> u64 {
        let p = self.begins.len();
        if p == 0 { 0 } else { self.ends.access(p) }
    }

    /// Locate the string that fully contains the k-mer starting at
    /// `kmer_pos`. Since intervals are sorted by begin and overlap by fewer
    /// than k characters, it can only be the string with the smallest end
    /// `>= kmer_pos + k`. Returns `(string_id, string_begin, string_end)`,
    /// or `None` when no string contains the k-mer.
    #[inline]
    pub fn locate_kmer(&self, kmer_pos: u64, k: u64) -> Option<(u64, u64, u64)> {
        let (string_id, _, end) = self.ends.locate_with_end(kmer_pos + k - 1)?;
        let begin = self.begins.index_value(string_id as usize) as u64;
        (kmer_pos >= begin).then_some((string_id, begin, end))
    }

    /// Locate a string containing a given absolute position.
    /// Returns `(string_id, string_begin)` or None if out of bounds.
    /// With overlapping intervals a position can belong to two strings; the
    /// one with the smallest end is returned.
    #[inline]
    pub fn locate(&self, absolute_pos: u64) -> Option<(u64, u64)> {
        self.locate_with_end(absolute_pos).map(|(id, begin, _)| (id, begin))
    }

    /// Locate a string containing a given absolute position, returning
    /// `(string_id, string_begin, string_end)`.
    #[inline]
    pub fn locate_with_end(&self, absolute_pos: u64) -> Option<(u64, u64, u64)> {
        let (string_id, _, end) = self.ends.locate_with_end(absolute_pos)?;
        let begin = self.begins.index_value(string_id as usize) as u64;
        (absolute_pos >= begin).then_some((string_id, begin, end))
    }

    /// Get memory usage in bits (excludes padding bytes)
    pub fn num_bits(&self) -> u64 {
        ((self.strings.len() - 8) as u64) * 8 + (self.offsets_bytes() as u64) * 8
    }

    /// Get the byte size of the packed strings data (excludes padding)
    pub fn strings_bytes(&self) -> usize {
        self.strings.len() - 8
    }

    /// Access the raw 2-bit-packed strings buffer. Intended for crates that
    /// need to clone the SPSS contents into their own representation.
    #[inline]
    pub fn strings_raw(&self) -> &[u8] {
        &self.strings
    }

    /// Number of string-endpoint entries (`num_strings + 1`), as for a plain
    /// SPSS whose strings partition the text.
    #[inline]
    pub fn offsets_len(&self) -> usize {
        self.begins.len() + 1
    }

    /// Get the i-th string endpoint (0..=num_strings), in bases: entry 0 is
    /// the begin of string 0 and entry i (i > 0) the end of string i - 1.
    /// For a plain SPSS this is the classic offsets array.
    #[inline]
    pub fn offset_at(&self, i: usize) -> u64 {
        if i == 0 {
            if self.begins.is_empty() { 0 } else { self.begins.index_value(0) as u64 }
        } else {
            self.ends.access(i)
        }
    }

    /// Get the byte size of the string-boundary structures (begins + ends +
    /// cumulative k-mer counts)
    pub fn offsets_bytes(&self) -> usize {
        std::mem::size_of_val(self.begins.as_slice())
            + self.ends.num_bytes() as usize
            + self.cum.num_bytes() as usize
    }

    /// Get the length of a specific string in bases
    pub fn string_length(&self, string_id: u64) -> usize {
        let (begin, end) = self.string_offsets(string_id as u32);
        (end - begin) as usize
    }

    /// Decode a k-mer from a specific position in a string
    ///
    /// Uses word-level loads from the packed buffer for efficiency.
    #[inline]
    pub fn decode_kmer<const K: usize>(&self, string_id: u64, kmer_pos: usize) -> Kmer<K>
    where
        Kmer<K>: KmerBits,
    {
        let begin = self.string_offset(string_id);
        let start_base = (begin as usize) + kmer_pos;

        let byte_offset = start_base / 4;
        let bit_shift = (start_base % 4) * 2;
        let needed_bits = K * 2;

        // Read enough bytes to cover K bases + potential misalignment
        // K bases need ceil((K*2 + bit_shift) / 8) bytes
        let needed_bytes = (needed_bits + bit_shift).div_ceil(8);

        if needed_bytes <= 8 {
            // Single u64 read suffices (K <= ~28-31 bases, common case)
            let mut buf = [0u8; 8];
            let avail = self.strings.len().saturating_sub(byte_offset).min(8);
            buf[..avail].copy_from_slice(&self.strings[byte_offset..byte_offset + avail]);
            let raw = u64::from_le_bytes(buf);
            let shifted = raw >> bit_shift;
            let mask = if needed_bits >= 64 { u64::MAX } else { (1u64 << needed_bits) - 1 };
            Kmer::<K>::new(<Kmer<K> as KmerBits>::from_u64(shifted & mask))
        } else {
            // Need u128 for K > 28 or large bit_shift.
            // For K=63 (needed_bits=126), when bit_shift > 2 we need more
            // than 128 bits from the byte stream. Load 16 bytes as u128
            // plus an extra byte to cover the overflow.
            let mut buf = [0u8; 17];
            let avail = self.strings.len().saturating_sub(byte_offset).min(17);
            buf[..avail].copy_from_slice(&self.strings[byte_offset..byte_offset + avail]);
            let raw = u128::from_le_bytes(buf[..16].try_into().unwrap());
            let shifted = if bit_shift > 0 {
                let extra = buf[16] as u128;
                (raw >> bit_shift) | (extra << (128 - bit_shift))
            } else {
                raw
            };
            let mask = if needed_bits >= 128 { u128::MAX } else { (1u128 << needed_bits) - 1 };
            Kmer::<K>::new(<Kmer<K> as KmerBits>::from_u128(shifted & mask))
        }
    }

    /// Decode the m-mer (m <= 31) at an absolute base position.
    ///
    /// Used by the lookup path to verify that a bucket's first offset really
    /// holds the queried minimizer (see `Dictionary::lookup_with_minimizer`).
    /// Same branchless unaligned read as the K<=31 arm of
    /// [`Self::decode_kmer_at`]; the 8-byte tail padding makes it safe.
    #[inline(always)]
    pub fn decode_mmer_at(&self, absolute_pos: usize, m: usize) -> u64 {
        debug_assert!(m >= 1 && 2 * m < 64);
        let byte_offset = absolute_pos / 4;
        let bit_shift = (absolute_pos % 4) * 2;
        let (raw, extra) = unsafe {
            let base = self.strings.as_ptr().add(byte_offset);
            // SAFETY: as in decode_kmer_at — `strings` carries 8 bytes of
            // tail padding and a decodable m-mer starts before the logical end.
            (std::ptr::read_unaligned(base as *const u64), *base.add(8) as u64)
        };
        let result = (raw >> bit_shift) | ((extra << 1) << (63 - bit_shift));
        result & ((1u64 << (2 * m)) - 1)
    }

    /// Decode a k-mer at an absolute base position in the concatenated strings.
    #[inline(always)]
    pub fn decode_kmer_at<const K: usize>(&self, absolute_pos: usize) -> Kmer<K>
    where
        Kmer<K>: KmerBits,
    {
        let byte_offset = absolute_pos / 4;
        let bit_shift = (absolute_pos % 4) * 2;
        let needed_bits = K * 2;

        if K <= 31 {
            // K<=31: need at most 62+6=68 bits, so a u64 at `byte_offset` plus the
            // one byte the top bits can spill into.
            //
            // The spill is folded in unconditionally. `bit_shift` is `(pos % 4) * 2`,
            // so it is one of {0, 2, 4, 6} and the spill only actually matters for
            // {4, 6} — i.e. a branch that is taken for half of all positions, on a
            // value derived from the read position, so it does not predict. Making
            // it branchless is worth more than the byte load it saves.
            //
            // `(extra << 1) << (63 - bit_shift)` equals the wanted
            // `extra << (64 - bit_shift)` for bit_shift in {2, 4, 6}, and is
            // identically zero for bit_shift == 0 (the `<< 1` clears the only bit
            // that would survive `<< 63`). That dodges the shift-by-64 UB that
            // makes the direct expression unusable without the branch.
            let (raw, extra) = unsafe {
                let base = self.strings.as_ptr().add(byte_offset);
                // SAFETY: `strings` carries 8 bytes of tail padding (see `new` and
                // the deserializer), and a decodable k-mer starts before the logical
                // end, so `byte_offset + 8` is always within the allocation. This is
                // the same maximum address the previous branching form could touch.
                (std::ptr::read_unaligned(base as *const u64), *base.add(8) as u64)
            };
            let result = (raw >> bit_shift) | ((extra << 1) << (63 - bit_shift));
            let mask = (1u64 << needed_bits) - 1;
            Kmer::<K>::new(<Kmer<K> as KmerBits>::from_u64(result & mask))
        } else {
            let needed_bytes = (needed_bits + bit_shift).div_ceil(8);
            if needed_bytes <= 8 {
                let raw = if byte_offset + 8 <= self.strings.len() {
                    unsafe {
                        std::ptr::read_unaligned(
                            self.strings.as_ptr().add(byte_offset) as *const u64
                        )
                    }
                } else {
                    let mut buf = [0u8; 8];
                    let avail = self.strings.len() - byte_offset;
                    buf[..avail].copy_from_slice(&self.strings[byte_offset..byte_offset + avail]);
                    u64::from_le_bytes(buf)
                };
                let shifted = raw >> bit_shift;
                let mask = if needed_bits >= 64 { u64::MAX } else { (1u64 << needed_bits) - 1 };
                Kmer::<K>::new(<Kmer<K> as KmerBits>::from_u64(shifted & mask))
            } else {
                // For K=63 (needed_bits=126), when bit_shift > 2 we need more
                // than 128 bits from the byte stream. Load u128 + extra byte.
                let (raw, extra_byte) = if byte_offset + 17 <= self.strings.len() {
                    let r = unsafe {
                        std::ptr::read_unaligned(
                            self.strings.as_ptr().add(byte_offset) as *const u128
                        )
                    };
                    (r, self.strings[byte_offset + 16])
                } else if byte_offset + 16 <= self.strings.len() {
                    let r = unsafe {
                        std::ptr::read_unaligned(
                            self.strings.as_ptr().add(byte_offset) as *const u128
                        )
                    };
                    let extra = if byte_offset + 16 < self.strings.len() {
                        self.strings[byte_offset + 16]
                    } else {
                        0u8
                    };
                    (r, extra)
                } else {
                    let mut buf = [0u8; 17];
                    let avail = self.strings.len() - byte_offset;
                    buf[..avail].copy_from_slice(&self.strings[byte_offset..byte_offset + avail]);
                    (u128::from_le_bytes(buf[..16].try_into().unwrap()), buf[16])
                };
                let shifted = if bit_shift > 0 {
                    (raw >> bit_shift) | ((extra_byte as u128) << (128 - bit_shift))
                } else {
                    raw
                };
                let mask = if needed_bits >= 128 { u128::MAX } else { (1u128 << needed_bits) - 1 };
                Kmer::<K>::new(<Kmer<K> as KmerBits>::from_u128(shifted & mask))
            }
        }
    }

    /// Serialize the SPSS to a writer using a custom binary format.
    ///
    /// Format:
    /// - k: u64 (LE)
    /// - m: u64 (LE)
    /// - rel_offset_bits: u64 (LE)
    /// - strings_len: u64 (LE)
    /// - strings: [u8; strings_len]
    /// - begins: epserde BitFieldVec format
    /// - ends: Elias-Fano binary format
    /// - cum: Elias-Fano binary format
    pub fn serialize_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        use epserde::ser::Serialize;

        writer.write_all(&(self.k as u64).to_le_bytes())?;
        writer.write_all(&(self.m as u64).to_le_bytes())?;
        writer.write_all(&(self.rel_offset_bits as u64).to_le_bytes())?;
        // Write logical length (exclude 8-byte padding)
        let logical_len = self.strings.len() - 8;
        writer.write_all(&(logical_len as u64).to_le_bytes())?;
        writer.write_all(&self.strings[..logical_len])?;
        // SAFETY: BitFieldVec serialization preserves data layout
        unsafe {
            self.begins.serialize(writer).map_err(io::Error::other)?;
        }
        self.ends.write_to(writer)?;
        self.cum.write_to(writer)?;
        Ok(())
    }

    /// Deserialize an SPSS from a reader.
    pub fn deserialize_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        use epserde::deser::Deserialize;

        let mut buf8 = [0u8; 8];
        reader.read_exact(&mut buf8)?;
        let k = u64::from_le_bytes(buf8) as usize;
        reader.read_exact(&mut buf8)?;
        let m = u64::from_le_bytes(buf8) as usize;
        reader.read_exact(&mut buf8)?;
        let rel_offset_bits = u64::from_le_bytes(buf8) as u32;
        reader.read_exact(&mut buf8)?;
        let strings_len = u64::from_le_bytes(buf8) as usize;
        let mut strings = vec![0u8; strings_len + 8]; // +8 padding for safe unaligned reads
        reader.read_exact(&mut strings[..strings_len])?;
        // SAFETY: Reading data serialized by our own serialize method
        let begins = unsafe { BitFieldVec::deserialize_full(reader) }.map_err(io::Error::other)?;
        let ends = EliasFanoOffsets::read_from(reader)?;
        let cum = EliasFanoOffsets::read_from(reader)?;
        Ok(Self { strings, begins, ends, cum, rel_offset_bits, k, m })
    }
}

impl std::fmt::Debug for SpectrumPreservingStringSet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SpectrumPreservingStringSet")
            .field("k", &self.k)
            .field("m", &self.m)
            .field("num_strings", &self.num_strings())
            .field("total_bases", &self.total_bases())
            .field("num_kmers", &self.num_kmers())
            .field("num_bits", &self.num_bits())
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::offsets::OffsetsVector;

    /// The `K <= 31` arm of `decode_kmer_at` folds the spill byte in without a
    /// branch. Pin that it computes exactly what the branching form did, over
    /// every `(bit_shift, K)` combination that can actually occur — `bit_shift`
    /// is `(pos % 4) * 2`, and K is odd in `3..=31` on this arm.
    #[test]
    fn branchless_spill_matches_branching_form() {
        fn branching(raw: u64, extra: u8, bit_shift: usize, needed_bits: usize) -> u64 {
            let mut result = raw >> bit_shift;
            if bit_shift > 2 {
                result |= (extra as u64) << (64 - bit_shift);
            }
            result & ((1u64 << needed_bits) - 1)
        }
        fn branchless(raw: u64, extra: u8, bit_shift: usize, needed_bits: usize) -> u64 {
            let result = (raw >> bit_shift) | (((extra as u64) << 1) << (63 - bit_shift));
            result & ((1u64 << needed_bits) - 1)
        }

        let mut s: u64 = 0x243F_6A88_85A3_08D3;
        for _ in 0..20_000 {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            let extra = (s >> 11) as u8;
            for bit_shift in [0usize, 2, 4, 6] {
                for k in (3..=31).step_by(2) {
                    assert_eq!(
                        branching(s, extra, bit_shift, k * 2),
                        branchless(s, extra, bit_shift, k * 2),
                        "raw={s:#x} extra={extra:#x} bit_shift={bit_shift} k={k}"
                    );
                }
            }
        }
    }

    /// Helper: 2-bit pack a text (case-insensitive).
    fn pack_text(text: &str) -> Vec<u8> {
        let mut packed = Vec::new();
        for (i, &b) in text.as_bytes().iter().enumerate() {
            let bits = match b {
                b'A' | b'a' => 0u8,
                b'C' | b'c' => 1u8,
                b'G' | b'g' => 3u8,
                b'T' | b't' => 2u8,
                _ => panic!("invalid base"),
            };
            let byte_idx = i / 4;
            let bit_off = (i % 4) * 2;
            if byte_idx >= packed.len() {
                packed.push(0u8);
            }
            packed[byte_idx] |= bits << bit_off;
        }
        packed
    }

    /// Helper: build a test SPSS from DNA strings (concatenated, partition).
    fn build_test_spss(k: usize, m: usize, strings: &[&str]) -> SpectrumPreservingStringSet {
        let text: String = strings.concat();
        let packed = pack_text(&text);
        let mut offsets = OffsetsVector::new();
        let mut total_bases: u64 = 0;
        for s in strings {
            total_bases += s.len() as u64;
            offsets.push(total_bases);
        }
        SpectrumPreservingStringSet::from_parts(packed, offsets, k, m)
    }

    #[test]
    fn test_spss_creation() {
        let spss = SpectrumPreservingStringSet::new(31, 13);
        assert_eq!(spss.k(), 31);
        assert_eq!(spss.m(), 13);
        assert_eq!(spss.num_strings(), 0);
    }

    #[test]
    fn test_spss_two_strings() {
        let spss = build_test_spss(31, 13, &[
            "ACGTACGTACGTACGTACGTACGTACGTACG",  // 31 chars
            "TGCATGCATGCATGCATGCATGCATGCATGCA", // 32 chars
        ]);
        assert_eq!(spss.num_strings(), 2);
        assert_eq!(spss.num_kmers(), 1 + 2);
    }

    #[test]
    fn test_spss_string_offsets() {
        let spss = build_test_spss(31, 13, &[
            "ACGTACGTACGTACGTACGTACGTACGTACG",  // 31 chars
            "TGCATGCATGCATGCATGCATGCATGCATGC", // 31 chars
        ]);

        let (begin1, end1) = spss.string_offsets(0);
        let (begin2, end2) = spss.string_offsets(1);

        assert_eq!(begin1, 0);
        assert_eq!(end1 - begin1, 31);
        assert_eq!(begin2, 31);
        assert_eq!(end2 - begin2, 31);
        assert_eq!(spss.first_kmer_id(0), 0);
        assert_eq!(spss.first_kmer_id(1), 1);
    }

    #[test]
    fn test_spss_total_bases() {
        let spss = build_test_spss(31, 13, &[
            "ACGTACGTACGTACGTACGTACGTACGTACG",  // 31 chars
            "TGCATGCATGCATGCATGCATGCATGCATGC", // 31 chars
        ]);
        assert_eq!(spss.total_bases(), 62);
    }

    #[test]
    fn test_spss_debug() {
        let spss = build_test_spss(31, 13, &[
            "ACGTACGTACGTACGTACGTACGTACGTACG",
        ]);
        let debug_str = format!("{:?}", spss);
        assert!(debug_str.contains("SpectrumPreservingStringSet"));
        assert!(debug_str.contains("k: 31"));
    }

    #[test]
    fn test_spss_overlapping_intervals() {
        // Text of 20 bases; k=7, m=3; intervals [0, 10) and [7, 17):
        // overlap of 3 < k, 4 k-mers each.
        let text = "ACGTACGTACGTACGTACGT";
        let packed = pack_text(text);
        let spss = SpectrumPreservingStringSet::from_interval_vecs(
            packed,
            vec![0, 7],
            vec![10, 17],
            7,
            3,
        );
        assert_eq!(spss.num_strings(), 2);
        assert_eq!(spss.string_offsets(0), (0, 10));
        assert_eq!(spss.string_offsets(1), (7, 17));
        assert_eq!(spss.num_kmers(), 8);
        assert_eq!(spss.first_kmer_id(0), 0);
        assert_eq!(spss.first_kmer_id(1), 4);
        assert_eq!(spss.total_bases(), 17);

        // locate_kmer: each k-mer occurrence lies in at most one interval.
        assert_eq!(spss.locate_kmer(0, 7), Some((0, 0, 10)));
        assert_eq!(spss.locate_kmer(3, 7), Some((0, 0, 10)));
        assert_eq!(spss.locate_kmer(7, 7), Some((1, 7, 17)));
        assert_eq!(spss.locate_kmer(10, 7), Some((1, 7, 17)));
        // A k-mer at position 4 spans [4, 11): in neither interval.
        assert_eq!(spss.locate_kmer(4, 7), None);
        assert_eq!(spss.locate_kmer(5, 7), None);
        assert_eq!(spss.locate_kmer(6, 7), None);
        // Out of range.
        assert_eq!(spss.locate_kmer(11, 7), None);

        // Encoded offsets round-trip, distinguishing the overlap positions:
        // a minimizer occurrence at absolute position 7 can be recorded by a
        // super-kmer of either interval.
        let e0 = spss.encode_offset(7, 0, 0);
        let e1 = spss.encode_offset(7, 7, 1);
        assert_ne!(e0, e1);
        let d0 = spss.decode_offset(e0);
        let d1 = spss.decode_offset(e1);
        assert_eq!((d0.absolute, d0.string_id, d0.string_begin, d0.string_end), (7, 0, 0, 10));
        assert_eq!((d1.absolute, d1.string_id, d1.string_begin, d1.string_end), (7, 1, 7, 17));
    }

    #[test]
    fn test_spss_interval_serialization_roundtrip() {
        let text = "ACGTACGTACGTACGTACGT";
        let packed = pack_text(text);
        let spss = SpectrumPreservingStringSet::from_interval_vecs(
            packed,
            vec![0, 7],
            vec![10, 17],
            7,
            3,
        );

        let mut buf = Vec::new();
        spss.serialize_to(&mut buf).unwrap();
        let spss2 = SpectrumPreservingStringSet::deserialize_from(&mut &buf[..]).unwrap();

        assert_eq!(spss2.k(), 7);
        assert_eq!(spss2.m(), 3);
        assert_eq!(spss2.num_strings(), 2);
        assert_eq!(spss2.num_kmers(), 8);
        assert_eq!(spss2.string_offsets(0), (0, 10));
        assert_eq!(spss2.string_offsets(1), (7, 17));
        assert_eq!(spss2.rel_offset_bits(), spss.rel_offset_bits());
        let d = spss2.decode_offset(spss2.encode_offset(9, 7, 1));
        assert_eq!((d.absolute, d.string_id), (9, 1));
        for pos in 0..17usize {
            assert_eq!(
                spss.decode_mmer_at(pos, 3),
                spss2.decode_mmer_at(pos, 3)
            );
        }
    }
}
