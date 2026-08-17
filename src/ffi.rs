//! C ABI for in-process SNPmer detection.
//!
//! Lets a C/C++ host (e.g. the hifiasm candidate-overlap fork) obtain the
//! gated SNPmer set directly, without the `snpmers.tsv` round-trip. The set is
//! exactly what `main` would dump: pairs of allele k-mers that passed the
//! both-alleles + binomial + Fisher strand-balance + high-frequency gates in
//! `get_snpmers_inplace_sort`.
//!
//! Bit order: each k-mer is a 2-bit-packed `u64` in **myloasm order**, i.e.
//! base at read position `i` occupies bits `[2i, 2i+1]` (first base in the low
//! bits), with A=0,C=1,G=2,T=3 (see `types::BYTE_TO_SEQ`). This is the same
//! packing `main`'s TSV dump decodes. hifiasm rolls its k-mers in the opposite
//! order (first base in the high bits), so the host must re-pack on ingest;
//! this is documented on the C side.

use std::ffi::CStr;
use std::os::raw::{c_char, c_int};
use std::panic::{self, AssertUnwindSafe};

use clap::Parser;

use crate::chain;
use crate::cli::Cli;
use crate::kmer_comp;
use crate::overlap;
use crate::seq_parse;
use crate::types::KmerGlobalInfo;

/// One SNPmer site: the two allele k-mers (myloasm bit order) and their
/// strand-summed counts. Layout is `#[repr(C)]` and must match the C struct.
#[repr(C)]
pub struct MyloSnpmer {
    pub allele0_kmer: u64,
    pub allele1_kmer: u64,
    pub allele0_count: u32,
    pub allele1_count: u32,
}

/// Result handle returned to C. Owns the allocation; free with
/// `myloasm_snpmers_free`. `ptr`/`len` describe a contiguous `MyloSnpmer`
/// array; `capacity` is retained so the exact `Vec` can be reconstructed.
#[repr(C)]
pub struct MyloSnpmerSet {
    pub ptr: *mut MyloSnpmer,
    pub len: usize,
    pub capacity: usize,
    pub k: c_int,
}

/// Detect SNPmers from the given read files.
///
/// - `paths`: array of `n_paths` NUL-terminated UTF-8 FASTA/FASTQ(.gz) paths.
/// - `kmer_size`: k (odd, < 24); pass 0 to use myloasm's default (21).
/// - `threads`: worker threads; pass 0 to use myloasm's default.
/// - `out`: on success, populated with the result set.
///
/// Returns 0 on success, non-zero on error (invalid args, panic, or no reads).
/// On error `out` is zeroed. Any Rust panic is caught and converted to an error
/// code so it never unwinds across the C boundary.
///
/// # Safety
/// `paths` must point to `n_paths` valid C strings; `out` must be a valid,
/// writable pointer.
#[no_mangle]
pub unsafe extern "C" fn myloasm_detect_snpmers(
    paths: *const *const c_char,
    n_paths: usize,
    kmer_size: c_int,
    threads: c_int,
    out: *mut MyloSnpmerSet,
) -> c_int {
    if out.is_null() {
        return 1;
    }
    // Zero the output up front so early returns leave it well-defined.
    std::ptr::write(
        out,
        MyloSnpmerSet {
            ptr: std::ptr::null_mut(),
            len: 0,
            capacity: 0,
            k: 0,
        },
    );
    if paths.is_null() || n_paths == 0 {
        return 2;
    }

    let result = panic::catch_unwind(AssertUnwindSafe(|| {
        // Collect input paths.
        let mut files: Vec<String> = Vec::with_capacity(n_paths);
        for i in 0..n_paths {
            let p = *paths.add(i);
            if p.is_null() {
                return Err(3);
            }
            match CStr::from_ptr(p).to_str() {
                Ok(s) => files.push(s.to_owned()),
                Err(_) => return Err(4),
            }
        }

        // Build a Cli with all defaults via clap, then override the few fields
        // detection reads. parse_from guarantees every default is populated, so
        // we never depend on hand-initialising the large Cli struct.
        let mut argv: Vec<String> = vec!["myloasm".to_string()];
        argv.extend(files.iter().cloned());
        let mut args = Cli::parse_from(argv);

        if kmer_size > 0 {
            args.kmer_size = kmer_size as usize;
        }
        if args.kmer_size % 2 == 0 {
            return Err(5); // k must be odd
        }
        if threads > 0 {
            args.threads = threads as usize;
        }

        // Count all k-mers, then detect + gate SNPmers. Neither call writes to
        // disk. rayon's global pool may already be initialised by a prior call;
        // ignore the resulting error rather than propagating it.
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .stack_size(16 * 1024 * 1024)
            .build_global();

        let big_kmer_map =
            seq_parse::read_to_split_kmers(args.kmer_size, args.threads, &args);
        let info =
            kmer_comp::get_snpmers_inplace_sort(big_kmer_map, args.kmer_size, &args);

        let k = args.kmer_size;
        let mut v: Vec<MyloSnpmer> = Vec::with_capacity(info.snpmer_info.len());
        for s in info.snpmer_info.iter() {
            // Reconstruct the two allele k-mers exactly as the TSV dump does:
            // the split k-mer with the middle base OR'd into the top position.
            let kmer0 = s.split_kmer | ((s.mid_bases[0] as u64) << (k - 1));
            let kmer1 = s.split_kmer | ((s.mid_bases[1] as u64) << (k - 1));
            v.push(MyloSnpmer {
                allele0_kmer: kmer0,
                allele1_kmer: kmer1,
                allele0_count: s.counts[0],
                allele1_count: s.counts[1],
            });
        }
        Ok((v, k))
    }));

    match result {
        Ok(Ok((mut v, k))) => {
            let len = v.len();
            let capacity = v.capacity();
            let ptr = v.as_mut_ptr();
            std::mem::forget(v); // ownership transferred to C
            std::ptr::write(
                out,
                MyloSnpmerSet {
                    ptr,
                    len,
                    capacity,
                    k: k as c_int,
                },
            );
            0
        }
        Ok(Err(code)) => code,
        Err(_) => 100, // caught panic
    }
}

/// Free a set returned by `myloasm_detect_snpmers`. Safe to call with a zeroed
/// set (null ptr). Must be called exactly once per successful detect.
///
/// # Safety
/// `set` must point to a `MyloSnpmerSet` previously written by
/// `myloasm_detect_snpmers` and not already freed.
#[no_mangle]
pub unsafe extern "C" fn myloasm_snpmers_free(set: *mut MyloSnpmerSet) {
    if set.is_null() {
        return;
    }
    let s = &mut *set;
    if !s.ptr.is_null() {
        drop(Vec::from_raw_parts(s.ptr, s.len, s.capacity));
        s.ptr = std::ptr::null_mut();
        s.len = 0;
        s.capacity = 0;
    }
}

// ---------------------------------------------------------------------------
// Per-read marker indexing (open syncmers + SNPmers) for the hifiasm
// fake-chain second pass.
// ---------------------------------------------------------------------------
//
// This exposes myloasm's own read indexing (get_twin_read_syncmer via
// twin_reads_from_snpmers): for every read it returns the open-syncmer
// positions and the SNPmer positions, each with the read's CANONICAL k-mer key
// at that position. The keys are the values Kmer48::to_u64() produces
// (canonical = min of fwd / rev-comp under the middle-masked comparison), so
// the SAME genomic locus yields the SAME key on any read regardless of strand.
// A C consumer can therefore match markers between two reads purely by key
// equality, with no bit-order or reverse-complement handling.
//
// Marker key layout (k=21): base i occupies bits [2i, 2i+1]; the middle base
// (index 10) is at bits [20,21]. Clipping to k=20 by masking the low 40 bits
// keeps the middle base, so SNPmer alleles stay distinct after the clip.

/// One marker on a read: its start position (0-based, forward strand) and the
/// canonical k-mer key at that position. `#[repr(C)]`; must match the C struct.
#[repr(C)]
pub struct MyloMarker {
    pub pos: u32,
    pub _pad: u32,
    pub key: u64,
}

/// Markers for one read. `syncmers` / `snpmers` point into the shared arenas
/// owned by `MyloReadIndex` (do not free individually). `name` points into the
/// shared name arena and is `name_len` bytes, NOT NUL-terminated; it is the
/// read's base id (first whitespace-delimited token), matching hifiasm's PAF
/// read names. `#[repr(C)]`.
#[repr(C)]
pub struct MyloReadMarkers {
    pub name: *const c_char,
    pub name_len: usize,
    pub syncmers: *const MyloMarker,
    pub n_syncmers: usize,
    pub snpmers: *const MyloMarker,
    pub n_snpmers: usize,
}

/// Result of `myloasm_index_reads`. Owns four contiguous allocations:
/// `reads` (the per-read descriptors), one syncmer arena and one snpmer arena
/// (all reads' markers concatenated; each read's slice is referenced by the
/// descriptor), and one name-byte arena. Free with `myloasm_read_index_free`.
/// `k` is the marker k-mer length used (21 by default). Capacities are retained
/// so the exact Vecs can be reconstructed on free. `#[repr(C)]`.
#[repr(C)]
pub struct MyloReadIndex {
    pub reads: *mut MyloReadMarkers,
    pub n_reads: usize,
    reads_cap: usize,

    sync_arena: *mut MyloMarker,
    sync_len: usize,
    sync_cap: usize,

    snp_arena: *mut MyloMarker,
    snp_len: usize,
    snp_cap: usize,

    name_arena: *mut u8,
    name_len: usize,
    name_cap: usize,

    pub k: c_int,
}

fn zeroed_read_index() -> MyloReadIndex {
    MyloReadIndex {
        reads: std::ptr::null_mut(),
        n_reads: 0,
        reads_cap: 0,
        sync_arena: std::ptr::null_mut(),
        sync_len: 0,
        sync_cap: 0,
        snp_arena: std::ptr::null_mut(),
        snp_len: 0,
        snp_cap: 0,
        name_arena: std::ptr::null_mut(),
        name_len: 0,
        name_cap: 0,
        k: 0,
    }
}

/// Detect SNPmers and index every read with open syncmers + SNPmers.
///
/// - `paths`: array of `n_paths` NUL-terminated FASTA/FASTQ(.gz) paths.
/// - `kmer_size`: k (odd, < 24); pass 0 for myloasm's default (21).
/// - `c`: syncmer compression (paper's c; s = k - min(c,11) + 1); 0 for default.
/// - `threads`: worker threads; 0 for default.
/// - `out`: on success, populated with the per-read marker index.
///
/// Returns 0 on success, non-zero on error. On error `out` is zeroed. Panics are
/// caught and turned into an error code so nothing unwinds across the boundary.
///
/// # Safety
/// `paths` must point to `n_paths` valid C strings; `out` must be a valid,
/// writable pointer. The returned index must be freed with
/// `myloasm_read_index_free`.
#[no_mangle]
pub unsafe extern "C" fn myloasm_index_reads(
    paths: *const *const c_char,
    n_paths: usize,
    kmer_size: c_int,
    c: c_int,
    threads: c_int,
    out: *mut MyloReadIndex,
) -> c_int {
    if out.is_null() {
        return 1;
    }
    std::ptr::write(out, zeroed_read_index());
    if paths.is_null() || n_paths == 0 {
        return 2;
    }

    let result = panic::catch_unwind(AssertUnwindSafe(|| {
        let mut files: Vec<String> = Vec::with_capacity(n_paths);
        for i in 0..n_paths {
            let p = *paths.add(i);
            if p.is_null() {
                return Err(3);
            }
            match CStr::from_ptr(p).to_str() {
                Ok(s) => files.push(s.to_owned()),
                Err(_) => return Err(4),
            }
        }

        let mut argv: Vec<String> = vec!["myloasm".to_string()];
        argv.extend(files.iter().cloned());
        let mut args = Cli::parse_from(argv);
        if kmer_size > 0 {
            args.kmer_size = kmer_size as usize;
        }
        if args.kmer_size % 2 == 0 {
            return Err(5); // k must be odd
        }
        if c > 0 {
            args.c = c as usize;
        }
        if threads > 0 {
            args.threads = threads as usize;
        }

        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .stack_size(16 * 1024 * 1024)
            .build_global();

        // Detect SNPmers (same gates as myloasm_detect_snpmers), then index
        // every read with syncmers + SNPmers. twin_reads_from_snpmers re-reads
        // the files and returns reads in worker-completion order, so we carry
        // each read's base id back to the host for name-based id mapping.
        let big = seq_parse::read_to_split_kmers(args.kmer_size, args.threads, &args);
        let mut info: KmerGlobalInfo =
            kmer_comp::get_snpmers_inplace_sort(big, args.kmer_size, &args);
        let twins = kmer_comp::twin_reads_from_snpmers(&mut info, &args);

        let k = args.kmer_size;

        // Flatten into shared arenas; each read records (offset,len) slices.
        let mut reads: Vec<MyloReadMarkers> = Vec::with_capacity(twins.len());
        let mut sync_arena: Vec<MyloMarker> = Vec::new();
        let mut snp_arena: Vec<MyloMarker> = Vec::new();
        let mut name_arena: Vec<u8> = Vec::new();

        // Record (offset,len) first; resolve to pointers after arenas stop
        // growing (growth reallocates and would dangle earlier pointers).
        struct Slice {
            name_off: usize,
            name_len: usize,
            sync_off: usize,
            n_sync: usize,
            snp_off: usize,
            n_snp: usize,
        }
        let mut slices: Vec<Slice> = Vec::with_capacity(twins.len());

        for tr in twins.iter() {
            let name_off = name_arena.len();
            name_arena.extend_from_slice(tr.base_id.as_bytes());
            let name_len = tr.base_id.len();

            let sync_off = sync_arena.len();
            for (pos, key) in tr.minimizers_vec() {
                sync_arena.push(MyloMarker {
                    pos,
                    _pad: 0,
                    key: key.to_u64(),
                });
            }
            let n_sync = sync_arena.len() - sync_off;

            let snp_off = snp_arena.len();
            for (pos, key) in tr.snpmers_vec() {
                snp_arena.push(MyloMarker {
                    pos,
                    _pad: 0,
                    key: key.to_u64(),
                });
            }
            let n_snp = snp_arena.len() - snp_off;

            slices.push(Slice {
                name_off,
                name_len,
                sync_off,
                n_sync,
                snp_off,
                n_snp,
            });
        }

        // Arenas are final; take stable base pointers.
        let sync_base = sync_arena.as_ptr();
        let snp_base = snp_arena.as_ptr();
        let name_base = name_arena.as_ptr();
        for s in slices.iter() {
            reads.push(MyloReadMarkers {
                name: name_base.add(s.name_off) as *const c_char,
                name_len: s.name_len,
                syncmers: sync_base.add(s.sync_off),
                n_syncmers: s.n_sync,
                snpmers: snp_base.add(s.snp_off),
                n_snpmers: s.n_snp,
            });
        }

        Ok((reads, sync_arena, snp_arena, name_arena, k))
    }));

    match result {
        Ok(Ok((mut reads, mut sync_arena, mut snp_arena, mut name_arena, k))) => {
            let idx = MyloReadIndex {
                reads: reads.as_mut_ptr(),
                n_reads: reads.len(),
                reads_cap: reads.capacity(),
                sync_arena: sync_arena.as_mut_ptr(),
                sync_len: sync_arena.len(),
                sync_cap: sync_arena.capacity(),
                snp_arena: snp_arena.as_mut_ptr(),
                snp_len: snp_arena.len(),
                snp_cap: snp_arena.capacity(),
                name_arena: name_arena.as_mut_ptr(),
                name_len: name_arena.len(),
                name_cap: name_arena.capacity(),
                k: k as c_int,
            };
            std::mem::forget(reads);
            std::mem::forget(sync_arena);
            std::mem::forget(snp_arena);
            std::mem::forget(name_arena);
            std::ptr::write(out, idx);
            0
        }
        Ok(Err(code)) => code,
        Err(_) => 100,
    }
}

// ---------------------------------------------------------------------------
// Anchor chaining (myloasm's dp_anchors_v2) for the hifiasm fake-chain pass.
// ---------------------------------------------------------------------------
//
// The hifiasm bridge matches markers in C++ (syncmers + SNPmers, canonical
// keys, m:n within the candidate interval), then calls this to run myloasm's
// own chaining DP over the merged anchor set. Only the chainer is exposed; the
// caller owns matching and anchor construction.
//
// Strand: encode each anchor's query position with the relative strand in bit
// 31 (0 = forward/same-strand, 1 = reverse-complement); `pos2` is the raw
// target position. Within one candidate the strand is uniform (hifiasm already
// decided it), so all anchors share the same bit.

/// One input anchor for chaining. `slot` is an opaque caller index carried
/// through to the output so the host can recover the anchor's type tag
/// (syncmer/SNPmer) and identity. `pos1` = query position with strand in bit 31;
/// `pos2` = raw target position. `#[repr(C)]`; must match the C struct.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct MyloAnchor {
    pub pos1: u32,
    pub pos2: u32,
    pub slot: u32,
    pub _pad: u32,
}

/// One anchor on the returned chain: the input `slot`, plus the (strand-stripped)
/// query/target positions. Ordered along the chain (increasing query position).
#[repr(C)]
#[derive(Clone, Copy)]
pub struct MyloChainAnchor {
    pub slot: u32,
    pub qpos: u32,
    pub tpos: u32,
    pub _pad: u32,
}

/// Chain a set of anchors with myloasm's DP and write the single best chain into
/// a caller-provided output buffer.
///
/// - `anchors` / `n_anchors`: input anchors (need not be sorted).
/// - `gap_cost`: c (pass 0 for the default 11).
/// - `match_score`: per-anchor score (pass 0 for the default 1).
/// - `band`: predecessor-iteration bound = max_mult * 20; pass 0 to let this
///   function use n_anchors (unbounded within the set).
/// - `min_chain_length`: pass 0 for the default 3.
/// - `out` / `out_cap`: caller buffer for the chained anchors (cap should be >=
///   n_anchors; the chain can never exceed the input count).
/// - `out_n`: receives the number of anchors written.
/// - `out_score`: receives the chain score.
/// - `out_is_reverse`: receives 1 if the best chain is reverse-strand, else 0.
///
/// Returns 0 on success (including "no chain": *out_n = 0), non-zero on error.
/// Panics are caught. Pure/thread-safe: no global state.
///
/// # Safety
/// `anchors` must point to `n_anchors` valid `MyloAnchor`; `out` must point to
/// `out_cap` writable `MyloChainAnchor`; the out scalars must be valid pointers.
#[no_mangle]
pub unsafe extern "C" fn myloasm_chain(
    anchors: *const MyloAnchor,
    n_anchors: usize,
    gap_cost: c_int,
    match_score: c_int,
    band: usize,
    min_chain_length: usize,
    out: *mut MyloChainAnchor,
    out_cap: usize,
    out_n: *mut usize,
    out_score: *mut c_int,
    out_is_reverse: *mut c_int,
) -> c_int {
    if out_n.is_null() || out_score.is_null() || out_is_reverse.is_null() {
        return 1;
    }
    *out_n = 0;
    *out_score = 0;
    *out_is_reverse = 0;
    if n_anchors == 0 {
        return 0;
    }
    if anchors.is_null() || (out.is_null() && out_cap != 0) {
        return 2;
    }

    let result = panic::catch_unwind(AssertUnwindSafe(|| {
        let input = std::slice::from_raw_parts(anchors, n_anchors);
        let matches: Vec<chain::Anchor> = input
            .iter()
            .map(|a| chain::Anchor {
                i: Some(a.slot),
                j: None,
                pos1: a.pos1,
                pos2: a.pos2,
            })
            .collect();

        let gc = if gap_cost > 0 { gap_cost } else { 11 };
        let ms = if match_score > 0 { match_score } else { 1 };
        let bnd = if band > 0 { band } else { n_anchors };
        let mcl = if min_chain_length > 0 {
            min_chain_length
        } else {
            chain::CHAIN_MIN_LENGTH
        };

        chain::best_chain(
            &matches,
            gc,
            ms,
            bnd,
            chain::CHAIN_MAX_GAP,
            chain::CHAIN_MAX_SKIP,
            chain::CHAIN_DOUBLE_GAP,
            mcl,
        )
    }));

    match result {
        Ok(Some((score, chain_anchors, is_reverse))) => {
            let n = chain_anchors.len().min(out_cap);
            for (k, a) in chain_anchors.iter().take(n).enumerate() {
                // pos1 here already has the strand bit stripped (best_chain masks
                // it); pos2 is the raw target position.
                *out.add(k) = MyloChainAnchor {
                    slot: a.i.unwrap_or(0),
                    qpos: a.pos1,
                    tpos: a.pos2,
                    _pad: 0,
                };
            }
            *out_n = n;
            *out_score = score;
            *out_is_reverse = if is_reverse { 1 } else { 0 };
            0
        }
        Ok(None) => 0, // no chain of >= min length
        Err(_) => 100,
    }
}

/// Free an index returned by `myloasm_index_reads`. Safe on a zeroed index.
///
/// # Safety
/// `idx` must point to a `MyloReadIndex` previously written by
/// `myloasm_index_reads` and not already freed.
#[no_mangle]
pub unsafe extern "C" fn myloasm_read_index_free(idx: *mut MyloReadIndex) {
    if idx.is_null() {
        return;
    }
    let s = &mut *idx;
    if !s.reads.is_null() {
        drop(Vec::from_raw_parts(s.reads, s.n_reads, s.reads_cap));
    }
    if !s.sync_arena.is_null() {
        drop(Vec::from_raw_parts(s.sync_arena, s.sync_len, s.sync_cap));
    }
    if !s.snp_arena.is_null() {
        drop(Vec::from_raw_parts(s.snp_arena, s.snp_len, s.snp_cap));
    }
    if !s.name_arena.is_null() {
        drop(Vec::from_raw_parts(s.name_arena, s.name_len, s.name_cap));
    }
    std::ptr::write(idx, zeroed_read_index());
}

// ---------------------------------------------------------------------------
// Native all-vs-all overlap detection.
//
// Runs myloasm's own overlapper (restored in src/overlap.rs from the pre-strip
// assembler): index every read's syncmers, find candidate pairs by shared
// minimizers, chain with the shared DP, refine with SNPmers, and emit dovetail
// + containment overlaps. This is the myloasm-native replacement for the
// hifiasm candidate + fakechain path.
// ---------------------------------------------------------------------------

/// One overlap between two reads. Coordinates are forward-strand, half-open, on
/// each read's own sequence. `reverse` = read2 is reverse-complemented relative
/// to read1. `#[repr(C)]`.
#[repr(C)]
pub struct MyloOverlap {
    pub read_i: u32,
    pub read_j: u32,
    pub start1: u32,
    pub end1: u32,
    pub start2: u32,
    pub end2: u32,
    pub shared_minimizers: u32,
    pub shared_snpmers: u32,
    pub diff_snpmers: u32,
    /// 1 = read_j is reverse-complemented w.r.t. read_i, 0 = same strand.
    pub reverse: u8,
    /// 1 = one read is contained in the other, 0 = dovetail.
    pub contained: u8,
    pub _pad: [u8; 2],
}

/// Result of `myloasm_detect_overlaps`. `overlaps`/`n_overlaps` is the overlap
/// array; `names`/`name_offsets`/`n_reads` map each read index (read_i/read_j)
/// to its base id (name_offsets has n_reads+1 entries, byte ranges into names).
/// Free with `myloasm_overlaps_free`. `#[repr(C)]`.
#[repr(C)]
pub struct MyloOverlaps {
    pub overlaps: *mut MyloOverlap,
    pub n_overlaps: usize,
    overlaps_cap: usize,

    pub names: *mut c_char,
    pub name_offsets: *mut usize,
    pub n_reads: usize,
    names_cap: usize,
    name_offsets_cap: usize,
}

fn zeroed_overlaps() -> MyloOverlaps {
    MyloOverlaps {
        overlaps: std::ptr::null_mut(),
        n_overlaps: 0,
        overlaps_cap: 0,
        names: std::ptr::null_mut(),
        name_offsets: std::ptr::null_mut(),
        n_reads: 0,
        names_cap: 0,
        name_offsets_cap: 0,
    }
}

#[no_mangle]
pub unsafe extern "C" fn myloasm_detect_overlaps(
    paths: *const *const c_char,
    n_paths: usize,
    kmer_size: c_int,
    c: c_int,
    threads: c_int,
    out: *mut MyloOverlaps,
) -> c_int {
    if out.is_null() {
        return 1;
    }
    std::ptr::write(out, zeroed_overlaps());
    if paths.is_null() || n_paths == 0 {
        return 2;
    }

    let result = panic::catch_unwind(AssertUnwindSafe(|| {
        let mut files: Vec<String> = Vec::with_capacity(n_paths);
        for i in 0..n_paths {
            let p = *paths.add(i);
            if p.is_null() {
                return Err(3);
            }
            match CStr::from_ptr(p).to_str() {
                Ok(s) => files.push(s.to_owned()),
                Err(_) => return Err(4),
            }
        }

        let mut argv: Vec<String> = vec!["myloasm".to_string()];
        argv.extend(files.iter().cloned());
        let mut args = Cli::parse_from(argv);
        if kmer_size > 0 {
            args.kmer_size = kmer_size as usize;
        }
        if args.kmer_size % 2 == 0 {
            return Err(5); // k must be odd
        }
        if c > 0 {
            args.c = c as usize;
        }
        if threads > 0 {
            args.threads = threads as usize;
        }

        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .stack_size(16 * 1024 * 1024)
            .build_global();

        // Same read -> twin-reads pipeline as myloasm_index_reads.
        let big = seq_parse::read_to_split_kmers(args.kmer_size, args.threads, &args);
        let mut info: KmerGlobalInfo =
            kmer_comp::get_snpmers_inplace_sort(big, args.kmer_size, &args);
        let twins = kmer_comp::twin_reads_from_snpmers(&mut info, &args);

        // All reads are "outer" reads: run all-vs-all over the full set.
        let outer: Vec<usize> = (0..twins.len()).collect();
        let configs = overlap::get_overlaps_outer_reads_twin(
            &twins, &outer, &args, None, None,
        );

        // Pack overlaps.
        let mut overlaps: Vec<MyloOverlap> = Vec::with_capacity(configs.len());
        for oc in configs.iter() {
            overlaps.push(MyloOverlap {
                read_i: oc.read_i as u32,
                read_j: oc.read_j as u32,
                start1: oc.start1 as u32,
                end1: oc.end1 as u32,
                start2: oc.start2 as u32,
                end2: oc.end2 as u32,
                shared_minimizers: oc.shared_mini as u32,
                shared_snpmers: oc.shared_snpmer as u32,
                diff_snpmers: oc.diff_snpmer as u32,
                reverse: if oc.reverse { 1 } else { 0 },
                contained: if oc.contained { 1 } else { 0 },
                _pad: [0; 2],
            });
        }

        // Pack read names (base id) so the host can map read_i/read_j -> dinara
        // reads, exactly like myloasm_index_reads.
        let mut name_arena: Vec<u8> = Vec::new();
        let mut name_offsets: Vec<usize> = Vec::with_capacity(twins.len() + 1);
        name_offsets.push(0);
        for tr in twins.iter() {
            name_arena.extend_from_slice(tr.base_id.as_bytes());
            name_offsets.push(name_arena.len());
        }

        Ok((overlaps, name_arena, name_offsets, twins.len()))
    }));

    match result {
        Ok(Ok((mut overlaps, mut name_arena, mut name_offsets, n_reads))) => {
            overlaps.shrink_to_fit();
            name_arena.shrink_to_fit();
            name_offsets.shrink_to_fit();

            let o = MyloOverlaps {
                overlaps: overlaps.as_mut_ptr(),
                n_overlaps: overlaps.len(),
                overlaps_cap: overlaps.capacity(),
                names: name_arena.as_mut_ptr() as *mut c_char,
                name_offsets: name_offsets.as_mut_ptr(),
                n_reads,
                names_cap: name_arena.capacity(),
                name_offsets_cap: name_offsets.capacity(),
            };
            std::mem::forget(overlaps);
            std::mem::forget(name_arena);
            std::mem::forget(name_offsets);
            std::ptr::write(out, o);
            0
        }
        Ok(Err(code)) => code,
        Err(_) => 100,
    }
}

#[no_mangle]
pub unsafe extern "C" fn myloasm_overlaps_free(ov: *mut MyloOverlaps) {
    if ov.is_null() {
        return;
    }
    let s = &mut *ov;
    if !s.overlaps.is_null() {
        drop(Vec::from_raw_parts(s.overlaps, s.n_overlaps, s.overlaps_cap));
    }
    if !s.names.is_null() {
        drop(Vec::from_raw_parts(
            s.names as *mut u8,
            s.name_offsets_last(),
            s.names_cap,
        ));
    }
    if !s.name_offsets.is_null() {
        // name_offsets has n_reads+1 entries.
        drop(Vec::from_raw_parts(
            s.name_offsets,
            s.n_reads + 1,
            s.name_offsets_cap,
        ));
    }
    std::ptr::write(ov, zeroed_overlaps());
}

impl MyloOverlaps {
    // Total name-arena byte length = last entry of name_offsets. Read before
    // name_offsets is freed.
    unsafe fn name_offsets_last(&self) -> usize {
        if self.name_offsets.is_null() || self.n_reads == 0 {
            0
        } else {
            *self.name_offsets.add(self.n_reads)
        }
    }
}
