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

use crate::cli::Cli;
use crate::kmer_comp;
use crate::seq_parse;

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
