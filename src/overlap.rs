//! All-vs-all read overlapper restored from the pre-strip myloasm assembler
//! (commit b450f39^: src/mapping.rs + src/twin_graph.rs). This is the native
//! myloasm overlap path -- index every read's syncmers, find candidate pairs by
//! shared minimizers, chain with the shared DP (chain::dp_anchors_v2), refine
//! with SNPmers, and emit dovetail/containment overlaps -- with no dependency on
//! the stripped graph/unitig/polishing modules.
//!
//! Deviations from the original:
//!   - Debug overlap/containment file dumps are sinked (flate2 was removed).
//!   - End-extension (polishing::alignment::extend_ends_chain) is not applied;
//!     overlap coordinates come straight from the minimizer chain endpoints.

#![allow(clippy::all)]

use crate::chain::dp_anchors_v2;
use crate::cli::Cli;
use crate::constants::{
    IDENTITY_THRESHOLDS, MAX_GAP_CHAINING, MAX_MULTIPLICITY_KMER, MIN_CHAIN_SCORE_COMPARE,
    OVERLAP_HANG_LENGTH,
};
// Anchor comes from types (same struct ChainInfo/TwinOverlap/CompareTwinReadOptions
// use). chain::Anchor is an identical copy used only by the shared DP; convert at
// the dp_anchors_v2 boundary.
use crate::types::*;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use rayon::prelude::*;
use rust_lapper::Interval;
use serde::{Deserialize, Serialize};
use statrs::distribution::{Binomial, DiscreteCDF};
use std::io::Write;
use std::path::PathBuf;
use std::sync::Mutex;

#[inline]
fn to_chain_anchor(a: &Anchor) -> crate::chain::Anchor {
    crate::chain::Anchor { i: a.i, j: a.j, pos1: a.pos1, pos2: a.pos2 }
}
#[inline]
fn from_chain_anchor(a: &crate::chain::Anchor) -> Anchor {
    Anchor { i: a.i, j: a.j, pos1: a.pos1, pos2: a.pos2 }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct HitInfo {
    //pub index: u32,
    pub position: u32,
    /// Bit 31 = canonical-strand flag of the reference minimizer (same encoding as
    /// `anchor.pos1`): 1 = canonical form is the reverse complement.
    /// Bits 30-0 = actual contig id.
    pub contig_id_strand: u32,
}


pub struct Anchors {
    pub anchors: Vec<Anchor>,
    pub max_mult: usize,
}

/// Pre-built index of one read's SNPmers (masked), reused across many comparisons.
/// Build once per outer read; query with each inner read's overlap-region SNPmers.
pub struct SnpmerIndex {
    pub index: FxHashMap<Kmer48, Vec<HitInfo>>,
}

impl SnpmerIndex {
    pub fn build(snpmers: &[(u32, FlagKmer48)], k: usize) -> Self {
        let mask = !(3u64 << (k - 1));
        let mut index: FxHashMap<Kmer48, Vec<HitInfo>> = FxHashMap::default();
        for &(pos, snpmer) in snpmers {
            let masked = Kmer48::from_u64(snpmer.kmer().to_u64() & mask);
            // contig_id = 0 (single read); strand in bit 31
            let contig_id_strand = (snpmer.strand() as u32) << 31;
            index.entry(masked).or_default().push(HitInfo {
                position: pos,
                contig_id_strand,
            });
        }
        SnpmerIndex { index }
    }
}

pub fn find_exact_matches_with_full_index(
    seq1: &[(u32, FlagKmer48)],
    index: &FxHashMap<Kmer48, Vec<HitInfo>>,
    _reference_seqs_owned: Option<&FxHashMap<usize, TwinRead>>,
    _reference_seqs_ref: Option<&FxHashMap<usize, &TwinRead>>,
) -> FxHashMap<u32, Anchors> {
    let mut max_mult = 0;
    let mut matches = FxHashMap::default();

    for (pos, flag_kmer) in seq1.iter() {
        let s1 = flag_kmer.strand() as u32;
        if let Some(indices) = index.get(&flag_kmer.kmer()) {
            if indices.len() > max_mult {
                max_mult = indices.len();
            }
            for hit in indices {
                let s2 = hit.contig_id_strand >> 31;
                let contig = hit.contig_id_strand & 0x7FFF_FFFF;
                let rel_strand = s1 ^ s2;
                let anchor = AnchorBuilder {
                    pos1: (rel_strand << 31) | *pos,
                    pos2: hit.position,
                };
                matches.entry(contig).or_insert(vec![]).push(anchor);
            }
        }
    }

    matches
        .into_iter()
        .map(|(k, v)| {
            let mut anchors = v
                .into_iter()
                .map(|anchor| Anchor {
                    i: None,
                    j: None,
                    pos1: anchor.pos1,
                    pos2: anchor.pos2,
                })
                .collect::<Vec<_>>();
            // Sort by (encoded pos1 as u32, pos2) so forward anchors precede reverse.
            anchors.sort_by_key(|a| (a.pos1, a.pos2));
            (k, Anchors { anchors, max_mult })
        })
        .collect()
}

fn find_exact_matches_indexes(
    seq1: &[(u32, FlagKmer48)],
    seq2: &[(u32, FlagKmer48)],
) -> (Vec<Anchor>, usize) {
    let mut max_mult = 0;
    let mut matches = Vec::new();

    // Index seq2: key = plain Kmer48 (no strand flag), value = (index, pos, is_rc).
    let mut index_map: FxHashMap<Kmer48, Vec<(usize, u32, bool)>> = FxHashMap::default();
    for (j, &(pos, kmer)) in seq2.iter().enumerate() {
        index_map
            .entry(kmer.kmer())
            .or_insert_with(Vec::new)
            .push((j, pos, kmer.strand()));
    }

    for (i, &(pos1, kmer)) in seq1.iter().enumerate() {
        let s1 = kmer.strand() as u32;
        if let Some(indices) = index_map.get(&kmer.kmer()) {
            if indices.len() > max_mult {
                max_mult = indices.len();
            }
            for &(j, pos2, s2) in indices {
                // rel_strand = 0 → same strand → forward chain (pos2 increases)
                // rel_strand = 1 → opposite strands → reverse chain (pos2 decreases)
                let rel_strand = s1 ^ (s2 as u32);
                let anchor = Anchor {
                    i: Some(i as u32),
                    j: Some(j as u32),
                    pos1: (rel_strand << 31) | pos1,
                    pos2,
                };
                matches.push(anchor);
            }
        }
    }

    (matches, max_mult)
}

fn find_exact_matches_with_cached_read1_index(
    splitmers1: &[(u32, FlagKmer48)],
    splitmers2: &[(u32, FlagKmer48)],
    cached: &SnpmerIndex,
    start1: usize,
    end1: usize,
) -> (Vec<Anchor>, usize) {
    // Small pos→splitmers1_idx map so we can recover `i` from a read1 position hit.
    // splitmers1 is already position-sorted, so this is O(|splitmers1|).
    let pos_to_i: FxHashMap<u32, u32> = splitmers1
        .iter()
        .enumerate()
        .map(|(i, &(pos, _))| (pos, i as u32))
        .collect();

    let mut anchors: Vec<Anchor> = Vec::new();
    let mut max_mult = 0usize;

    for (j, &(pos2, kmer2)) in splitmers2.iter().enumerate() {
        let masked_kmer = kmer2.kmer();
        if let Some(hits) = cached.index.get(&masked_kmer) {
            if hits.len() > max_mult {
                max_mult = hits.len();
            }
            for hit in hits {
                let read1_pos = hit.position as usize;
                if read1_pos < start1 || read1_pos > end1 {
                    continue;
                }
                let i = match pos_to_i.get(&hit.position) {
                    Some(&i) => i,
                    None => continue,
                };
                let s1 = hit.contig_id_strand >> 31; // read1 strand (bit 31 of contig_id_strand)
                let s2 = kmer2.strand() as u32;
                let rel_strand = s1 ^ s2;
                anchors.push(Anchor {
                    i: Some(i),
                    j: Some(j as u32),
                    pos1: (rel_strand << 31) | hit.position,
                    pos2,
                });
            }
        }
    }

    (anchors, max_mult)
}

fn find_optimal_chain(
    anchors: &Vec<Anchor>,
    match_score: i32,
    gap_cost: i32,
    band_opt: Option<usize>,
    tr_options: &CompareTwinReadOptions,
) -> Vec<ChainInfo> {
    let band;
    let matches = anchors;
    let max_gap = tr_options.max_gap;
    let double_gap = tr_options.double_gap;
    if band_opt.is_none() {
        band = 50;
    } else {
        band = band_opt.unwrap();
    }

    if anchors.is_empty() {
        return vec![];
    }

    // Single strand-aware DP pass — both forward and reverse chains in one call.
    // dp_anchors_v2 works on chain::Anchor; convert to/from types::Anchor.
    let matches_chain: Vec<crate::chain::Anchor> =
        matches.iter().map(to_chain_anchor).collect();
    let mut all_chains = dp_anchors_v2(
        &matches_chain,
        gap_cost,
        match_score,
        band,
        tr_options.max_skip,
        max_gap,
        double_gap,
        tr_options.min_chain_length,
    );

    if all_chains.is_empty() {
        return vec![];
    }

    // if tr_options.debug{
    //     log::debug!("Found {} chains before filtering", all_chains.len());
    // }

    let max_score = all_chains.iter().map(|x| x.0).max().unwrap();
    let mut chains = vec![];
    let mut reference_intervals: Vec<Interval<u32, bool>> = vec![];
    let mut query_intervals: Vec<Interval<u32, bool>> = vec![];
    all_chains.sort_unstable_by(|a, b| b.0.cmp(&a.0));

    for (score, chain, reverse) in all_chains {
        let large_indel = false;
        let cond1 = score as f64
            > tr_options.supplementary_threshold_ratio.unwrap_or(0.25) * max_score as f64;
        let cond2 = score as f64 > tr_options.supplementary_threshold_score.unwrap_or(f64::MAX);
        if cond1 || cond2 {
            let l = chain.first().unwrap().pos2;
            let r = chain.last().unwrap().pos2;
            let interval = Interval {
                start: l.min(r),
                stop: l.max(r),
                val: true,
            };

            if reference_intervals.iter().any(|x| {
                let intersect = x.intersect(&interval);
                intersect as f64 / (interval.stop - interval.start) as f64 > 0.25
            }) {
                if tr_options.force_ref_nonoverlap {
                    continue;
                }
            }

            reference_intervals.push(interval);

            let l_q = chain.first().unwrap().pos1;
            let r_q = chain.last().unwrap().pos1;
            let interval_q = Interval {
                start: l_q.min(r_q),
                stop: l_q.max(r_q),
                val: true,
            };

            if query_intervals.iter().any(|x| {
                let intersect = x.intersect(&interval_q);
                intersect as f64 / (interval_q.stop - interval_q.start) as f64 > 0.25
            }) {
                if tr_options.force_query_nonoverlap {
                    continue;
                } else {
                    let secondary_ratio = tr_options.secondary_threshold.unwrap_or(0.50);
                    if (score as f64) < secondary_ratio * (max_score as f64) {
                        continue;
                    }
                }
            }

            query_intervals.push(interval_q);

            let chain: Vec<Anchor> = chain.iter().map(from_chain_anchor).collect();
            chains.push(ChainInfo {
                chain,
                reverse: reverse,
                score: score,
                large_indel: large_indel,
            });
        }
    }

    // if tr_options.debug{
    // log::debug!("Filtering chains took {:?}, {} chains remain after filtering", timer.elapsed(), chains.len());
    // }

    return chains;
}

pub fn compare_twin_reads(
    seq1: &TwinRead,
    seq2: &TwinRead,
    mini_anchors: Option<&Anchors>,
    snpmer_anchors: Option<&Anchors>,
    cached_read1_snpmer_index: Option<&SnpmerIndex>,
    i: usize,
    j: usize,
    options: &CompareTwinReadOptions,
    args: &Cli,
) -> Vec<TwinOverlap> {
    let mut mini_chain_infos;
    let time = std::time::Instant::now();
    if let Some(anchors) = mini_anchors {
        mini_chain_infos = find_optimal_chain(
            &anchors.anchors,
            args.c as i32,
            1,
            Some(anchors.max_mult * 20),
            options,
        );
    } else {
        let anchors;
        if let Some(seq1_minimizers) = options.read1_mininimizers.as_ref() {
            anchors = find_exact_matches_indexes(seq1_minimizers, &seq2.minimizers_vec_strand());
        } else {
            anchors = find_exact_matches_indexes(
                &seq1.minimizers_vec_strand(),
                &seq2.minimizers_vec_strand(),
            );
        }
        mini_chain_infos = find_optimal_chain(
            &anchors.0,
            args.c as i32,
            1,
            Some((anchors.1 * 20).min(MAX_MULTIPLICITY_KMER)),
            options,
        );
    }

    if options.maximal_only {
        mini_chain_infos.retain(|mini_chain_info| {
            let mini_chain = &mini_chain_info.chain;
            // Use positions from the Anchor struct directly
            let l1 = mini_chain[0].pos1;
            let r1 = mini_chain[mini_chain.len() - 1].pos1;
            let l2 = mini_chain[0].pos2;
            let r2 = mini_chain[mini_chain.len() - 1].pos2;
            let start1 = l1.min(r1);
            let end1 = l1.max(r1) + seq1.k as u32 - 1;
            let start2 = l2.min(r2);
            let end2 = l2.max(r2) + seq2.k as u32 - 1;
            let shared_minimizers = mini_chain.len();
            let end_fuzz_pair = seq1.overlap_hang_length.unwrap();
            let end_fuzz = end_fuzz_pair.0.max(end_fuzz_pair.1);
            let max_mapping = check_maximal_overlap(
                start1 as usize,
                end1 as usize,
                start2 as usize,
                end2 as usize,
                seq1.base_length,
                seq2.base_length,
                mini_chain_info.reverse,
                end_fuzz,
            );
            if !max_mapping {
                return false;
            }
            if (shared_minimizers as f64)
                < (seq1.base_length as f64 / args.c as f64 / args.absolute_minimizer_cut_ratio)
            {
                return false;
            }
            return true;
        });
    }

    mini_chain_infos.retain(|mini_chain_info| mini_chain_info.score >= MIN_CHAIN_SCORE_COMPARE);

    let mini_chain_infos = mini_chain_infos;

    if mini_chain_infos.is_empty() {
        return vec![];
    }

    let mut twin_overlaps = vec![];
    let k = seq1.k as usize;

    let temp_vec;
    let mut snpmer_vec = &vec![];

    let temp_vec2;
    let mut snpmer_vec_2 = &vec![];

    // We have to populate snpmers_vec() which may be large for contigs
    if options.compare_snpmers {
        if let Some(snpmers_vec_1) = options.read1_snpmers.as_ref() {
            snpmer_vec = snpmers_vec_1;
        } else {
            temp_vec = seq1.snpmers_vec_strand();
            snpmer_vec = &temp_vec;
        }

        temp_vec2 = seq2.snpmers_vec_strand();
        snpmer_vec_2 = &temp_vec2;
    }

    let time_snpchain = std::time::Instant::now();
    for mini_chain_info in mini_chain_infos {
        let mini_chain = &mini_chain_info.chain;

        let mut shared_snpmer = usize::MAX;
        let mut diff_snpmer = usize::MAX;

        if options.compare_snpmers {
            shared_snpmer = 0;
            diff_snpmer = 0;

            // Use positions from the Anchor struct directly
            let l1 = mini_chain[0].pos1 as usize;
            let r1 = mini_chain[mini_chain.len() - 1].pos1 as usize;
            let l2 = mini_chain[0].pos2 as usize;
            let r2 = mini_chain[mini_chain.len() - 1].pos2 as usize;
            let start1 = l1.min(r1);
            let end1 = l1.max(r1) + k - 1;
            let start2 = l2.min(r2);
            let end2 = l2.max(r2) + k - 1;

            let mask = !(3 << (k - 1));

            let mut splitmers1: Vec<(u32, FlagKmer48)> = vec![];
            let mut ind_redirect1 = vec![];

            for (i, &(pos, snpmer)) in snpmer_vec.iter().enumerate() {
                if pos as usize >= start1 && pos as usize <= end1 {
                    ind_redirect1.push(i);
                    let masked = FlagKmer48::new(
                        Kmer48::from_u64(snpmer.kmer().to_u64() & mask),
                        snpmer.strand(),
                    );
                    splitmers1.push((pos, masked));
                }
            }

            let mut splitmers2: Vec<(u32, FlagKmer48)> = vec![];
            let mut ind_redirect2 = vec![];

            for (i, &(pos, snpmer)) in snpmer_vec_2.iter().enumerate() {
                if pos as usize >= start2 && pos as usize <= end2 {
                    ind_redirect2.push(i);
                    let masked = FlagKmer48::new(
                        Kmer48::from_u64(snpmer.kmer().to_u64() & mask),
                        snpmer.strand(),
                    );
                    splitmers2.push((pos, masked));
                }
            }

            let split_chain_opt;
            let mut split_options = options.clone();
            split_options.min_chain_length = 2;
            split_options.double_gap = 2_000_000;
            if let Some(anchors) = snpmer_anchors {
                split_chain_opt = find_optimal_chain(
                    &anchors.anchors,
                    50,
                    1,
                    Some((anchors.max_mult * 10).min(50)),
                    &split_options,
                )
                .into_iter()
                .max_by_key(|x| x.score);
            } else if let Some(cached) = cached_read1_snpmer_index {
                let anchors = find_exact_matches_with_cached_read1_index(
                    &splitmers1,
                    &splitmers2,
                    cached,
                    start1,
                    end1,
                );
                let chains = find_optimal_chain(
                    &anchors.0,
                    50,
                    1,
                    Some((anchors.1 * 10).min(50)),
                    &split_options,
                );
                split_chain_opt = chains.into_iter().max_by_key(|x| x.score);
            } else {
                let anchors = find_exact_matches_indexes(&splitmers1, &splitmers2);
                let chains = find_optimal_chain(
                    &anchors.0,
                    50,
                    1,
                    Some((anchors.1 * 10).min(50)),
                    &split_options,
                );
                split_chain_opt = chains.into_iter().max_by_key(|x| x.score);
            }

            //If mini chain goes opposite from split chain, probably split chain
            //is not reliable, so set shared and diff = 0.
            if let Some(split_chain) = split_chain_opt.as_ref() {
                if split_chain.reverse == mini_chain_info.reverse || split_chain.chain.len() == 1 {
                    let snpmer_kmers_seq1 = &snpmer_vec;
                    let snpmer_kmers_seq2 = &snpmer_vec_2;
                    for anchor in split_chain.chain.iter() {
                        let i = anchor.i;
                        let i = ind_redirect1[i.unwrap() as usize];
                        let j = anchor.j;
                        let j = ind_redirect2[j.unwrap() as usize];

                        //if seq1.snpmer_kmers[i as usize] == seq2.snpmer_kmers[j as usize] {
                        if snpmer_kmers_seq1[i as usize].1.kmer()
                            == snpmer_kmers_seq2[j as usize].1.kmer()
                        {
                            shared_snpmer += 1;
                        } else {
                            diff_snpmer += 1;
                        }
                    }
                }
            }

            //Only if log level is trace
            if log::log_enabled!(log::Level::Trace) && true {
                if diff_snpmer < 10 && shared_snpmer > 100 {
                    let mut positions_read1_snpmer_diff = vec![];
                    let mut positions_read2_snpmer_diff = vec![];

                    let mut kmers_read1_diff = vec![];
                    let mut kmers_read2_diff = vec![];

                    let snpmer_kmers_seq1 = seq1.snpmer_kmers();
                    let snpmer_kmers_seq2 = seq2.snpmer_kmers();
                    let snpmer_positions_seq1 = seq1.snpmer_positions();
                    let snpmer_positions_seq2 = seq2.snpmer_positions();

                    for anchor in split_chain_opt.unwrap().chain.iter() {
                        let i = anchor.i;
                        let i = ind_redirect1[i.unwrap() as usize];
                        let j = anchor.j;
                        let j = ind_redirect2[j.unwrap() as usize];
                        if snpmer_kmers_seq1[i as usize] != snpmer_kmers_seq2[j as usize] {
                            positions_read1_snpmer_diff.push(snpmer_positions_seq1[i as usize]);
                            positions_read2_snpmer_diff.push(snpmer_positions_seq2[j as usize]);

                            let kmer1 = decode_kmer48(snpmer_kmers_seq1[i as usize], seq1.k as u8);
                            let kmer2 = decode_kmer48(snpmer_kmers_seq2[j as usize], seq2.k as u8);

                            kmers_read1_diff.push(kmer1);
                            kmers_read2_diff.push(kmer2);
                        }
                    }
                    log::trace!(
                        "{}--{:?} {}--{:?}, snp_diff:{} snp_shared:{}, kmers1:{:?}, kmers2:{:?}",
                        &seq1.id,
                        positions_read1_snpmer_diff,
                        &seq2.id,
                        positions_read2_snpmer_diff,
                        diff_snpmer,
                        shared_snpmer,
                        kmers_read1_diff,
                        kmers_read2_diff
                    );
                }
            }
        }

        // Use positions from the Anchor struct directly
        let l1 = mini_chain[0].pos1;
        let r1 = mini_chain[mini_chain.len() - 1].pos1;
        let l2 = mini_chain[0].pos2;
        let r2 = mini_chain[mini_chain.len() - 1].pos2;
        let start1 = l1.min(r1);
        let end1 = l1.max(r1) + k as u32 - 1;
        let start2 = l2.min(r2);
        let end2 = l2.max(r2) + k as u32 - 1;
        let shared_minimizers = mini_chain.len();
        let mut mini_chain_return = None;
        if options.retain_chain {
            mini_chain_return = Some(mini_chain_info.chain);
        }
        let twinol = TwinOverlap {
            i1: i,
            i2: j,
            start1: start1 as usize,
            end1: end1 as usize,
            start2: start2 as usize,
            end2: end2 as usize,
            shared_minimizers,
            shared_snpmers: shared_snpmer,
            diff_snpmers: diff_snpmer,
            snpmers_in_both: (seq1.snpmer_count(), seq2.snpmer_count()),
            chain_reverse: mini_chain_info.reverse,
            chain_score: mini_chain_info.score,
            minimizer_chain: mini_chain_return,
            large_indel: mini_chain_info.large_indel,
        };
        twin_overlaps.push(twinol);
    }
    if options.debug {
        log::debug!(
            "Compared SNPmers for read {} to read {} in {:?}",
            seq1.id,
            seq2.id,
            time_snpchain.elapsed()
        );
        log::debug!(
            "Finished comparing read {} to read {} in {:?}, found {} overlaps",
            seq1.id,
            seq2.id,
            time.elapsed(),
            twin_overlaps.len()
        );
    }
    return twin_overlaps;
}

pub fn id_est(shared_minimizers: usize, diff_snpmers: usize, c: u64, large_indel: bool) -> f64 {
    let diff_snps = diff_snpmers as f64;
    let shared_minis = shared_minimizers as f64;
    let alpha = diff_snps as f64 / shared_minis as f64 / c as f64;
    let theta = alpha / (1. + alpha);
    let mut id_est = 1. - theta;

    if large_indel {
        //Right now it's 0.5% penalty for large indels, but we don't use largeindels.
        let penalty = IDENTITY_THRESHOLDS.last().unwrap() - IDENTITY_THRESHOLDS.first().unwrap();
        let penalty = penalty / 2.;
        id_est -= penalty;
    }

    return id_est;
}

pub fn get_minimizer_index(
    tr_owned: Option<&FxHashMap<usize, TwinRead>>,
    tr_ref: Option<&FxHashMap<usize, &TwinRead>>,
) -> FxHashMap<Kmer48, Vec<HitInfo>> {
    let mut mini_index = FxHashMap::default();
    if let Some(twinreads) = tr_owned {
        let mut sorted_keys = twinreads.keys().collect::<Vec<_>>();
        sorted_keys.sort();
        for (&id, tr) in sorted_keys.iter().map(|&x| (x, &twinreads[x])) {
            for (pos, mini) in tr.minimizers_vec_strand().into_iter() {
                let hit = HitInfo {
                    contig_id_strand: (mini.strand() as u32) << 31 | id as u32,
                    position: pos,
                };
                mini_index.entry(mini.kmer()).or_insert(vec![]).push(hit);
            }
        }
    } else if let Some(twinreads) = tr_ref {
        let mut sorted_keys = twinreads.keys().collect::<Vec<_>>();
        sorted_keys.sort();
        for (&id, tr) in sorted_keys.iter().map(|&x| (x, &twinreads[x])) {
            for (pos, mini) in tr.minimizers_vec_strand().into_iter() {
                let hit = HitInfo {
                    contig_id_strand: (mini.strand() as u32) << 31 | id as u32,
                    position: pos,
                };
                mini_index.entry(mini.kmer()).or_insert(vec![]).push(hit);
            }
        }
    } else {
        panic!("No minimizer index provided");
    }

    if mini_index.len() == 0 {
        return mini_index;
    }

    let mut minimizer_to_hit_count = mini_index.iter().map(|(_, v)| v.len()).collect::<Vec<_>>();

    minimizer_to_hit_count.sort_by(|a, b| b.cmp(&a));
    let threshold = minimizer_to_hit_count[minimizer_to_hit_count.len() / 100_000];
    log::trace!(
        "Minimizer index size: {}. Threshold: {}",
        minimizer_to_hit_count.len(),
        threshold
    );

    // Only threshold when necessary
    if mini_index.len() > 500_000 {
        mini_index.retain(|_, v| v.len() < threshold);
    }

    mini_index.shrink_to_fit();

    return mini_index;
}

pub fn check_maximal_overlap(
    start1: usize,
    end1: usize,
    start2: usize,
    end2: usize,
    len1: usize,
    len2: usize,
    reverse: bool,
    endpoint_fuzz: usize,
) -> bool {
    let edge_fuzz = endpoint_fuzz;

    //Can not extend to the left (cond1) and cannot extend to the right (cond2)
    //  ------->             OR          --------->
    // ---------->                            --------->
    if !reverse {
        if (start1 < edge_fuzz || start2 < edge_fuzz)
            && (len1 < edge_fuzz + end1 || len2 < edge_fuzz + end2)
        {
            return true;
        }
    }
    //  ------->             OR          --------->            OR          <----------
    //      <------                    <--------------                            -------->
    else {
        let max1right = len1 < edge_fuzz + end1;
        let max2right = len2 < edge_fuzz + end2;
        let max1left = start1 < edge_fuzz;
        let max2left = start2 < edge_fuzz;

        let ol_plus_minus = max1right && max2right;
        let ol_minus_plus = max1left && max2left;
        let contained1 = max1left && max1right;
        let contained2 = max2left && max2right;

        if ol_plus_minus || ol_minus_plus || contained1 || contained2 {
            return true;
        }
    }

    return false;
}

#[derive(Debug, Clone, PartialEq, Default, Hash, Eq, Serialize, Deserialize, PartialOrd, Ord)]
pub struct OverlapConfig {
    pub read_i: usize,
    pub read_j: usize,
    pub diff_snpmer: usize,
    pub hang1: usize,
    pub hang2: usize,
    pub forward1: bool,
    pub forward2: bool,
    pub overlap1_len: usize,
    pub overlap2_len: usize,
    pub shared_mini: usize,
    pub shared_snpmer: usize,
    pub contained: bool,
    pub large_indel: bool,
    // Raw aligned intervals on each read (forward coordinates), copied from the
    // TwinOverlap. Added for the FFI so the host gets explicit q/t boxes; `reverse`
    // is the relative strand (read_j reverse-complemented w.r.t. read_i).
    pub start1: usize,
    pub end1: usize,
    pub start2: usize,
    pub end2: usize,
    pub reverse: bool,
}

pub fn get_overlaps_outer_reads_twin(
    twin_reads: &[TwinRead],
    outer_read_indices: &[usize],
    args: &Cli,
    overlap_file_path: Option<&PathBuf>,
    contained_file_path: Option<&PathBuf>,
) -> Vec<OverlapConfig> {
    let _ = (overlap_file_path, contained_file_path);
    let bufwriter = Mutex::new(Box::new(std::io::sink()) as Box<dyn Write + Send>);
    let bufwriter_contained = Mutex::new(Box::new(std::io::sink()) as Box<dyn Write + Send>);

    // Contained reads may still lurk in the outer reads, remove again for sure.
    let contained_reads_again = Mutex::new(FxHashSet::default());
    let overlaps = Mutex::new(vec![]);

    let outer_index_batches = outer_read_indices.chunks(args.read_map_batch_size);

    for outer_index_batch in outer_index_batches {
        let outer_twin_reads_batch = outer_index_batch
            .iter()
            .map(|x| (*x, &twin_reads[*x]))
            .collect::<FxHashMap<_, _>>();
        let inverted_index = get_minimizer_index(None, Some(&outer_twin_reads_batch));
        outer_read_indices.into_par_iter().for_each(|&i| {
            let read = &twin_reads[i];
            let mini_anchors = find_exact_matches_with_full_index(&read.minimizers_vec_strand(), &inverted_index, None, Some(&outer_twin_reads_batch));

            let comparison_options = CompareTwinReadOptions{
                compare_snpmers: true,
                retain_chain: false,
                //force_one_to_one_alignments: true,
                force_query_nonoverlap: true,
                force_ref_nonoverlap: true,
                supplementary_threshold_score: Some(500.),
                supplementary_threshold_ratio: Some(0.25),
                secondary_threshold: None,
                read1_mininimizers: None, // indexed anchors
                read1_snpmers: Some(read.snpmers_vec_strand()),
                max_gap: MAX_GAP_CHAINING * 3/2,
                double_gap: 25_00,
                max_skip: 10,
                min_chain_length: 3,
                maximal_only: false, // Already use "dovetail possibility" as a similar filter
                debug: false
            };

            let read1_snpmer_idx = comparison_options.read1_snpmers.as_ref()
                .map(|v| SnpmerIndex::build(v, args.kmer_size));

            for (outer_ref_id, anchors) in mini_anchors.into_iter(){

                //Only compare once. I think we get slightly different results if we
                //compare in both directinos, but this forces consistency.
                if i <= outer_ref_id as usize {
                    continue;
                }

                let read2 = &twin_reads[outer_ref_id as usize];

                if !dovetail_possibility(&anchors, &read, &read2){
                    continue;
                }

                let twlaps = compare_twin_reads(&read, &read2, Some(&anchors), None, read1_snpmer_idx.as_ref(), i, outer_ref_id as usize, &comparison_options, args);

                if twlaps.len() > 1{
                    log::trace!("Multiple overlaps for {}:{} and {}:{}", &read.id, i,  &read2.id, outer_ref_id);
                    for twlap in twlaps.iter(){
                        log::trace!("{}-{} Overlap: {}-{} {}-{}, reverse {}", i, outer_ref_id, twlap.start1, twlap.end1, twlap.start2, twlap.end2, twlap.chain_reverse);
                    }
                }

                let mut possible_containment = false;
                //Check for contained read
                for twlap in twlaps.iter(){
                    let mut twlap_contain = false;
                    let mut smaller_read_index = i;
                    let mut snpmer_threshold = args.snpmer_threshold_strict;
                    if r1_contained_r2(&twlap, &read, &read2, true, args.c, args.hifi){
                        twlap_contain = true;
                        possible_containment = true;
                        snpmer_threshold = read.snpmer_id_threshold.unwrap_or(100.);
                    }
                    else if r1_contained_r2(&twlap, &read2, &read, true, args.c, args.hifi){
                        twlap_contain = true;
                        possible_containment = true;
                        smaller_read_index = outer_ref_id as usize;
                        snpmer_threshold = read2.snpmer_id_threshold.unwrap_or(100.);
                    }
                    if twlap_contain {
                        if !args.no_containment_removal && same_strain(twlap.shared_minimizers, twlap.diff_snpmers, twlap.shared_snpmers, args.c as u64, snpmer_threshold, args.snpmer_error_rate_strict, twlap.large_indel){
                            contained_reads_again.lock().unwrap().insert(smaller_read_index);

                            writeln!(bufwriter_contained.lock().unwrap(), "{} ({}) {} ({}), SMALLER: {}, LEN1: {} RANGE1: {}-{}, LEN2:{} RANGE2: {}-{}, SNP_DIFF: {}, SNP_SHARE: {}, MINI: {}",
                                read.id, i, read2.id, outer_ref_id, smaller_read_index,
                                twin_reads[i].base_length, twlap.start1, twlap.end1,
                                twin_reads[outer_ref_id as usize].base_length, twlap.start2, twlap.end2,
                                twlap.diff_snpmers, twlap.shared_snpmers, twlap.shared_minimizers).unwrap();
                        }
                    }
                }

                if possible_containment && !args.no_containment_removal {
                    if twlaps.len() == 1 {
                        let twlap = &twlaps[0];
                        let contained_overlap_config = OverlapConfig{
                            hang1: 0,
                            hang2: 0,
                            //Forward doesn't make sense for contained reads
                            forward1: true,
                            forward2: true,
                            overlap1_len: twlap.end1 - twlap.start1,
                            overlap2_len: twlap.end2 - twlap.start2,
                            read_i: twlap.i1,
                            read_j: twlap.i2,
                            start1: twlap.start1,
                            end1: twlap.end1,
                            start2: twlap.start2,
                            end2: twlap.end2,
                            reverse: twlap.chain_reverse,
                            shared_mini: twlap.shared_minimizers,
                            shared_snpmer: twlap.shared_snpmers,
                            diff_snpmer: twlap.diff_snpmers,
                            contained: true,
                            large_indel: twlap.large_indel
                        };
                        log::trace!("Contained read {} in read {}", twlap.i2, twlap.i1);
                        overlaps.lock().unwrap().push(contained_overlap_config);
                    }
                    continue
                }

                let best_overlaps = comparison_to_overlap(twlaps, &twin_reads, args, &bufwriter);
                for best_ol in best_overlaps{
                    overlaps.lock().unwrap().push(best_ol);
                }
            }
        });
    }

    let contained_reads = contained_reads_again.into_inner().unwrap();
    let mut ol = overlaps.into_inner().unwrap();
    ol.retain(|x| !contained_reads.contains(&x.read_i) && !contained_reads.contains(&x.read_j));
    ol.sort_by_key(|x| {
        (
            -((x.overlap1_len + x.overlap2_len) as i32) + (x.hang1 + x.hang2) as i32,
            x.read_i,
            x.read_j,
        )
    });

    return ol;
}

pub fn comparison_to_overlap<T>(
    twlaps: Vec<TwinOverlap>,
    twin_reads: &[TwinRead],
    args: &Cli,
    writer: &Mutex<T>,
) -> Vec<OverlapConfig>
where
    T: Write + Send,
{
    let mut overlap_possib_out = vec![];
    let mut overlap_possib_in = vec![];

    for twlap in twlaps {
        let overlap_possibility = overlap_config_from_twlap(&twlap, twin_reads, args, writer);
        if let Some(ol) = overlap_possibility {
            let first_index = ol.read_i == twlap.i1;
            if first_index && ol.forward1 {
                overlap_possib_out.push(ol);
            } else if first_index && !ol.forward1 {
                overlap_possib_in.push(ol);
            } else if !first_index && ol.forward2 {
                overlap_possib_in.push(ol);
            } else if !first_index && !ol.forward2 {
                overlap_possib_out.push(ol);
            }
        }
    }

    let best_overlap_forward = overlap_possib_out
        .into_iter()
        .max_by_key(|x| (x.overlap1_len + x.overlap2_len) as i32 - (x.hang1 + x.hang2) as i32);
    let best_overlap_backward = overlap_possib_in
        .into_iter()
        .max_by_key(|x| (x.overlap1_len + x.overlap2_len) as i32 - (x.hang1 + x.hang2) as i32);
    let mut best_overlaps = vec![];

    if best_overlap_forward.is_some() && best_overlap_backward.is_none() {
        best_overlaps.push(best_overlap_forward.unwrap());
    } else if best_overlap_forward.is_none() && best_overlap_backward.is_some() {
        best_overlaps.push(best_overlap_backward.unwrap());
    } else if best_overlap_forward.is_some() && best_overlap_backward.is_some() {
        // Check if circular overlap is concordant.
        let ol_f = best_overlap_forward.unwrap();
        let ol_b = best_overlap_backward.unwrap();

        let reverse_f = ol_f.forward1 != ol_f.forward2;
        let reverse_b = ol_b.forward1 != ol_b.forward2;

        //Consistent
        if reverse_f == reverse_b {
            best_overlaps.push(ol_f);
            best_overlaps.push(ol_b);
        }
        //Inconsistent; take best
        else {
            let best_ol = [ol_f, ol_b]
                .into_iter()
                .max_by_key(|x| {
                    (x.overlap1_len + x.overlap2_len) as i32 - (x.hang1 + x.hang2) as i32
                })
                .unwrap();
            best_overlaps.push(best_ol);
        }
    }

    return best_overlaps;
}

fn overlap_config_from_twlap<T>(
    twlap: &TwinOverlap,
    twin_reads: &[TwinRead],
    args: &Cli,
    writer: &Mutex<T>,
) -> Option<OverlapConfig>
where
    T: Write + Send,
{
    let mut overlap_possibilities = vec![];

    let i = twlap.i1;
    let j = twlap.i2;
    let mut exist_overlap = false;

    let read1 = &twin_reads[i];
    let read2 = &twin_reads[j];

    //check if end-to-end overlap
    let identity = id_est(
        twlap.shared_minimizers,
        twlap.diff_snpmers,
        args.c as u64,
        twlap.large_indel,
    );
    let same_strain_lax = args.no_same_strain_filter
        || same_strain(
            twlap.shared_minimizers,
            twlap.diff_snpmers,
            twlap.shared_snpmers,
            args.c as u64,
            args.snpmer_threshold_lax,
            args.snpmer_error_rate_lax,
            twlap.large_indel,
        );

    let aln_len1 = twlap.end1 - twlap.start1;
    let aln_len2 = twlap.end2 - twlap.start2;

    if aln_len1.max(aln_len2) < args.min_ol {
        return None;
    }

    let (hang1_start, hang1_end) = read1.overlap_hang_length.unwrap();
    let (hang2_start, hang2_end) = read2.overlap_hang_length.unwrap();

    let hang_start = hang1_start.max(hang2_start);
    let hang_end = hang1_end.max(hang2_end);

    //let (ext_s1, ext_e1, ext_s2, ext_e2) = alignment::extend_ends_chain(&twin_reads[twlap.i1].dna_seq, &twin_reads[twlap.i2].dna_seq, twlap, args);

    //let aln_len1 = ext_e1 - ext_s1 + 1;
    //let aln_len2 = ext_e2 - ext_s2 + 1;

    //let mini_chain = twlap.minimizer_chain.as_ref().unwrap();

    // println!("EXT {}, {}, {}, {} READ {}", ext_s1, ext_e1, ext_s2, ext_e2, &read1.id);
    // println!("OLRANGE {}, {}, {}, {} READ {}", twlap.start1, twlap.end1, twlap.start2, twlap.end2, &read1.id);
    // println!("HANG {}, {}, {}, {} READ {}", hang1_start, hang1_end, hang2_start, hang2_end, &read1.id);

    // let (hang1_start, hang1_end) = (OVERLAP_HANG_LENGTH, OVERLAP_HANG_LENGTH);
    // let (hang2_start, hang2_end) = (OVERLAP_HANG_LENGTH, OVERLAP_HANG_LENGTH);

    if twlap.chain_reverse {
        if twlap.start1 < hang_start && twlap.start2 < hang_start {
            let ol_config = OverlapConfig {
                hang1: twlap.start1,
                hang2: twlap.start2,
                forward1: false,
                forward2: true,
                overlap1_len: aln_len1,
                overlap2_len: aln_len2,
                read_i: i,
                read_j: j,
                start1: twlap.start1,
                end1: twlap.end1,
                start2: twlap.start2,
                end2: twlap.end2,
                reverse: twlap.chain_reverse,
                shared_mini: twlap.shared_minimizers,
                shared_snpmer: twlap.shared_snpmers,
                diff_snpmer: twlap.diff_snpmers,
                contained: false,
                large_indel: twlap.large_indel,
            };
            if same_strain_lax {
                overlap_possibilities.push(ol_config);
            }
            exist_overlap = true;
        } else if twlap.end1 + hang_end > read1.base_length
            && twlap.end2 + hang_end > read2.base_length
        {
            let ol_config = OverlapConfig {
                //hang1: read1.base_length - twlap.end1 - 1,
                //hang2: read2.base_length - twlap.end2 - 1,
                hang1: read1.base_length - twlap.end1 - 1,
                hang2: read2.base_length - twlap.end2 - 1,
                forward1: true,
                forward2: false,
                overlap1_len: aln_len1,
                overlap2_len: aln_len2,
                read_i: i,
                read_j: j,
                start1: twlap.start1,
                end1: twlap.end1,
                start2: twlap.start2,
                end2: twlap.end2,
                reverse: twlap.chain_reverse,
                shared_mini: twlap.shared_minimizers,
                shared_snpmer: twlap.shared_snpmers,
                diff_snpmer: twlap.diff_snpmers,
                contained: false,
                large_indel: twlap.large_indel,
            };
            if same_strain_lax {
                overlap_possibilities.push(ol_config);
            }
            exist_overlap = true;
        }
    } else {
        if twlap.start1 < hang_start && twlap.end2 + hang_end > read2.base_length {
            let ol_config = OverlapConfig {
                //hang1: twlap.start1,
                //hang2: read2.base_length - twlap.end2 - 1,
                hang1: twlap.start1,
                hang2: read2.base_length - twlap.end2 - 1,
                forward1: false,
                forward2: false,
                overlap1_len: aln_len1,
                overlap2_len: aln_len2,
                read_i: i,
                read_j: j,
                start1: twlap.start1,
                end1: twlap.end1,
                start2: twlap.start2,
                end2: twlap.end2,
                reverse: twlap.chain_reverse,
                shared_mini: twlap.shared_minimizers,
                shared_snpmer: twlap.shared_snpmers,
                diff_snpmer: twlap.diff_snpmers,
                contained: false,
                large_indel: twlap.large_indel,
            };
            if same_strain_lax {
                overlap_possibilities.push(ol_config);
            }
            exist_overlap = true;
        } else if twlap.end1 + hang_end > read1.base_length && twlap.start2 < hang_start {
            let ol_config = OverlapConfig {
                //hang1: read1.base_length - twlap.end1 - 1,
                //hang2: twlap.start2,
                hang1: read1.base_length - twlap.end1 - 1,
                hang2: twlap.start2,
                forward1: true,
                forward2: true,
                overlap1_len: aln_len1,
                overlap2_len: aln_len2,
                read_i: i,
                read_j: j,
                start1: twlap.start1,
                end1: twlap.end1,
                start2: twlap.start2,
                end2: twlap.end2,
                reverse: twlap.chain_reverse,
                shared_mini: twlap.shared_minimizers,
                shared_snpmer: twlap.shared_snpmers,
                diff_snpmer: twlap.diff_snpmers,
                contained: false,
                large_indel: twlap.large_indel,
            };
            if same_strain_lax {
                overlap_possibilities.push(ol_config);
            }
            exist_overlap = true;
        }
    }

    if exist_overlap {
        let mut bufwriter = writer.lock().unwrap();
        let mut possibilties_string = String::new();
        for possib in overlap_possibilities.iter() {
            possibilties_string.push_str(&format!(
                "{}:{}-{}:{} HANG {} {}",
                possib.read_i,
                possib.forward1,
                possib.read_j,
                possib.forward2,
                possib.hang1,
                possib.hang2
            ));
        }
        writeln!(bufwriter,
            "{} {} {} {} fsv:{} SHARE:{} DIFF:{} MINI: {}, LEN1:{} {}-{} LEN2:{} {}-{}, REVERSE: {}, Possibilties: {}",
            i,
            j,
            &read1.id.split_ascii_whitespace().next().unwrap(),
            &read2.id.split_ascii_whitespace().next().unwrap(),
            identity * 100.,
            twlap.shared_snpmers,
            twlap.diff_snpmers,
            twlap.shared_minimizers,
            read1.base_length,
            twlap.start1,
            twlap.end1,
            read2.base_length,
            twlap.start2,
            twlap.end2,
            twlap.chain_reverse,
            possibilties_string
        ).unwrap();
    }

    let best_overlap = overlap_possibilities.into_iter().max_by_key(|x| {
        (x.overlap1_len as i64 + x.overlap2_len as i64)
            - (x.hang1 as i64 + x.hang2 as i64)
            - (x.hang1 as i64 - x.hang2 as i64).abs()
    });
    return best_overlap;
}

pub fn r1_contained_r2(
    twin_overlap: &TwinOverlap,
    read1: &TwinRead,
    read2: &TwinRead,
    same_strain: bool,
    c: usize,
    hifi: bool,
) -> bool {
    let ol_len = twin_overlap.end1 - twin_overlap.start1;
    //if ol_len as f64 + (30. * c as f64) > 0.95 * (read1.base_length as f64) && same_strain && read1.base_length < read2.base_length {
    let slack = if hifi { 10. } else { 20. };
    if ol_len as f64 + (slack * c as f64) > 0.95 * read1.base_length as f64
        && same_strain
        && read1.base_length < read2.base_length
    {
        return true;
    }
    return false;
}

pub fn binomial_test(n: u64, k: u64, p: f64) -> f64 {
    // n: number of trials
    // k: number of successes
    // p: probability of success

    // Create a binomial distribution
    let binomial = Binomial::new(p, n).unwrap();

    // Calculate the probability of observing k or more successes
    let p_value = 1.0 - binomial.cdf(k);

    p_value
}

pub fn same_strain(
    minimizers: usize,
    snp_diff: usize,
    snp_shared: usize,
    c: u64,
    snpmer_threshold: f64,
    snpmer_error_rate: f64,
    large_indel: bool,
) -> bool {
    assert!(snpmer_threshold > 1.0); // Some ambiguity with percentages vs fractions...
    let identity = id_est(minimizers, snp_diff, c, large_indel);
    let high_id;
    if identity >= snpmer_threshold / 100. {
        high_id = true;
    } else {
        high_id = false;
    }
    let p_val = binomial_test(
        (snp_diff + snp_shared) as u64,
        snp_diff as u64,
        snpmer_error_rate,
    );
    let miscalled_snpmers;
    if p_val > 0.05 {
        miscalled_snpmers = true;
    } else {
        miscalled_snpmers = false;
    }

    return high_id || miscalled_snpmers;
}

fn dovetail_possibility(anchors: &Anchors, read1: &TwinRead, read2: &TwinRead) -> bool {
    let read1_length = read1.base_length;
    let read2_length = read2.base_length;

    let mut read1_possible = false;
    let mut read2_possible = false;

    //TODO change 750 to an adaptive threshold based on solid k-mers and error rates?
    for anchor in anchors.anchors.iter() {
        if anchor.pos1 < OVERLAP_HANG_LENGTH as u32
            || anchor.pos1 + OVERLAP_HANG_LENGTH as u32 > read1_length as u32
        {
            read1_possible = true;
        }
        if anchor.pos2 < OVERLAP_HANG_LENGTH as u32
            || anchor.pos2 + OVERLAP_HANG_LENGTH as u32 > read2_length as u32
        {
            read2_possible = true;
        }

        if read2_possible && read1_possible {
            return true;
        }
    }

    return false;
}
