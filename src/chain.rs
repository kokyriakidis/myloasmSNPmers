//! Anchor chaining, restored verbatim from upstream myloasm's mapping module
//! (dp_anchors_v2 + dp_inner), plus a thin `best_chain` helper for the hifiasm
//! fake-chain FFI.
//!
//! Only the chainer is restored here (not the full mapping pipeline): the
//! hifiasm bridge matches markers in C++ and calls this DP to obtain the ordered
//! colinear chain. The DP is strand-aware in a single pass: an anchor's strand
//! is carried in bit 31 of `pos1` (0 = forward, 1 = reverse); `pos2` is the raw
//! target position and the reverse encoding (`!pos2`) is applied internally.
//!
//! Kept byte-for-byte compatible with upstream so chaining behaviour matches
//! myloasm exactly. Chain defaults (see CompareTwinReadOptions / constants):
//!   match_score = 1, gap_cost = c (11), band = max_mult * 20,
//!   max_gap = 200, max_skip = 10, double_gap = 10_000, min_chain_length = 3.

/// One anchor. `i`/`j` are opaque caller indices (we use `i` to map a chained
/// anchor back to the C-side merged-anchor slot so its type tag survives).
/// `pos1` carries the query position with strand in bit 31; `pos2` is the raw
/// target position.
#[derive(Clone, Copy, Debug)]
pub struct Anchor {
    pub i: Option<u32>,
    pub j: Option<u32>,
    pub pos1: u32,
    pub pos2: u32,
}

/// Chain a set of anchors. Returns, per chain, (score, ordered anchors,
/// is_reverse). Verbatim from upstream dp_anchors_v2.
pub fn dp_anchors_v2(
    matches: &[Anchor],
    gap_cost: i32,
    match_score: i32,
    max_iter: usize,
    max_skip: usize,
    max_gap: usize,
    double_gap: usize,
    min_chain_length: usize,
) -> Vec<(i32, Vec<Anchor>, bool)> {
    if matches.is_empty() {
        return vec![];
    }
    let n = matches.len();

    let mut sorted_indices: Vec<usize> = (0..n).collect();
    sorted_indices.sort_unstable_by_key(|&i| {
        let pos1_u32 = matches[i].pos1 as i32;
        let pos2_raw = matches[i].pos2 as i32;
        let enc_pos2 = if pos1_u32 >> 31 == 1 {
            !pos2_raw
        } else {
            pos2_raw
        };
        (pos1_u32, enc_pos2)
    });

    let pos1s: Vec<i32> = sorted_indices
        .iter()
        .map(|&i| matches[i].pos1 as i32)
        .collect();
    let pos2s: Vec<i32> = sorted_indices
        .iter()
        .map(|&i| {
            let pos2_raw = matches[i].pos2;
            if matches[i].pos1 >> 31 == 1 {
                (!pos2_raw) as i32
            } else {
                pos2_raw as i32
            }
        })
        .collect();

    let mut f: Vec<i32> = vec![match_score; n];
    let mut p: Vec<i32> = vec![-1i32; n];
    let mut t: Vec<i32> = vec![-1i32; n];

    let double_gap_i32 = double_gap as i32;
    let max_gap_i32 = max_gap as i32;
    let mut st: usize = 0;
    let mut max_ii: i32 = -1;

    dp_inner(
        &pos1s,
        &pos2s,
        &mut f,
        &mut p,
        &mut t,
        gap_cost,
        match_score,
        max_iter,
        max_skip,
        double_gap_i32,
        max_gap_i32,
        &mut st,
        &mut max_ii,
        n,
    );

    // Chain reconstruction. Repurpose t[] as claimed-anchor marker.
    t.fill(0);

    let mut chains = Vec::new();

    let mut best_indices_ordered = (0..n as i32)
        .filter(|&i| f[i as usize] > min_chain_length as i32 * match_score / 2)
        .map(|i| (f[i as usize], i))
        .collect::<Vec<_>>();
    best_indices_ordered.sort_by_key(|&(score, _)| -score);

    for (score, best_index) in best_indices_ordered {
        let bi = best_index as usize;
        if t[bi] != 0 {
            continue;
        }

        let chain_is_reverse = matches[sorted_indices[bi]].pos1 >> 31 == 1;

        let mut chain = Vec::new();
        let mut idx = best_index;
        let mut good_chain = true;
        while idx >= 0 {
            let u = idx as usize;
            if t[u] != 0 {
                good_chain = false;
                break;
            }
            t[u] = 1;
            let orig = &matches[sorted_indices[u]];
            chain.push(Anchor {
                pos1: orig.pos1 & 0x7FFF_FFFF,
                ..*orig
            });
            idx = p[u];
        }

        if chain.len() < min_chain_length {
            break;
        }

        if good_chain {
            chain.reverse();
            chains.push((score, chain, chain_is_reverse));
        }
    }

    chains
}

/// Monomorphised DP kernel (upstream dp_inner). The `REV` const generic is kept
/// for byte-for-byte fidelity but is always instantiated `false`; strand is
/// handled by the pos1/pos2 encoding.
#[inline(never)]
fn dp_inner(
    pos1s: &[i32],
    pos2s: &[i32],
    f: &mut [i32],
    p: &mut [i32],
    t: &mut [i32],
    gap_cost: i32,
    match_score: i32,
    max_iter: usize,
    max_skip: usize,
    double_gap_i32: i32,
    max_gap_i32: i32,
    st: &mut usize,
    max_ii: &mut i32,
    n: usize,
) {
    const REV: bool = false;
    for i in 0..n {
        let s1 = pos1s[i];
        let s2 = pos2s[i];
        let i_i32 = i as i32;

        while *st < i && s1 > double_gap_i32 + pos1s[*st] {
            *st += 1;
        }
        let lo = if i >= max_iter {
            (i - max_iter).max(*st)
        } else {
            *st
        };

        let mut max_f = match_score;
        let mut max_j = -1i32;
        let mut n_skip = 0usize;
        let mut end_j = lo;
        let strand_i = s1 >> 31;

        'inner: for j in (lo..i).rev() {
            end_j = j;
            let e1 = pos1s[j];
            let e2 = pos2s[j];

            let dist1 = s1 - e1;
            let dist2 = if REV { e2 - s2 } else { s2 - e2 };
            let same_strand = strand_i == (e1 >> 31);
            if dist1 <= 0 || dist2 <= 0 || dist2 > double_gap_i32 || !same_strand {
                continue;
            }

            let gap_penalty = (dist1 - dist2).abs();
            if gap_penalty > max_gap_i32 {
                continue;
            }

            let kmer_overlap_score = dist1.min(dist2).min(match_score);
            let sc = f[j] + kmer_overlap_score - gap_cost * gap_penalty;

            if sc > max_f {
                max_f = sc;
                max_j = j as i32;
                if n_skip > 0 {
                    n_skip -= 1;
                }
            } else if t[j] == i_i32 {
                n_skip += 1;
                if n_skip > max_skip {
                    break 'inner;
                }
            }

            let pj = p[j];
            if pj >= 0 {
                t[pj as usize] = i_i32;
            }
        }

        let rescan = *max_ii < 0 || s1 > double_gap_i32 + pos1s[*max_ii as usize];
        if rescan {
            let mut best = i32::MIN;
            *max_ii = -1;
            for j in (lo..i).rev() {
                if f[j] > best {
                    best = f[j];
                    *max_ii = j as i32;
                }
            }
        }
        let mii = *max_ii;
        if mii >= 0 && (mii as usize) < end_j {
            let m = mii as usize;
            let e1 = pos1s[m];
            let e2 = pos2s[m];
            let dist1 = s1 - e1;
            let dist2 = if REV { e2 - s2 } else { s2 - e2 };
            if dist1 > 0 && dist2 > 0 && dist2 <= double_gap_i32 {
                let gap_penalty = (dist1 - dist2).abs();
                if gap_penalty <= max_gap_i32 {
                    let kmer_overlap_score = dist1.min(dist2).min(match_score);
                    let sc = f[m] + kmer_overlap_score - gap_cost * gap_penalty;
                    if sc > max_f {
                        max_f = sc;
                        max_j = mii;
                    }
                }
            }
        }

        f[i] = max_f;
        p[i] = max_j;

        if *max_ii < 0 || f[*max_ii as usize] < f[i] {
            *max_ii = i_i32;
        }
    }
}

/// Chain defaults matching upstream CompareTwinReadOptions / constants.
pub const CHAIN_MAX_GAP: usize = 200; // MAX_GAP_CHAINING
pub const CHAIN_MAX_SKIP: usize = 10;
pub const CHAIN_DOUBLE_GAP: usize = 10_000;
pub const CHAIN_MIN_LENGTH: usize = 3;

/// Return the single highest-scoring chain (score, ordered anchors, is_reverse),
/// or None if no chain of >= min_chain_length was found. hifiasm already handed
/// us one candidate pair, so we take the best chain and skip the supplementary /
/// secondary multi-mapping filtering that find_optimal_chain applies.
///
/// `band` is myloasm's band (max_mult * 20); `gap_cost` is c (11) by default;
/// `match_score` is 1.
pub fn best_chain(
    matches: &[Anchor],
    gap_cost: i32,
    match_score: i32,
    band: usize,
    max_gap: usize,
    max_skip: usize,
    double_gap: usize,
    min_chain_length: usize,
) -> Option<(i32, Vec<Anchor>, bool)> {
    let chains = dp_anchors_v2(
        matches,
        gap_cost,
        match_score,
        band,
        max_skip,
        max_gap,
        double_gap,
        min_chain_length,
    );
    chains.into_iter().max_by_key(|c| c.0)
}
