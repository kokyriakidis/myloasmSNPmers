// SNPmer-only fork: only the modules on the SNPmer detection path are kept.
// Removed downstream modules: graph, map_processing, mapping, mphmap, polishing,
// polishing_mod, skani_dereplicate, small_genomes, twin_graph, unitig,
// unitig_utils.
pub mod chain;
pub mod cli;
pub mod constants;
pub mod ffi;
pub mod kmc_reader;
pub mod kmer_comp;
pub mod seeding;
pub mod seq_parse;
pub mod types;
pub mod utils;

//pub mod cbloom;
//
//#[cfg(target_arch = "x86_64")]
//pub mod avx2_seeding;
//#[cfg(target_arch = "x86_64")]
//pub mod avx2_chaining;

// Use of a mod or pub mod is not actually necessary.
pub mod built_info {
    // The file has been placed there by the build script.
    // include!(concat!(env!("OUT_DIR"), "/built.rs"));
}
