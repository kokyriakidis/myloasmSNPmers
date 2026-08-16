mod timing;
use timing::PipelineTimer;

use bincode;
use clap::Parser;
use flexi_logger::style;
use flexi_logger::{DeferredNow, Duplicate, FileSpec, Record};
use myloasm::cli;
use myloasm::constants::*;
use myloasm::kmer_comp;
use myloasm::seq_parse;
use myloasm::types;
use myloasm::utils::*;
use std::fs::File;
use std::io::BufReader;
use std::io::IsTerminal;
use std::io::BufWriter;
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;
use sysinfo::System;
use tikv_jemallocator::Jemalloc;

#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

// SNPmer-only fork: this main runs k-mer counting + SNPmer detection and dumps
// the SNPmers, then stops. All downstream assembly (twin reads, overlaps,
// unitig/twin graph, mapping, polishing, dereplication) has been removed.
fn main() {
    let mut args = cli::Cli::parse();

    let output_dir = initialize_setup(&mut args);

    log::info!("Starting SNPmer detection...");
    let mut timer = PipelineTimer::new();

    // Process k-mers, count k-mers, and get SNPmers.
    let kmer_info = timer.measure("get_kmers_and_snpmers", || {
        get_kmers_and_snpmers(&args, &output_dir)
    });
    log_memory_usage(true, "Obtained SNPmers");

    // Dump the detected SNPmers. Each SNPmer is a pair of k-mers identical except
    // at the middle base (position (k-1)/2), both passing the binomial + Fisher
    // strand-balance gates in get_snpmers_inplace_sort. Output is a TSV with the
    // two reconstructed allele k-mers and their strand-summed counts so the set
    // can be inspected / compared against alignment-based het calls.
    {
        use std::io::Write;
        let k = args.kmer_size;
        let decode_kmer = |kmer: u64, k: usize| -> String {
            let mut s = String::with_capacity(k);
            for i in 0..k {
                let val = (kmer >> (2 * i)) & 3;
                s.push(match val {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    3 => 'T',
                    _ => unreachable!(),
                });
            }
            s
        };
        let out_path = Path::new(&args.output_dir).join("snpmers.tsv");
        let mut w = BufWriter::new(
            File::create(&out_path).expect("Could not create snpmers.tsv"),
        );
        writeln!(
            w,
            "split_kmer\tmid_pos\tallele0_base\tallele1_base\tallele0_kmer\tallele1_kmer\tallele0_count\tallele1_count"
        )
        .unwrap();
        let base_char = |b: u8| -> char {
            match b {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                3 => 'T',
                _ => 'N',
            }
        };
        let mid_pos = (k - 1) / 2;
        for s in kmer_info.snpmer_info.iter() {
            let kmer0 = s.split_kmer | ((s.mid_bases[0] as u64) << (k - 1));
            let kmer1 = s.split_kmer | ((s.mid_bases[1] as u64) << (k - 1));
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                s.split_kmer,
                mid_pos,
                base_char(s.mid_bases[0]),
                base_char(s.mid_bases[1]),
                decode_kmer(kmer0, k),
                decode_kmer(kmer1, k),
                s.counts[0],
                s.counts[1],
            )
            .unwrap();
        }
        w.flush().unwrap();
        log::info!(
            "Wrote {} SNPmers to {} (k={}).",
            kmer_info.snpmer_info.len(),
            out_path.display(),
            k
        );
    }
}

fn my_own_format_colored(
    w: &mut dyn std::io::Write,
    now: &mut DeferredNow,
    record: &Record,
) -> Result<(), std::io::Error> {
    let mut paintlevel = record.level();
    if paintlevel == log::Level::Info {
        paintlevel = log::Level::Debug;
    }
    write!(
        w,
        "({}) {} [{}] {}",
        now.format(TS_DASHES_BLANK_COLONS_DOT_BLANK),
        style(paintlevel).paint(record.level().to_string()),
        record.module_path().unwrap_or(""),
        &record.args()
    )
}

fn my_own_format(
    w: &mut dyn std::io::Write,
    now: &mut DeferredNow,
    record: &Record,
) -> Result<(), std::io::Error> {
    write!(
        w,
        "({}) {} [{}] {}",
        now.format(TS_DASHES_BLANK_COLONS_DOT_BLANK),
        record.level(),
        record.module_path().unwrap_or(""),
        &record.args()
    )
}

fn initialize_setup(args: &mut cli::Cli) -> PathBuf {
    if args.markdown_help {
        let markdown_options = clap_markdown::MarkdownOptions::default();
        markdown_options.show_table_of_contents(true);
        clap_markdown::print_help_markdown::<cli::Cli>();
        std::process::exit(0);
    }

    for file in &args.input_files {
        if !Path::new(file).exists() && file != MAGIC_EXIST_STRING {
            eprintln!(
                "ERROR [myloasm] Input file {} does not exist. Exiting.",
                file
            );
            std::process::exit(1);
        }
    }

    let output_dir = Path::new(args.output_dir.as_str());

    if !output_dir.exists() {
        std::fs::create_dir_all(output_dir).expect("Could not create output directory. Exiting.");
    } else {
        if !output_dir.is_dir() {
            eprintln!(
                "ERROR [myloasm] Output directory specified by `-o` exists and is not a directory."
            );
            std::process::exit(1);
        }
    }

    let binary_temp_dir = output_dir.join("binary_temp");
    if !binary_temp_dir.exists() {
        std::fs::create_dir_all(&binary_temp_dir)
            .expect("Could not create temp directory for binary files");
    } else {
        if !binary_temp_dir.is_dir() {
            panic!("Could not create temp directory for binary files. Exiting.");
        }
    }

    // Initialize logger with CLI-specified level
    let log_spec = format!("{},skani=info", args.log_level_filter().to_string());
    let filespec = FileSpec::default()
        .directory(output_dir)
        .basename("myloasm");

    if std::io::stdout().is_terminal() && std::io::stderr().is_terminal() {
        flexi_logger::Logger::try_with_str(log_spec)
            .expect("Something went wrong with logging")
            .log_to_file(filespec) // write logs to file
            .duplicate_to_stderr(Duplicate::Info) // print warnings and errors also to the console
            .format(my_own_format_colored) // use a simple colored format
            .format_for_files(my_own_format)
            .start()
            .expect("Something went wrong with creating log file");
    }

    else{
        flexi_logger::Logger::try_with_str(log_spec)
            .expect("Something went wrong with logging")
            .log_to_file(filespec) // write logs to file
            .duplicate_to_stderr(Duplicate::Info) // print warnings and errors also to the console
            .format(my_own_format) // use a simple colored format
            .format_for_files(my_own_format)
            .start()
            .expect("Something went wrong with creating log file");
    }

    let cli_args: Vec<String> = std::env::args().collect();
    log::info!("COMMAND: {}", cli_args.join(" "));
    log::info!("VERSION: {}", env!("CARGO_PKG_VERSION"));
    log::info!(
        "SYSTEM NAME: {}",
        System::name().unwrap_or(format!("Unknown"))
    );
    log::info!(
        "SYSTEM HOST NAME: {}",
        System::host_name().unwrap_or(format!("Unknown"))
    );
    //log::debug!("BINARY BUILD DATE: {}",  built_info::BUILT_TIME_UTC);
    // The built info is available in the `built` module

    // Validate k-mer size
    if args.kmer_size % 2 == 0 {
        log::error!("K-mer size must be odd");
        std::process::exit(1);
    }
    // Initialize thread pool, bigger stack size because sorting k-mers fails otherwise...
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .stack_size(16 * 1024 * 1024)
        .build_global()
        .unwrap();

    if args.nano_r9 {
        args.snpmer_error_rate_lax = 0.05;
        args.contain_subsample_rate = 2;
        args.kmer_size = 17;
        args.c = 7;
        args.absolute_minimizer_cut_ratio = 50.;
        args.relative_minimizer_cut_ratio = 10.;
        args.min_reads_contig = 2;
    }

    if args.hifi {
        log::info!("HiFi mode enabled. Setting -c to {}.", args.c);
    }

    if let Some(compression) = args.compression{
        args.c = compression;
    }

    return output_dir.to_path_buf();
}

fn get_kmers_and_snpmers(args: &cli::Cli, output_dir: &PathBuf) -> types::KmerGlobalInfo {
    let saved_input = args.input_files == [MAGIC_EXIST_STRING];

    let binary_temp_dir = output_dir.join("binary_temp");
    let snpmer_info_path = binary_temp_dir.join("snpmer_info.bin");

    let kmer_info;
    if saved_input {
        if !snpmer_info_path.exists() {
            log::error!("No input files provided. See --help for usage.");
            std::process::exit(1);
        }
    }

    if saved_input && snpmer_info_path.exists() {
        kmer_info =
            bincode::deserialize_from(BufReader::new(File::open(snpmer_info_path).unwrap()))
                .unwrap();
        log::info!("Loaded snpmer info from file.");
    } else {
        let start = Instant::now();
        let big_kmer_map;
        if args.kmc_db.is_some() {
            log::info!(
                "Using precomputed KMC database at {}",
                args.kmc_db.as_ref().unwrap()
            );
            big_kmer_map = seq_parse::read_kmers_from_kmc_db(
                args.kmer_size,
                args.threads,
                args.kmc_db.as_ref().unwrap(),
                &args,
            );
        } else {
            big_kmer_map = seq_parse::read_to_split_kmers(args.kmer_size, args.threads, &args);
        }
        log::info!(
            "Time elapsed in for counting k-mers is: {:?}",
            start.elapsed()
        );

        let start = Instant::now();
        //kmer_info = kmer_comp::get_snpmers(big_kmer_map, args.kmer_size, &args);
        kmer_info = kmer_comp::get_snpmers_inplace_sort(big_kmer_map, args.kmer_size, &args);
        log::info!(
            "Time elapsed in for parsing snpmers is: {:?}",
            start.elapsed()
        );

        if !args.clean_dir {
            bincode::serialize_into(
                BufWriter::new(File::create(snpmer_info_path).unwrap()),
                &kmer_info,
            )
            .unwrap();
        }
    }
    return kmer_info;
}

