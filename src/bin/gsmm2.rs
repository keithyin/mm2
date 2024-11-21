use std::thread;

use clap::Parser;
use gskits::fastx_reader::read_fastx;
use gskits::{
    fastx_reader::fasta_reader::FastaFileReader,
    samtools::{samtools_bai, sort_by_coordinates},
};
use mm2::{
    align_worker,
    bam_writer::{write_bam_worker, BamOupArgs},
    build_aligner,
    cli::{self, ReadsToRefAlignArgs},
    query_seq_sender, targets_to_targetsidx,
};

fn alignment(preset: &str, align_threads: Option<usize>, args: &ReadsToRefAlignArgs) {
    let target_filename = args
        .io_args
        .target
        .as_ref()
        .expect("target need to be provided");

    let fa_iter = FastaFileReader::new(target_filename.to_string());
    let targets = read_fastx(fa_iter);
    let target2idx = targets_to_targetsidx(&targets);

    let aligners = build_aligner(
        preset,
        &args.index_args,
        &args.map_args,
        &args.align_args,
        &args.oup_args,
        &targets,
    );

    /*
    1. query_seq_sender
    2. align
    3. write to bam
    */
    thread::scope(|s| {
        let aligners = &aligners;
        let target2idx = &target2idx;

        let (qs_sender, qs_recv) = crossbeam::channel::bounded(1000);
        s.spawn(move || {
            query_seq_sender(&args.io_args.query, qs_sender);
        });

        let num_threads = align_threads.unwrap_or(num_cpus::get_physical());
        let (align_res_sender, align_res_recv) = crossbeam::channel::bounded(1000);
        for _ in 0..num_threads {
            let qs_recv_ = qs_recv.clone();
            let align_res_sender_ = align_res_sender.clone();
            s.spawn(move || align_worker(qs_recv_, align_res_sender_, aligners, target2idx));
        }
        drop(qs_recv);
        drop(align_res_sender);

        write_bam_worker(
            align_res_recv,
            target2idx,
            &args.io_args.get_oup_path(),
            &BamOupArgs::from(&args.oup_args),
            "gsmm2",
            env!("CARGO_PKG_VERSION"),
            true,
        );
    });
    sort_by_coordinates(&args.io_args.get_oup_path(), align_threads);
    samtools_bai(&args.io_args.get_oup_path(), true, align_threads).unwrap();
}

fn main() {
    let args = cli::Cli::parse();

    let preset = &args.preset;
    let align_threads = args.threads.clone();

    match args.commands {
        cli::Commands::Align(ref args) => {
            alignment(preset, align_threads, args);
        }

        _ => panic!("not implemented yet"),
    }
}
