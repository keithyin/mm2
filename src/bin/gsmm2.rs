use std::{collections::HashMap, thread};

use clap::Parser;
use mm2::{
    aligner, bam_writer, build_aligner,
    cli::{self, ReadsToRefAlignArgs},
    fille_reader::read_fasta,
    query_seq_sender,
};

fn alignment(preset: &str, args: &ReadsToRefAlignArgs) {
    let target_filename = args
        .io_args
        .target
        .as_ref()
        .expect("target need to be provided");

    let targets = read_fasta(target_filename).unwrap();

    let mut target2idx = HashMap::new();
    targets.iter().enumerate().for_each(|(idx, target)| {
        target2idx.insert(target.qname.clone(), (idx, target.sequence.len()));
    });

    let aligners = build_aligner(
        preset,
        &args.index_args,
        &args.map_args,
        &args.align_args,
        &targets,
    );

    thread::scope(|s| {
        let aligners = &aligners;
        let target2idx = &target2idx;

        let (qs_sender, qs_recv) = crossbeam::channel::bounded(1000);
        s.spawn(move || {
            query_seq_sender(&args.io_args.query, qs_sender);
        });

        let num_threads = num_cpus::get_physical();
        let (align_res_sender, align_res_recv) = crossbeam::channel::bounded(1000);
        for _ in 0..num_threads {
            let qs_recv_ = qs_recv.clone();
            let align_res_sender_ = align_res_sender.clone();
            s.spawn(move || aligner(qs_recv_, align_res_sender_, aligners));
        }
        drop(qs_recv);
        drop(align_res_sender);

        bam_writer(align_res_recv, target2idx, &args.io_args.get_oup_path());
    });
}

fn main() {
    let args = cli::Cli::parse();

    let preset = &args.preset;

    match args.commands {
        cli::Commands::R2R(ref args) => {
            alignment(preset, args);
        }

        _ => panic!("not implemented yet"),
    }
}
