use clap::Parser;
use mm2::cli;


fn main() {
    let args = cli::Cli::parse();
    println!("{:?}", args);

}