
use std::env;

pub fn command_line_str() -> String {
    let args: Vec<String> = env::args().collect();

    args.join(" ")
}