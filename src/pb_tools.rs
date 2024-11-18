use std::time::Duration;

use indicatif::{ProgressBar, ProgressStyle};

pub const DEFAULT_INTERVAL: Duration = Duration::from_millis(500);

#[allow(unused)]
pub fn get_bar_pb(message: String, interval: Duration, total: u64) -> ProgressBar {
    let pb = indicatif::ProgressBar::new(total);
    pb.set_style(indicatif::ProgressStyle::with_template("{msg}::> {wide_bar} {human_pos}/{human_len} speed:{per_sec} elapsed:{elapsed} eta:{eta}").unwrap());
    pb.set_message(message);
    pb.enable_steady_tick(interval);
    pb
}


#[allow(unused)]
pub fn get_spin_pb(message: String, interval: Duration) -> ProgressBar {
    let pb = ProgressBar::new_spinner();
    pb.set_style(ProgressStyle::with_template("{msg} {spinner} {human_pos}  {per_sec} {elapsed}").unwrap());
    pb.enable_steady_tick(interval);
    pb.set_message(message);
    pb
}