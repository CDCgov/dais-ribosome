use dais_ribosome::annotation::AnnotationModule;
use jiff::Timestamp;
use std::collections::HashSet;

const PROG: &str = "RIBOSOME";

/// Traditional IFX logger formatting
///
/// Modified from SSWSORT
pub fn time_stamp(message: &str, use_stderr: bool) {
    let shlvl = std::env::var("SHLVL").ok().and_then(|v| v.parse::<usize>().ok()).unwrap_or(0);
    let pad = "  ".repeat(shlvl.saturating_sub(1));

    if use_stderr {
        eprintln!(
            "[{now}] {pad}{PROG} :: {message}",
            now = Timestamp::now().strftime("%Y-%m-%d %k:%M:%S")
        );
    } else {
        println!(
            "[{now}] {pad}{PROG} :: {message}",
            now = Timestamp::now().strftime("%Y-%m-%d %k:%M:%S")
        );
    }
}

pub fn ts(message: &str) {
    time_stamp(message, false);
}

pub fn print_unimplemented_ctypes(set: HashSet<String>, module: &AnnotationModule<'_>) {
    if set.is_empty() {
        return;
    }

    let mut unimplemented_ctypes: Vec<String> = set.clone().into_iter().collect();
    unimplemented_ctypes.sort();

    let mut msg = "no specification yet for:".to_owned();
    for ctype in unimplemented_ctypes {
        msg.push(' ');
        msg.push_str(&ctype);
    }
    ts(&msg);

    let found = module.suggest_modules_for_compound_types(set);

    if !found.is_empty() {
        for (module, ctypes) in found {
            let mut msg = format!("NOTE, you can use module '{module}' for:");
            for ctype in ctypes {
                msg.push(' ');
                msg.push_str(&ctype);
            }
            ts(&msg);
        }
    }
}
