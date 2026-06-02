use dais_ribosome::{AnnotationModule, errors::UnimplementedCtype};
use jiff::Timestamp;
use std::collections::{HashMap, HashSet};

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

pub fn print_unimplemented_ctypes(set: HashSet<UnimplementedCtype>, module: &AnnotationModule<'_>) {
    if set.is_empty() {
        return;
    }

    // Get unimplemented ctypes in alphabetical order

    let mut unimplemented_ctypes: Vec<String> = set.clone().into_iter().map(|x| x.0).collect();
    unimplemented_ctypes.sort();

    let mut msg = "no specification yet for:".to_owned();
    for ctype in &unimplemented_ctypes {
        msg.push(' ');
        msg.push_str(ctype);
    }
    ts(&msg);

    // Get unimplemented ctypes grouped by module, in arbitrary order

    let mut ctype_by_module: HashMap<_, Vec<_>> = HashMap::new();
    for ctype in unimplemented_ctypes {
        if let Some(other_module) = module.find_in_other_module(&ctype) {
            ctype_by_module.entry(other_module).or_default().push(ctype);
        }
    }

    if !ctype_by_module.is_empty() {
        for (module, ctypes) in ctype_by_module {
            let mut msg = format!("NOTE, you can use module '{module}' for:");
            for ctype in ctypes {
                msg.push(' ');
                msg.push_str(&ctype);
            }
            ts(&msg);
        }
    }
}
