use rayon::ThreadPoolBuilder;
use std::num::NonZero;

/// Initialize the global Rayon thread pool and return the thread count.
///
/// For single core execution this is a no-op.
///
/// Uses `--threads` CLI arg, then `IFX_LOCAL_PROCS` env var, then
/// physical core count via [`num_cpus::get_physical`].
pub fn init_thread_pool(threads: Option<NonZero<usize>>) -> usize {
    let t = if let Some(n) = threads {
        n.get()
    } else if let Some(v) = std::env::var("IFX_LOCAL_PROCS").ok().and_then(|v| v.parse::<usize>().ok()) {
        v
    } else {
        num_cpus::get_physical()
    };

    if t > 1 {
        ThreadPoolBuilder::new().num_threads(t).build_global().unwrap();
    }

    t
}
