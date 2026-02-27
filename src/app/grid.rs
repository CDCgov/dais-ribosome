use crate::app::log;
use dais_ribosome::config::current_exe;
use std::{
    env,
    fs::File,
    io::{BufReader, BufWriter},
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

/// Supported grid engine schedulers.
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub(crate) enum GridScheduler {
    Sge,
    Slurm,
}

/// Grid array-task environment parsed from SGE or Slurm.
#[derive(Clone, Debug)]
pub(crate) struct Grid {
    pub task_id:       usize,
    pub task_first:    usize,
    pub task_last:     usize,
    pub task_stepsize: usize,
}

impl Grid {
    /// Parses array-task variables from SGE or Slurm environment.
    ///
    /// Returns an error if no scheduler is detected. Panics if required
    /// variables are invalid.
    pub fn task_vars_from_env() -> std::io::Result<Self> {
        match (env::var("SGE_TASK_ID"), env::var("SLURM_ARRAY_TASK_ID")) {
            (Ok(sge_task_id), _) => Ok(Self::from_sge_env(sge_task_id)),
            (_, Ok(slurm_array_task_id)) => Ok(Self::from_slurm_env(slurm_array_task_id)),
            _ => Err(std::io::Error::other("No supported grid scheduler detected (SGE or Slurm).")),
        }
    }

    fn from_sge_env(sge_task_id: String) -> Self {
        Self {
            task_id:       sge_task_id.parse::<usize>().expect("Invalid SGE_TASK_ID!"),
            task_first:    env::var("SGE_TASK_FIRST")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .unwrap_or(1),
            task_last:     env::var("SGE_TASK_LAST")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .expect("Invalid SGE_TASK_LAST!"),
            task_stepsize: env::var("SGE_TASK_STEPSIZE")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .unwrap_or(1),
        }
    }

    fn from_slurm_env(slurm_array_task_id: String) -> Self {
        Self {
            task_id:       slurm_array_task_id.parse::<usize>().expect("Invalid SLURM_ARRAY_TASK_ID!"),
            task_first:    env::var("SLURM_ARRAY_TASK_MIN")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .unwrap_or(1),
            task_last:     env::var("SLURM_ARRAY_TASK_MAX")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .expect("Invalid SLURM_ARRAY_TASK_MAX!"),
            task_stepsize: 1,
        }
    }
}

/// Picks which grid scheduler to use on the submission node by checking for
/// `qsub` (SGE) or `sbatch` (Slurm) in `PATH`.
///
/// ## Errors
///
/// If neither command exists, then an error is returned.
fn pick_submission_scheduler() -> std::io::Result<GridScheduler> {
    if command_exists("qsub") {
        Ok(GridScheduler::Sge)
    } else if command_exists("sbatch") {
        Ok(GridScheduler::Slurm)
    } else {
        Err(std::io::Error::other(
            "No grid scheduler found (SGE or Slurm). Use local execution only!",
        ))
    }
}

/// Detects whether a command exists by using `which`.
fn command_exists(cmd: &str) -> bool {
    Command::new("which")
        .arg(cmd)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

/// Builds a partition filename of the form `{stem}_{id:03}.{extension}`.
pub(crate) fn get_partition_filename(path: &Path, id: usize, extension: &str) -> String {
    if let Some(stem) = path.file_stem() {
        format!("{stem}_{id:03}.{extension}", stem = stem.display())
    } else {
        format!("ribosome_output_{id:03}.{extension}")
    }
}

/// Submits an array job (`qsub -sync yes` or `sbatch --wait`) and collates
/// partition files into final outputs on completion.
pub fn submit_job_sync(
    task_count: usize, module: &str, input_path: &Path, output_paths: Vec<(PathBuf, &str)>,
) -> std::io::Result<()> {
    log::ts("started, submitting grid job");

    let current_exe = current_exe()?;

    // Log path derived from first output
    let log_path = if let Some((first_out, _)) = output_paths.first() {
        let mut lp = first_out.clone();
        let log_file = lp.file_stem().map_or("ribosome_log.txt".into(), |s| {
            let mut out = s.to_owned();
            out.push("_log.txt");
            out
        });
        lp.set_file_name(log_file);
        lp
    } else {
        PathBuf::from("ribosome_log.txt")
    };

    // Build the base arguments that the child process will receive
    // The child will be invoked with --is-grid-task
    let mut child_args = vec![input_path.to_string_lossy().to_string()];

    // Add output paths as positional args
    for (path, _) in &output_paths {
        child_args.push(path.to_string_lossy().to_string());
    }

    child_args.extend(["--module".to_string(), module.to_string(), "--is-grid-task".to_string()]);

    let output_cmd = match pick_submission_scheduler()? {
        GridScheduler::Sge => {
            let mut cmd = Command::new("qsub");
            cmd.args(["-t", &format!("1-{task_count}"), "-sync", "yes", "-cwd", "-j", "yes", "-o"])
                .arg(&log_path)
                .args(["-b", "yes"])
                .arg(&current_exe)
                .args(&child_args)
                .stdout(Stdio::piped())
                .stderr(Stdio::piped());

            eprintln!("Executing: {cmd:#?}");
            cmd.output()?
        }
        GridScheduler::Slurm => {
            let mut cmd = Command::new("sbatch");
            cmd.args([
                "--wait",
                &format!("--array=1-{task_count}"),
                "--output",
                &log_path.to_string_lossy(),
                "--wrap",
            ])
            .arg(format!(
                "{exe} {args}",
                exe = current_exe.to_string_lossy(),
                args = child_args.join(" "),
            ))
            .stdout(Stdio::piped())
            .stderr(Stdio::piped());

            eprintln!("Executing: {cmd:#?}");
            cmd.output()?
        }
    };

    if !output_cmd.status.success() {
        let stderr = String::from_utf8_lossy(&output_cmd.stderr);
        return Err(std::io::Error::other(stderr));
    }

    log::ts("collating data");

    // Collate partition files into final outputs
    for (output_path, extension) in &output_paths {
        // A copy of the collated output path, which we will mutate to equal
        // each partition output path
        let mut partition_path = output_path.clone();

        let mut writer = BufWriter::new(File::create(output_path)?);

        for id in 1..=task_count {
            partition_path.set_file_name(get_partition_filename(output_path, id, extension));

            // TODO: If it doesn't exist, then we have issues. Should we throw
            // an error? Rather than silently lose data?
            if partition_path.exists() {
                let mut reader = BufReader::new(File::open(&partition_path)?);
                std::io::copy(&mut reader, &mut writer)?;
                std::fs::remove_file(&partition_path)?;
            }
        }
    }

    log::ts("finished");

    Ok(())
}
