//! Code for handling possible grid execution of a binary.
//!
//! This code assumes that `clap` is being used for argument parsing, and that
//! SGE is used for grid execution.
//!
//! ## Tutorial
//!
//! 1. Define the CLI for the binary with `clap`, implementing/deriving
//!    [`Parser`] on them as normal. Add the arguments `is_grid_task: bool` and
//!    `submit_grid_job: Option<usize>` for indicating whether the execution is
//!    part of an existing grid job, and whether grid execution is requested
//!    (with the given number of tasks). These arguments can also be renamed if
//!    desired.
//! 2. Implement [`GridCompatibleCli`] on the struct, implementing the relevant
//!    getter functions and possibly overriding the associated constants if the
//!    arguments were renamed.
//! 3. If needed, define a parsed arguments struct containing any
//!    transformations of the raw CLI. This transformed struct should contain
//!    all the output paths, but should not create any writers yet. If no
//!    transformations or validation need to occur, then skip this step.
//! 4. Implement [`GridCompatibleArgs`] on the parsed arguments struct. This
//!    includes getters for the output paths, as well as the logic for
//!    converting from the CLI to the parsed arguments. If no parsed arguments
//!    struct is being used, then implement [`GridCompatibleArgs`] on the CLI,
//!    with no logic performed in [`from_cli`].
//! 5. In your binary, parse the arguments using [`parse_maybe_grid`] or
//!    [`parse_from_maybe_grid`] (instead of clap's [`parse`] and
//!    [`parse_from`]). In addition to parsing, these methods automatically
//!    fetch the relevant grid information, and they also mutate the output
//!    paths to include the task ID if applicable. The parsed arguments and any
//!    grid information are returned in a tuple.
//! 6. After performing any other validation or logic, handle the
//!    [`GridInfo::Requested`] variant by calling [`submit_job_sync`].
//! 7. To handle the [`GridInfo::Task`] variant (for an actively running task),
//!    call [`GridTask::run_task`], passing all processing logic as a closure.
//!    [`GridTaskInfo::select_inputs`] can be used within the closure to select
//!    the items to handle from an iterator over the inputs.
//! 8. If the returned grid information was `None`, perform processing as
//!    normal.
//!
//! [`from_cli`]: GridCompatibleArgs::from_cli
//! [`parse_maybe_grid`]: GridCompatibleArgs::parse_maybe_grid
//! [`parse_from_maybe_grid`]: GridCompatibleArgs::parse_from_maybe_grid
//! [`parse`]: Parser::parse
//! [`parse_from`]: Parser::parse_from
//! [`submit_job_sync`]: GridRequestedInfo::submit_job_sync`

use crate::log;
use clap::{Arg, ArgAction, ArgMatches, CommandFactory, FromArgMatches, Parser, parser::ValueSource};
use std::{
    cell::OnceCell,
    collections::{BTreeMap, btree_map::Entry},
    env::{self, current_exe},
    ffi::{OsStr, OsString},
    fs::{File, read_to_string},
    io::{BufReader, BufWriter, ErrorKind, Write},
    path::{Path, PathBuf},
    process::{Command, Stdio},
};
use zoe::{
    data::err::{ErrorWithContext, Fail, ResultWithErrorContext, WithErrorContext},
    prelude::OrFail,
};

/// An extension trait for [`Parser`] enabling the grid-execution arguments to
/// be detected.
///
/// This requires the following two arguments to be present in the CLI:
///
/// - The argument ID `is_grid_task`, which holds a `bool` for whether or not
///   the current execution is part of an array job.
/// - The argument ID `submit_grid_job`, which holds an `Option<usize>`. This
///   signals that an array job is requested, and takes the size of the array
///   job as an argument.
///
/// Naming the arguments `is_grid_task` and `submit_grid_job` ensures the IDs
/// are present. If a different name is used, then the constants
/// [`IS_GRID_TASK_ID`] and [`SUBMIT_GRID_JOB_ID`] can be overridden.
///
/// [`IS_GRID_TASK_ID`]: GridCompatibleCli::IS_GRID_TASK_ID
/// [`SUBMIT_GRID_JOB_ID`]: GridCompatibleCli::SUBMIT_GRID_JOB_ID
pub trait GridCompatibleCli: Parser + Sized {
    /// The clap ID for indicating that the task is running as part of a grid
    /// task.
    const IS_GRID_TASK_ID: &'static str = "is_grid_task";

    /// The clap ID for indicating that the task should be run on the grid as an
    /// array job with the specified size.
    const SUBMIT_GRID_JOB_ID: &'static str = "submit_grid_job";

    /// A getter function for the parsed argument corresponding to
    /// [`IS_GRID_TASK_ID`].
    ///
    /// [`IS_GRID_TASK_ID`]: GridCompatibleCli::IS_GRID_TASK_ID
    fn is_grid_task(&self) -> bool;

    /// A getter function for the parsed argument corresponding to
    /// [`SUBMIT_GRID_JOB_ID`].
    ///
    /// [`SUBMIT_GRID_JOB_ID`]: GridCompatibleCli::SUBMIT_GRID_JOB_ID
    fn submit_grid_job(&self) -> Option<usize>;
}

/// A trait to be implemented on a parsed arguments struct, which can be
/// constructed from a [`GridCompatibleCli`] and any extracted [`GridInfo`].
///
/// The [`GridCompatibleArgs`] must contain all the output paths, and an
/// iterator over them must be provided under [`outputs`]. When parsing and
/// collating, task IDs are incorporated into these paths using [`add_id`] and
/// [`add_id_tmp`], which have defaults provided.
///
/// The [`GridCompatibleArgs`] should not contain any writers, since the output
/// paths may be mutated.
///
/// [`outputs`]: GridCompatibleArgs::outputs
/// [`add_id`]: GridCompatibleArgs::add_id
/// [`add_id_tmp`]: GridCompatibleArgs::add_id_tmp
pub trait GridCompatibleArgs: Sized {
    /// The underlying CLI used to form the parsed arguments.
    type Cli: GridCompatibleCli;

    /// Constructs the parsed arguments from the CLI.
    ///
    /// After implementing this, [`parse_maybe_grid`] and
    /// [`parse_from_maybe_grid`] can be used similarly to [`Parser::parse`] and
    /// [`Parser::parse_from`].
    ///
    /// This function must ensure that all the output paths are stored in a way
    /// that is accessible to [`outputs`]. This function should not open any
    /// writers, since the output paths may be mutated to include task IDs.
    ///
    /// [`outputs`]: GridCompatibleArgs::outputs
    /// [`parse_maybe_grid`]: GridCompatibleArgs::parse_maybe_grid
    /// [`parse_from_maybe_grid`]: GridCompatibleArgs::parse_from_maybe_grid
    fn from_cli(cli: Self::Cli) -> std::io::Result<Self>;

    /// Similar to [`from_cli`], but provides a mutable reference to the
    /// [`ArgMatches`] generated by clap, so that they can be edited.
    ///
    /// For example, suppose that a command line utility has an option for a
    /// random seed, and that grid submission requires that the seed is the same
    /// for all tasks. Then upon confirming that a seed was not set in `cli`,
    /// [`from_cli_with_matches`] can inject the flag as if it were provided by
    /// the user with:
    ///
    /// ```rust
    /// let command = Command::new("temp").arg(Arg::new("seed").long("seed").num_args(1));
    /// matches.push(command.get_matches_from(["temp", "--seed", "42"]));
    /// ```
    ///
    /// If this is overridden, then [`from_cli`] can be implemented to call to
    /// this with an empty vector for `matches`.
    ///
    /// [`from_cli`]: GridCompatibleArgs::from_cli
    /// [`from_cli_with_matches`]: GridCompatibleArgs::from_cli_with_matches
    fn from_cli_with_matches(cli: Self::Cli, matches: &mut Vec<ArgMatches>) -> std::io::Result<Self> {
        let _ = matches;
        Self::from_cli(cli)
    }

    /// The output files that need to be collated.
    ///
    /// This function must not mutate `self`. Instead, it should directly return
    /// mutable references to the output paths.
    fn outputs_mut(&mut self) -> impl Iterator<Item = &mut PathBuf>;

    /// The path for where to write the `stdout` and `stderr` streams of the
    /// job.
    ///
    /// The task ID will be added to this path using [`add_id`] prior to
    /// submission.
    ///
    /// ## Errors
    ///
    /// See the implementor for any specific errors. Since path manipulation can
    /// be fallible, this method is fallible to prevent implementors from
    /// needing to panic.
    ///
    /// [`add_id`]: GridCompatibleArgs::add_id
    fn log_path(&self) -> std::io::Result<PathBuf>;

    /// Adds the specified `id` to the `path` to get the output file for a
    /// particular task.
    ///
    /// By default, this adds the ID before the extension, but this can be
    /// overridden.
    #[must_use]
    fn add_id(path: &Path, id: &str) -> PathBuf {
        let mut path = path.to_path_buf();

        let Some(extension) = path.extension().map(|s| s.to_os_string()) else {
            let Some(file_name) = path.file_name() else {
                return path.join(id);
            };

            let mut file_name = file_name.to_os_string();
            file_name.push("_");
            file_name.push(id);
            path.set_file_name(file_name);
            return path;
        };

        path.set_extension("");

        let Some(file_name) = path.file_name() else {
            // This should be unreachable, but let's handle it anyways
            return path.join(id).with_extension(extension);
        };

        let mut file_name = file_name.to_os_string();
        file_name.push("_");
        file_name.push(id);
        path.set_file_name(file_name);
        path.add_extension(extension);

        path
    }

    /// Adds the specified `id` to the `path` to get the output file for a
    /// particular task, as well as altering the path to indicate that the file
    /// has not finished being written to.
    ///
    /// By default, this uses [`add_id`] then adds a `.partial` suffix,
    /// but this can be overridden.
    ///
    /// [`add_id`]: GridCompatibleArgs::add_id
    fn add_id_tmp(path: &Path, id: &str) -> PathBuf {
        let mut path = Self::add_id(path, id);
        path.add_extension("partial");
        path
    }

    /// Similar to [`Parser::parse`], but also fetches the relevant grid
    /// execution information, as well as mutating any output paths to include
    /// the task ID if appropriate.
    ///
    /// ## Errors
    ///
    /// In the case of grid execution, an error is returned if no scheduler is
    /// detected or the required environmental variables are unset or invalid.
    /// See [`from_cli`] for any other errors.
    ///
    /// [`from_cli`]: GridCompatibleArgs::from_cli
    fn parse_maybe_grid() -> std::io::Result<(Self, Option<GridInfo>)> {
        let GridParser {
            parser: cli,
            mut arg_matches,
        } = GridParser::<Self::Cli>::parse();

        let is_grid_task = cli.is_grid_task();
        let submit_grid_job = cli.submit_grid_job();

        let mut args = Self::from_cli_with_matches(cli, &mut arg_matches)?;

        let grid_info = GridInfo::new(&mut args, &arg_matches, is_grid_task, submit_grid_job)?;

        Ok((args, grid_info))
    }

    /// Similar to [`Parser::parse_from`], but also fetches the relevant grid
    /// execution information, as well as mutating any output paths to include
    /// the task ID if appropriate.
    ///
    /// ## Errors
    ///
    /// In the case of grid execution, an error is returned if no scheduler is
    /// detected or the required environmental variables are unset or invalid.
    /// See [`from_cli`] for any other errors.
    ///
    /// [`from_cli`]: GridCompatibleArgs::from_cli
    fn parse_from_maybe_grid<I, T>(itr: I) -> std::io::Result<(Self, Option<GridInfo>)>
    where
        I: IntoIterator<Item = T>,
        T: Into<OsString> + Clone, {
        let GridParser {
            parser: cli,
            mut arg_matches,
        } = GridParser::<Self::Cli>::parse_from(itr);

        let is_grid_task = cli.is_grid_task();
        let submit_grid_job = cli.submit_grid_job();

        let mut args = Self::from_cli_with_matches(cli, &mut arg_matches)?;

        let grid_info = GridInfo::new(&mut args, &arg_matches, is_grid_task, submit_grid_job)?;

        Ok((args, grid_info))
    }
}

impl GridTaskInfo {
    /// Parses array-task variables from SGE environment variables.
    ///
    /// ## Errors
    ///
    /// An error is returned if no scheduler is detected or the required
    /// environmental variables are unset or invalid.
    fn new<T>(args: &mut T) -> std::io::Result<Self>
    where
        T: GridCompatibleArgs, {
        match env::var("SGE_TASK_ID") {
            Ok(sge_task_id) => Self::from_sge_env(sge_task_id, args),
            _ => Err(std::io::Error::other("No supported grid scheduler detected (SGE).")),
        }
    }

    /// Constructs a [`GridTaskInfo`] from SGE environmental variables.
    ///
    /// ## Errors
    ///
    /// `SGE_TASK_ID` and `SGE_TASK_LAST` must be set and be valid.
    fn from_sge_env<T>(sge_task_id: String, args: &mut T) -> std::io::Result<Self>
    where
        T: GridCompatibleArgs, {
        let task_id = sge_task_id
            .parse::<usize>()
            .with_context(format!("Invalid SGE_TASK_ID: {sge_task_id}"))?;

        let task_first = env::var("SGE_TASK_FIRST")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(1);

        let task_last = env::var("SGE_TASK_LAST").with_context("Invalid SGE_TASK_LAST")?;
        let task_last = task_last
            .parse::<usize>()
            .with_context(format!("Invalid SGE_TASK_LAST: {task_last}"))?;

        let task_stepsize = env::var("SGE_TASK_STEPSIZE")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(1);

        Ok(Self {
            task_id,
            task_first,
            task_last,
            task_stepsize,
            scheduler: GridScheduler::Sge,
            output_paths: args.outputs_mut().map(|x| (*x).clone()).collect(),
            add_id: T::add_id,
            add_id_tmp: T::add_id_tmp,
        })
    }
}

/// Picks which grid scheduler to use on the submission node by checking for
/// `qsub` (SGE) in `PATH`.
///
/// ## Errors
///
/// If the `qsub` command does not exist, then an error is returned.
fn pick_submission_scheduler() -> std::io::Result<GridScheduler> {
    if command_exists("qsub") {
        Ok(GridScheduler::Sge)
    } else {
        Err(std::io::Error::other("SGE not found. Use local execution only!"))
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

/// Information related to a request to use grid submission, which can be used
/// to submit the job.
pub struct GridRequestedInfo {
    /// The command for re-running the binary with the same arguments.
    command:      Vec<OsString>,
    /// The number of tasks to use when submitting the job.
    task_count:   usize,
    /// The scheduler to use for the grid job.
    scheduler:    GridScheduler,
    /// The output paths that will be formed via collation.
    output_paths: Vec<PathBuf>,
    /// The path to use for the scheduler's log file, without the placeholder
    /// for the ID added yet.
    log_path:     PathBuf,
    /// A pointer to [`GridCompatibleArgs::add_id`].
    add_id:       fn(&Path, &str) -> PathBuf,
    /// A pointer to [`GridCompatibleArgs::add_id_tmp`].
    add_id_tmp:   fn(&Path, &str) -> PathBuf,
}

impl GridRequestedInfo {
    /// Submits the job as an array task.
    ///
    /// This function will block on the completion of the job (corresponding to
    /// `-sync yes` in [`Sge`]). The job executes in the current working
    /// directory (`-cwd` in [`Sge`]). The `stdout` and `stderr` stream
    /// generated by the scheduler are merged and written to [`log_path`].
    ///
    /// ## Errors
    ///
    /// A custom error is returned with as much context as possible given the
    /// failure conditions, designed to be displayed directly using
    /// [`JobErrorOrFail::unwrap_or_fail`]. Any of the following can trigger an
    /// error:
    ///
    /// - The job submission command failing
    /// - Any of the tasks not running to completion, as evidenced by the output
    ///   files not being renamed properly
    /// - Any of the tasks having missing output files
    /// - Collation failing (e.g., due to IO errors)
    ///
    /// [`Sge`]: GridScheduler::Sge
    /// [`log_path`]: GridCompatibleArgs::log_path
    /// [`ErrorKind::Other`]: std::io::ErrorKind::Other
    pub fn submit_job_sync(&self) -> Result<(), SubmitError> {
        log::ts("started, submitting grid job");

        // For unclear reasons, when the output directories do not exist yet,
        // there is sometimes a race condition causing collation to fail even
        // though the directory ends up existing by the end of the program. To
        // be safe, create all directories now.
        for output in &self.output_paths {
            if let Some(parent) = output.parent() {
                std::fs::create_dir_all(parent)
                    .with_path_context("Failed to create output directory", parent)
                    .map_err(|e| SubmitErrorRepr::Io(e.into()))?
            }
        }

        let add_id = self.add_id;

        let output_cmd = match self.scheduler {
            GridScheduler::Sge => {
                let log_path = add_id(&self.log_path, "$TASK_ID");

                let mut cmd = Command::new("qsub");

                cmd.args([
                    "-t",
                    &format!("1-{task_count}", task_count = self.task_count),
                    "-sync",
                    "yes",
                    "-cwd",
                    "-j",
                    "yes",
                    "-o",
                ])
                .arg(log_path)
                .args(["-b", "yes"])
                .args(&self.command)
                .stdout(Stdio::piped())
                .stderr(Stdio::piped());

                self.print_submission_msg(self.task_count, &self.log_path, add_id);

                cmd.output()
                    .with_context("An error using qsub occurred")
                    .map_err(|e| SubmitErrorRepr::Io(e.into()))?
            }
        };

        let command_result = if !output_cmd.status.success() {
            let stderr = String::from_utf8_lossy(&output_cmd.stderr).trim().to_owned();
            if stderr.is_empty() { Err(None) } else { Err(Some(stderr)) }
        } else {
            Ok(())
        };

        // The exit code to return if there is an error (ensured to be non-zero)
        let exit_code = output_cmd.status.code().unwrap_or(1).max(1);

        log::ts("collating data");

        // Perform collation, accumulating task errors and the first IO error
        let mut collation_task_errors = TaskErrors::default();
        let mut collation_io_error = OnceCell::new();

        for output in &self.output_paths {
            match self.collate_output(output) {
                Ok(()) => {}
                Err(CollationError::Io(e)) => {
                    collation_io_error.get_or_init(|| {
                        SubmitErrorRepr::Io(
                            e.with_path_context("Failed to collate outputs for file", output.as_path())
                                .into(),
                        )
                    });
                }
                Err(CollationError::Task(e)) => {
                    collation_task_errors.merge(e);
                }
            }
        }

        // If we have confirmed that there is a specific task that failed, build
        // that error, which may involve slurping a log file.
        let task_error = collation_task_errors.to_string().map(|task_errors_str| {
            let command_result = command_result.clone();

            if let Some((first_task_failed, _)) = collation_task_errors
                .by_id
                .iter()
                .find(|(_, cause)| matches!(cause, TaskErrorCause::TaskFailed))
                && let log_path = add_id(&self.log_path, &format!("{first_task_failed}"))
                && let Ok(slurped_log) = read_to_string(log_path)
            {
                let slurped_log = slurped_log.trim();
                JobError {
                    command_error: command_result.err().flatten(),
                    task_errors: Some(task_errors_str),
                    log_mention: Some(format!("See log for task {first_task_failed}:\n\n{slurped_log}")),
                    exit_code,
                }
            } else {
                JobError {
                    command_error: command_result.err().flatten(),
                    task_errors: Some(task_errors_str),
                    log_mention: None,
                    exit_code,
                }
            }
        });

        // Collate the log files
        let log_result = self
            .collate_log()
            .with_path_context("Failed to collate log files into file", &self.log_path);

        // Propagate task error, adding missing log context if log collation
        // succeeded
        if let Some(mut task_error) = task_error {
            if log_result.is_ok() {
                task_error
                    .log_mention
                    .get_or_insert(format!("The cause may be in {}", self.log_path.display()));
            }
            return Err(task_error.into());
        }

        // Propoagate command error
        if let Err(command_error) = command_result {
            return Err(JobError {
                command_error,
                task_errors: None,
                log_mention: log_result
                    .is_ok()
                    .then(|| format!("The cause may be in {}", self.log_path.display())),
                exit_code,
            }
            .into());
        }

        // Propoagate collation IO error
        if let Some(collation_io_error) = collation_io_error.take() {
            return Err(collation_io_error.into());
        }

        // Propagate log collation error
        log_result.map_err(|e| SubmitErrorRepr::Io(e.into()))?;

        log::ts("finished");

        Ok(())
    }

    /// Prints the confirmation message that a grid job has been submitted,
    /// listing the number of tasks and the log file.
    fn print_submission_msg(&self, tasks: usize, log_path: &Path, add_id: fn(&Path, &str) -> PathBuf) {
        let log_path = add_id(log_path, "<ID>");
        log::ts(&format!(
            "submitted synchronous job with {tasks} tasks, log files at: '{log_path}'",
            log_path = log_path.display()
        ));
    }

    /// Collate the output for the given path (without the task ID added).
    ///
    /// ## Errors
    ///
    /// IO errors creating the collated file or reading the files for each task
    /// are propagated as [`CollationError::Io`]. If an output file for a task
    /// is not found (or was not renamed by the task to indicate successful
    /// completion), a [`CollationError::Task`] error is returned.
    fn collate_output(&self, output: &Path) -> Result<(), CollationError> {
        let mut writer = BufWriter::new(File::create(output).with_context("Failed to create collated file")?);

        let mut task_errors = BTreeMap::new();

        let add_id = self.add_id;
        let add_id_tmp = self.add_id_tmp;

        for id in 1..=self.task_count {
            let id_str = format!("{id}");

            let path = add_id(output, &id_str);

            let file = match File::open(&path) {
                Ok(file) => file,
                Err(e) if e.kind() == ErrorKind::NotFound => {
                    let tmp_file = add_id_tmp(output, &id_str);

                    let cause = if tmp_file.exists() {
                        TaskErrorCause::TaskFailed
                    } else {
                        TaskErrorCause::MissingOutputFile(path.clone())
                    };

                    task_errors.insert(id, cause);

                    continue;
                }
                Err(e) => {
                    return Err(e
                        .with_path_context(format!("Failed to open file for task {id}"), &path)
                        .into());
                }
            };

            let mut reader = BufReader::new(file);
            std::io::copy(&mut reader, &mut writer)
                .with_path_context(format!("Failed to copy output from file for task {id}"), &path)?;
            std::fs::remove_file(&path).with_path_context(format!("Failed to remove file for task {id}"), &path)?;
        }

        if task_errors.is_empty() {
            Ok(())
        } else {
            Err(CollationError::Task(TaskErrors { by_id: task_errors }))
        }
    }

    /// Collates the log files for each task.
    ///
    /// ## Errors
    ///
    /// IO errors creating the collated log file or reading the log files for
    /// each task are propagated. Missing log files are ignored.
    fn collate_log(&self) -> std::io::Result<()> {
        let mut some_printed = false;

        let mut writer = BufWriter::new(File::create(&self.log_path)?);

        let add_id = self.add_id;

        // TODO: Different starting index? Step?
        for id in 1..=self.task_count {
            let path = add_id(&self.log_path, &format!("{id}"));

            Self::write_header(&mut writer, id, &mut some_printed)?;

            let log_file = match File::open(&path) {
                Ok(file) => file,
                Err(e) if e.kind() == ErrorKind::NotFound => {
                    continue;
                }
                Err(e) => {
                    return Err(e.with_path_context(format!("Failed to open log for task {id}"), path).into());
                }
            };

            let mut reader = BufReader::new(log_file);

            std::io::copy(&mut reader, &mut writer)
                .with_path_context(format!("Failed to copy output from log for task {id}"), &path)?;

            std::fs::remove_file(&path).with_path_context(format!("Failed to remove log for task {id}"), &path)?;
        }

        Ok(())
    }

    /// Helper function for [`collate_log`] that writes a header to the log file
    /// for the given `id`.
    ///
    /// [`collate_log`]: GridRequestedInfo::collate_log
    fn write_header(writer: &mut BufWriter<File>, id: usize, some_printed: &mut bool) -> std::io::Result<()> {
        const HEADER: &str = "Output for task ID: ";

        let task_id_len = id.checked_ilog10().unwrap_or(0) as usize + 1;
        let header_len = HEADER.len() + task_id_len;

        if *some_printed {
            writeln!(writer)?;
        } else {
            *some_printed = true;
        }
        writeln!(writer, "{empty:=<header_len$}", empty = "")?;
        writeln!(writer, "{HEADER}{id}")?;
        writeln!(writer, "{empty:=<header_len$}", empty = "")?;
        writeln!(writer)
    }
}

pub struct SubmitError {
    inner: SubmitErrorRepr,
}

enum SubmitErrorRepr {
    Io(std::io::Error),
    Job(JobError),
}

impl From<SubmitErrorRepr> for SubmitError {
    fn from(inner: SubmitErrorRepr) -> Self {
        SubmitError { inner }
    }
}

impl From<JobError> for SubmitError {
    fn from(value: JobError) -> Self {
        SubmitError {
            inner: SubmitErrorRepr::Job(value),
        }
    }
}

struct JobError {
    command_error: Option<String>,
    task_errors:   Option<String>,
    log_mention:   Option<String>,
    exit_code:     i32,
}

impl SubmitError {
    pub fn fail(self) -> ! {
        match self.inner {
            SubmitErrorRepr::Io(err) => err.fail(),
            SubmitErrorRepr::Job(err) => {
                if let Ok(bin) = std::env::current_exe() {
                    eprintln!("Error in {b}:\n", b = bin.display());
                }

                eprintln!("The grid job failed");

                if let Some(command_error) = err.command_error {
                    eprintln!("\nError output from scheduler:\n\n{command_error}");
                }

                if let Some(task_errors) = err.task_errors {
                    eprintln!("\nThe following tasks encountered problems:\n\n{task_errors}");
                }

                if let Some(log_mention) = err.log_mention {
                    eprintln!("\n{log_mention}");
                }

                std::process::exit(err.exit_code);
            }
        }
    }
}

pub trait JobErrorOrFail<T> {
    fn unwrap_or_fail(self) -> T;
}

impl<T> JobErrorOrFail<T> for Result<T, SubmitError> {
    fn unwrap_or_fail(self) -> T {
        match self {
            Ok(val) => val,
            Err(e) => e.fail(),
        }
    }
}

pub enum TaskErrorCause {
    MissingOutputFile(PathBuf),
    TaskFailed,
}

#[derive(Default)]
pub struct TaskErrors {
    by_id: BTreeMap<usize, TaskErrorCause>,
}

impl TaskErrors {
    fn to_string(&self) -> Option<String> {
        const TASK_HEADER: &str = "Task";
        const STATUS_HEADER: &str = "Status";
        const PATH_HEADER: &str = "Path";
        const MISSING_OUTPUT_VAL: &str = "Missing output file";
        const TASK_FAILED_VAL: &str = "Task failed";
        const PADDING: &str = "    ";

        if self.by_id.is_empty() {
            return None;
        }

        let mut num_task_failed = 0;
        let mut num_missing_output = 0;

        for cause in self.by_id.values() {
            match cause {
                TaskErrorCause::MissingOutputFile(_) => num_missing_output += 1,
                TaskErrorCause::TaskFailed => num_task_failed += 1,
            }
        }

        let task_width = self.by_id.keys().max().copied().unwrap_or(0).max(TASK_HEADER.len());
        let status_width = if num_missing_output == 0 {
            STATUS_HEADER.len().max(TASK_FAILED_VAL.len())
        } else if num_task_failed == 0 {
            STATUS_HEADER.len().max(MISSING_OUTPUT_VAL.len())
        } else {
            STATUS_HEADER.len().max(TASK_FAILED_VAL.len()).max(MISSING_OUTPUT_VAL.len())
        };

        let mut out = format!("{TASK_HEADER:<task_width$}{PADDING}{STATUS_HEADER:<status_width$}");
        if num_missing_output > 0 {
            out.push_str(&format!("{PADDING}{PATH_HEADER}"));
        }
        out.push('\n');

        for (i, (id, cause)) in self.by_id.iter().enumerate() {
            out.push_str(&format!("{id:<task_width$}{PADDING}"));

            match cause {
                TaskErrorCause::MissingOutputFile(path) => {
                    out.push_str(&format!(
                        "{MISSING_OUTPUT_VAL:<status_width$}{PADDING}{path}",
                        path = path.display()
                    ));
                }
                TaskErrorCause::TaskFailed => {
                    out.push_str(&format!("{TASK_FAILED_VAL:<status_width$}"));
                }
            }

            if i < self.by_id.len() - 1 {
                out.push('\n');
            }
        }

        Some(out)
    }

    fn merge(&mut self, other: TaskErrors) {
        for (k, v) in other.by_id {
            match self.by_id.entry(k) {
                Entry::Occupied(mut entry) => {
                    if matches!(entry.get(), TaskErrorCause::TaskFailed) {
                        continue;
                    }
                    if matches!(v, TaskErrorCause::TaskFailed) {
                        entry.insert(v);
                    }
                }
                Entry::Vacant(entry) => {
                    entry.insert(v);
                }
            }
        }
    }
}

pub enum CollationError {
    Task(TaskErrors),
    Io(std::io::Error),
}

impl From<std::io::Error> for CollationError {
    fn from(value: std::io::Error) -> Self {
        CollationError::Io(value)
    }
}

impl From<ErrorWithContext> for CollationError {
    fn from(value: ErrorWithContext) -> Self {
        CollationError::Io(value.into())
    }
}

/// A wrapper around [`GridTaskInfo`] to ensure proper handling.
///
/// When inside a grid task (i.e., [`GridInfo::Task`] is present), the following
/// steps should be used:
///
/// 1. Match on the [`GridInfo`] to extract the [`GridTask`] (i.e., confirm that
///    we are inside a grid task).
/// 2. Call [`GridTask::run_task`], with all the processing logic inside the
///    passed closure. The closure takes one argument, which is the inner
///    [`GridTaskInfo`].
/// 3. Within the closure, select the inputs pertaining to the given task using
///    [`GridTaskInfo::select_inputs`].
/// 4. Upon any errors, the program should be terminated within the closure.
///    This can be done with *Zoe*'s [`Fail`] or [`OrFail`], or with
///    [`std::process::exit`].
/// 5. If the processing logic completes successfully, the closure should return
///    normally. [`GridTask::run_task`] will handle renaming the output files to
///    indicate successful completion. An exit code of 0 is returned.
///
/// [`Fail`]: zoe::data::err::Fail
pub struct GridTask {
    info: GridTaskInfo,
}

impl GridTask {
    /// Performs all the processing logic for the grid task as specified in
    /// closure `p`. The closure takes a single argument of type
    /// [`GridTaskInfo`], which can be used to select the proper inputs with
    /// [`GridTaskInfo::select_inputs`].
    ///
    /// See [`GridTask`] for more details.
    pub fn run_task<P>(self, p: P) -> !
    where
        P: FnOnce(&GridTaskInfo), {
        p(&self.info);

        let add_id = self.info.add_id;
        let add_id_tmp = self.info.add_id_tmp;

        for output_path in self.info.output_paths {
            let id = format!("{}", self.info.task_id);
            let tmp_path = add_id_tmp(&output_path, &id);
            let final_path = add_id(&output_path, &id);
            std::fs::rename(&tmp_path, &final_path).unwrap_or_die(&format!(
                "Failed to rename {tmp_path} to {final_path}",
                tmp_path = tmp_path.display(),
                final_path = final_path.display()
            ))
        }

        std::process::exit(0)
    }
}

/// Grid array-task environment information parsed from SGE.
#[derive(Clone, Debug)]
pub struct GridTaskInfo {
    pub task_id:       usize,
    pub task_first:    usize,
    pub task_last:     usize,
    pub task_stepsize: usize,
    pub scheduler:     GridScheduler,
    /// The output paths that will be formed via collation.
    output_paths:      Vec<PathBuf>,
    /// A pointer to [`GridCompatibleArgs::add_id`].
    add_id:            fn(&Path, &str) -> PathBuf,
    /// A pointer to [`GridCompatibleArgs::add_id_tmp`].
    add_id_tmp:        fn(&Path, &str) -> PathBuf,
}

impl GridTaskInfo {
    /// Downsamples the provided iterator so that only the inputs relevant to
    /// this task are included.
    ///
    /// Specifically, round-robin partitioning is used. For example:
    ///
    /// ```rust
    /// let grid_info = GridTaskInfo {
    ///     task_id:       2,
    ///     task_first:    1,
    ///     task_last:     10,
    ///     task_stepsize: 1,
    ///     scheduler:     GridScheduler::Sge,
    /// };
    ///
    /// let mut inputs_for_task = grid_info.select_inputs(0..50);
    ///
    /// assert_eq!(inputs_for_task.next(), Some(1));
    /// assert_eq!(inputs_for_task.next(), Some(11));
    /// assert_eq!(inputs_for_task.next(), Some(21));
    /// assert_eq!(inputs_for_task.next(), Some(31));
    /// assert_eq!(inputs_for_task.next(), Some(41));
    /// assert_eq!(inputs_for_task.next(), None);
    /// ```
    pub fn select_inputs<I>(&self, iter: I) -> impl Iterator<Item = I::Item>
    where
        I: Iterator, {
        // 0-based offset so we start on our partition
        let offset = (self.task_id - self.task_first) / self.task_stepsize;
        // modulus for interleaved partitioning
        let array_size = (self.task_last - self.task_first + 1).div_ceil(self.task_stepsize);

        iter.skip(offset).step_by(array_size)
    }
}

/// Supported grid engine schedulers.
#[non_exhaustive]
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub enum GridScheduler {
    Sge,
}

/// Any information pertaining to grid execution.
pub enum GridInfo {
    /// Information related to a request to use grid submission, which can be
    /// used to submit the job.
    Requested(GridRequestedInfo),
    /// Information about the task to which this execution belongs. This can be
    /// used to select the input records to handle.
    Task(GridTask),
}

impl GridInfo {
    /// Constructs a new [`GridInfo`] using `args`, `arg_matches`, and any
    /// environmental variables related to grid submission.
    ///
    /// ## Errors
    ///
    /// If `is_grid_task` is `true`, then parsing of the environmental variables
    /// must succeed. If `submit_grid_job` is `Some`, then the path to the
    /// current executable must be obtained without error, a scheduler must be
    /// located, and [`GridCompatibleArgs::log_path`] must succeed.
    fn new<T>(
        args: &mut T, arg_matches: &[ArgMatches], is_grid_task: bool, submit_grid_job: Option<usize>,
    ) -> std::io::Result<Option<Self>>
    where
        T: GridCompatibleArgs, {
        if is_grid_task {
            let grid_info = GridTaskInfo::new(args)?;

            // Update all the output paths to be temporary outputs for the task
            for path in args.outputs_mut() {
                *path = T::add_id_tmp(path, &format!("{}", grid_info.task_id));
            }

            Ok(Some(GridInfo::Task(GridTask { info: grid_info })))
        } else if let Some(task_count) = submit_grid_job {
            let mut grid_command = vec![
                current_exe()
                    .with_context("Failed to get the path to the current executable")?
                    .into_os_string(),
            ];

            let command = T::Cli::command();

            for arg in command.get_positionals().chain(command.get_opts()) {
                let id = arg.get_id().as_str();

                if id == T::Cli::SUBMIT_GRID_JOB_ID {
                    continue;
                } else if id == T::Cli::IS_GRID_TASK_ID {
                    push_arg_name(&mut grid_command, arg);
                }

                for arg_matches in arg_matches {
                    if arg_matches.try_contains_id(id).is_err()
                        || !matches!(arg_matches.value_source(id), Some(ValueSource::CommandLine))
                    {
                        continue;
                    }

                    match arg.get_action() {
                        ArgAction::Set | ArgAction::Append => {
                            if arg.is_positional() {
                                if let Some(values) = arg_matches.get_raw(id) {
                                    grid_command.extend(values.map(OsStr::to_os_string));
                                }
                            } else if let Some(values) = arg_matches.get_raw(id) {
                                for value in values {
                                    push_arg_name(&mut grid_command, arg);
                                    grid_command.push(value.to_os_string());
                                }
                            }
                        }
                        ArgAction::SetTrue | ArgAction::SetFalse => {
                            push_arg_name(&mut grid_command, arg);
                        }
                        ArgAction::Count => {
                            let n = arg_matches.get_count(id);
                            for _ in 0..n {
                                push_arg_name(&mut grid_command, arg);
                            }
                        }
                        _ => {}
                    }
                }
            }

            let scheduler = pick_submission_scheduler()?;

            Ok(Some(GridInfo::Requested(GridRequestedInfo {
                command: grid_command,
                task_count,
                scheduler,
                output_paths: args.outputs_mut().map(|x| (*x).clone()).collect(),
                log_path: args.log_path()?,
                add_id: T::add_id,
                add_id_tmp: T::add_id_tmp,
            })))
        } else {
            Ok(None)
        }
    }
}

/// A wrapper around a [`Parser`] that also stores cloned copies of the
/// [`ArgMatches`] used in parsing, for use in grid submission.
struct GridParser<P> {
    /// The inner parser being wrapped.
    parser:      P,
    /// Any [`ArgMatches`] passed via [`from_arg_matches`] or
    /// [`update_from_arg_matches`].
    ///
    /// [`from_arg_matches`]: FromArgMatches::from_arg_matches
    /// [`update_from_arg_matches`]: FromArgMatches::update_from_arg_matches
    arg_matches: Vec<ArgMatches>,
}

impl<P: GridCompatibleCli> FromArgMatches for GridParser<P> {
    fn from_arg_matches(matches: &ArgMatches) -> Result<Self, clap::Error> {
        <P as FromArgMatches>::from_arg_matches(matches).map(|parser| Self {
            parser,
            arg_matches: vec![matches.clone()],
        })
    }

    fn update_from_arg_matches(&mut self, matches: &ArgMatches) -> Result<(), clap::Error> {
        self.parser.update_from_arg_matches(matches)?;
        self.arg_matches.push(matches.clone());
        Ok(())
    }
}

impl<P: GridCompatibleCli> CommandFactory for GridParser<P> {
    #[inline]
    fn command() -> clap::Command {
        P::command()
    }

    #[inline]
    fn command_for_update() -> clap::Command {
        P::command_for_update()
    }
}

impl<P: GridCompatibleCli> Parser for GridParser<P> {}

/// A helper function for [GridInfo::new] which appends the specified keyword
/// argument to the command.
fn push_arg_name(grid_command: &mut Vec<OsString>, arg: &Arg) {
    if let Some(long) = arg.get_long() {
        grid_command.push(format!("--{long}").into());
    } else if let Some(short) = arg.get_short() {
        grid_command.push(format!("-{short}").into());
    }
}
