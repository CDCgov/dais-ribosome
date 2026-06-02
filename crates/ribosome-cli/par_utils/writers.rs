use std::{
    error::Error,
    fmt::Display,
    io::Write,
    sync::{
        Arc, OnceLock,
        mpsc::{RecvError, SendError},
    },
};
use zoe::data::err::GetCode;

/// A clonable writer supporting writing from multiple threads via an [`mpsc`]
/// channel.
///
/// A single dedicated thread is used for writing to the underlying writer to
/// avoid interleaved writes. Calling [`write`] is non-blocking, instead
/// queueing the bytes in the channel. Calling [`flush`] causes the thread to
/// block until all previous data in the channel has been written to the
/// underlying writer.
///
/// To ensure that errors are not silently ignored, for each linked
/// [`WriterThreaded`], ensure that either [`flush`] is called on it, or it is
/// dropped before a [`flush`] call on a different clone of the
/// [`WriterThreaded`].
///
/// [`flush`]: WriterThreaded::flush
/// [`write`]: WriterThreaded::write
/// [`mpsc`]: std::sync::mpsc
#[derive(Clone, Debug)]
pub struct WriterThreaded {
    /// The sending portion of the channel.
    sender: std::sync::mpsc::Sender<Msg>,

    /// The first IO error produced by the underlying writer in a
    /// commonly-accessible location. The outer [`Arc`] ensures that the memory
    /// is shared between all the threads. The [`OnceLock`] acts as an `Option`
    /// that is written to once. The [`SharedIoError`] contains another [`Arc`]
    /// to allow the inner error to be clone-able.
    writer_error: Arc<OnceLock<SharedIoError>>,
}

enum Msg {
    /// A line of data to be written to the underlying writer.
    Data(Vec<u8>),
    /// A request to flush the underlying writer and all previously received
    /// messages in the channel. The sender is used to communicate when the
    /// flush has been successfully completed.
    Flush(std::sync::mpsc::Sender<()>),
}

impl Write for WriterThreaded {
    /// Writes a buffer into the writer, returning how many bytes were queued
    /// for writing.
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self.sender.send(Msg::Data(buf.to_vec())) {
            Ok(()) => Ok(buf.len()),
            Err(SendError(_)) => Err(self.writer_error.get().map_or_else(
                || std::io::Error::new(std::io::ErrorKind::BrokenPipe, "writer thread is no longer running"),
                Into::into,
            )),
        }
    }

    /// Flushes the contents of the queue before this call was made. This
    /// function will block until the flush is successful.
    ///
    /// To ensure that errors are not silently ignored, for each linked
    /// [`WriterThreaded`], ensure that either [`flush`] is called on it, or it
    /// is dropped before a [`flush`] call on a different clone of the
    /// [`WriterThreaded`]. This is similar to a [`BufWriter`], where failing to
    /// call [`flush`] may result in errors getting silently ignored.
    ///
    /// [`BufWriter`]: std::io::BufWriter
    /// [`flush`]: WriterThreaded::flush
    fn flush(&mut self) -> std::io::Result<()> {
        let (sender, receiver) = std::sync::mpsc::channel();

        match self.sender.send(Msg::Flush(sender)) {
            Ok(()) => match receiver.recv() {
                Ok(()) => Ok(()),
                Err(RecvError) => Err(self.writer_error.get().map_or_else(
                    || std::io::Error::new(std::io::ErrorKind::BrokenPipe, "writer thread is no longer running"),
                    Into::into,
                )),
            },
            Err(SendError(_)) => Err(self.writer_error.get().map_or_else(
                || std::io::Error::new(std::io::ErrorKind::BrokenPipe, "writer thread is no longer running"),
                Into::into,
            )),
        }
    }

    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        self.write(buf)?;
        Ok(())
    }

    fn write_fmt(&mut self, args: std::fmt::Arguments<'_>) -> std::io::Result<()> {
        let bytes = if let Some(s) = args.as_str() {
            s.as_bytes().to_vec()
        } else {
            std::fmt::format(args).into_bytes()
        };

        self.write_all(&bytes)
    }
}

impl WriterThreaded {
    /// Constructs a [`WriterThreaded`] from a regular writer by moving it into
    /// a thread and creating a channel.
    #[inline]
    #[must_use]
    pub fn new<W>(mut writer: W) -> Self
    where
        W: Write + Send + 'static, {
        let (sender, receiver) = std::sync::mpsc::channel();

        let writer_error = Arc::new(OnceLock::new());

        let thread_writer_error = writer_error.clone();

        // We do not bind the thread handle, but the thread will still run to
        // completion (it will run until an error, or until all senders are
        // dropped)
        std::thread::spawn(move || {
            while let Ok(msg) = receiver.recv() {
                match msg {
                    Msg::Data(bytes) => {
                        let res = writer.write_all(&bytes);
                        if let Err(e) = res {
                            thread_writer_error.get_or_init(|| e.into());
                            return;
                        }
                    }
                    Msg::Flush(sender) => {
                        let res = writer.flush();
                        if let Err(e) = res {
                            thread_writer_error.get_or_init(|| e.into());
                            return;
                        }

                        // Signal to the thread that called flush that the flush
                        // has now been completed. This should never fail since
                        // the thread calls `recv` immediately after sending the
                        // Msg, but if it does for some pathological case, we
                        // may as well poison as much as possible.
                        let res = sender.send(());
                        if let Err(SendError(())) = res {
                            thread_writer_error
                                .get_or_init(|| std::io::Error::other("failed to confirm successful flush").into());
                            return;
                        }
                    }
                }
            }

            // Error is not propagated, similarly to how not calling flush on
            // BufWriter can cause silent errors when Drop is called.
            let _ = writer.flush();
        });

        Self { sender, writer_error }
    }
}

impl Drop for WriterThreaded {
    /// Drops the [`WriterThreaded`], blocking until [`flush`] finishes.
    ///
    /// [`flush`]: WriterThreaded::flush
    fn drop(&mut self) {
        // If it succeeds, then all contents from this writer have been flushed.
        // If it fails, then either an IO error occured and has been populated
        // in writer_error (so that a future explicit flush call will yield it),
        // or the writer thread unexpectedly ended without an IO error (in which
        // case a future explicit flush call will also yield this error, since
        // the thread cannot be restarted)
        let _ = self.flush();
    }
}

/// A [`std::io::Error`] for which [`Clone`] can be called (allowing many
/// threads to yield the same underlying error while preserving any source
/// errors).
#[derive(Clone, Debug)]
struct SharedIoError(Arc<std::io::Error>);

impl Display for SharedIoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl Error for SharedIoError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        self.0.source()
    }
}

impl GetCode for SharedIoError {
    fn get_code(&self) -> i32 {
        self.0.get_code()
    }
}

impl From<&SharedIoError> for std::io::Error {
    fn from(value: &SharedIoError) -> Self {
        std::io::Error::new(value.0.kind(), value.clone())
    }
}

impl From<SharedIoError> for std::io::Error {
    fn from(value: SharedIoError) -> Self {
        std::io::Error::new(value.0.kind(), value)
    }
}

impl From<std::io::Error> for SharedIoError {
    fn from(value: std::io::Error) -> Self {
        Self(Arc::new(value))
    }
}
