use crate::par_utils::writers::WriterThreaded;
use dais_ribosome::tsv::{
    DelRowView, DelWriter, DeletedSeqRowView, EmptySeqRowView, Finish, GenDelRowView, GenDelWriter, GenInsRowView,
    GenInsWriter, GenSeqRowView, GenSeqWriter, InsRowView, InsWriter, SeqDataView, SeqWriter,
};
use std::io::Write;

impl Finish for WriterThreaded {
    fn finish(mut self) -> std::io::Result<()> {
        self.flush()
    }
}

impl SeqWriter for WriterThreaded {
    fn write_seq_data(&mut self, data: &SeqDataView) -> std::io::Result<()> {
        write!(self, "{data}")
    }

    fn write_empty_seq_row(&mut self, row: &EmptySeqRowView<'_>) -> std::io::Result<()> {
        write!(self, "{row}")
    }

    fn write_deleted_seq_row(&mut self, row: &DeletedSeqRowView<'_>) -> std::io::Result<()> {
        write!(self, "{row}")
    }
}

impl InsWriter for WriterThreaded {
    fn write_ins_row(&mut self, row: &InsRowView) -> std::io::Result<()> {
        write!(self, "{row}")
    }
}

impl DelWriter for WriterThreaded {
    fn write_del_row(&mut self, row: &DelRowView) -> std::io::Result<()> {
        write!(self, "{row}")
    }
}

impl GenSeqWriter for WriterThreaded {
    fn write_gen_seq_row(&mut self, row: &GenSeqRowView) -> std::io::Result<()> {
        write!(self, "{row}")
    }
}

impl GenInsWriter for WriterThreaded {
    fn write_gen_ins_row(&mut self, row: &GenInsRowView) -> std::io::Result<()> {
        write!(self, "{row}")
    }
}

impl GenDelWriter for WriterThreaded {
    fn write_gen_del_row(&mut self, row: &GenDelRowView) -> std::io::Result<()> {
        write!(self, "{row}")
    }
}
