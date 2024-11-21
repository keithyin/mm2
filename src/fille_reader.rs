use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

fn fastx_header_line_to_header(header_line: &str) -> String {
    let header = if header_line[1..].contains(" ") {
        header_line[1..].split_once(" ").unwrap().0.to_string()
    } else {
        header_line[1..].to_string()
    };
    return header;
}

/// 表示一个 FASTA 记录的结构体
#[derive(Debug, PartialEq)]
pub struct QueryRecord {
    pub qname: String,
    pub sequence: String,

    pub ch: Option<usize>,
    pub np: Option<usize>,
}

impl QueryRecord {
    pub fn from_bam_record(record: &rust_htslib::bam::Record, qname_suffix: Option<&str>) -> Self {
        let mut qname = unsafe { String::from_utf8_unchecked(record.qname().to_vec()) };
        if let Some(suffix) = qname_suffix {
            qname.push_str(suffix);
        }

        let record_ext = gskits::gsbam::bam_record_ext::BamRecordExt::new(record);
        let seq = unsafe { String::from_utf8_unchecked(record.seq().as_bytes()) };

        Self {
            qname: qname,
            sequence: seq,
            ch: record_ext.get_ch(),
            np: record_ext.get_np(),
        }
    }
}

/// 解析 FASTA 文件，返回一个包含所有记录的 Vec
pub fn read_fasta<P: AsRef<Path>>(file_path: P) -> io::Result<Vec<QueryRecord>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    let mut records = Vec::new();
    let mut current_header = None;
    let mut current_sequence = String::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if let Some(header) = current_header.take() {
                records.push(QueryRecord {
                    qname: header,
                    sequence: current_sequence.clone(),
                    ch: None,
                    np: None,
                });
                current_sequence.clear();
            }
            current_header = Some(fastx_header_line_to_header(&line));
        } else {
            current_sequence.push_str(&line);
        }
    }

    // 添加最后一条记录
    if let Some(header) = current_header {
        records.push(QueryRecord {
            qname: header,
            sequence: current_sequence,
            ch: None,
            np: None,
        });
    }

    Ok(records)
}

pub struct FastaFileReader {
    fname: String,
    lines: Option<io::Lines<BufReader<File>>>,
    current_header: Option<String>,
}

impl FastaFileReader {
    pub fn new(fname: String) -> Self {
        Self {
            fname: fname,
            lines: None,
            current_header: None,
        }
    }
}

impl Iterator for FastaFileReader {
    type Item = QueryRecord;
    fn next(&mut self) -> Option<Self::Item> {
        if self.lines.is_none() {
            let file = File::open(&self.fname).expect(&format!("open file error: {}", self.fname));
            let lines = BufReader::new(file).lines();
            self.lines = Some(lines);
            let first_header = self.lines.as_mut().unwrap().next().unwrap().unwrap();
            assert!(first_header.starts_with(">"));

            self.current_header = Some(fastx_header_line_to_header(&first_header));
        }

        if self.current_header.is_none() {
            return None;
        }

        let mut fasta_record = QueryRecord {
            qname: self.current_header.take().unwrap(),
            sequence: String::new(),
            ch: None,
            np: None,
        };
        while let Some(line) = self.lines.as_mut().unwrap().next() {
            let line = line.unwrap();
            if line.starts_with(">") {
                self.current_header = Some(fastx_header_line_to_header(&line));
                return Some(fasta_record);
            }

            fasta_record.sequence.push_str(line.trim());
        }

        return Some(fasta_record);
    }
}

// header, seq, qual
#[derive(Debug)]
pub struct FastqRecord(pub String, pub String, pub String);

impl FastqRecord {
    pub fn eqaul(&self, other: &FastqRecord) -> bool {
        self.1.eq(&other.1) && self.2.eq(&other.2)
    }
}

pub struct FastqReaderIter<'a> {
    reader: &'a mut dyn BufRead,
}

impl<'a> FastqReaderIter<'a> {
    pub fn new(reader: &'a mut dyn BufRead) -> Self {
        Self { reader }
    }

    fn read_one_line(&mut self) -> Option<String> {
        let mut line = String::new();
        if let Ok(n) = self.reader.read_line(&mut line) {
            if n == 0 {
                return None;
            }
            line = line.trim().to_string();
        } else {
            return None;
        }
        return Some(line);
    }
}

impl<'a> Iterator for FastqReaderIter<'a> {
    type Item = FastqRecord;
    fn next(&mut self) -> Option<Self::Item> {
        let header = self.read_one_line();
        if header.is_none() || header.as_ref().unwrap().trim().len() == 0 {
            return None;
        }
        let mut header = header.unwrap();
        if !header.starts_with("@") {
            panic!("header:'{}' not a valid fastq header", header);
        }

        header = fastx_header_line_to_header(&header);

        let seq = self
            .read_one_line()
            .expect("not a valid FastqRecord")
            .trim()
            .to_string();
        let plus = self
            .read_one_line()
            .expect("not a valid FastqRecord")
            .trim()
            .to_string();

        if !plus.starts_with("+") {
            panic!("plus:'{}' not a valid fastq plus", plus);
        }

        let qual = self
            .read_one_line()
            .expect("not a valid FastqRecord")
            .trim()
            .to_string();


        Some(FastqRecord(header, seq, qual))
    }
}

#[cfg(test)]
mod test {
    use super::{read_fasta, FastaFileReader};

    #[test]
    fn test_read_fasta() {
        let fname = "test_data/seq.fasta";
        let records = read_fasta(fname).unwrap();
        assert_eq!(records[0].sequence, "ATCGTACGTACGTACGTAGC");
        assert_eq!(records[1].sequence, "GGGTTTAAACCCGGGTTT");
    }

    #[test]
    fn test_fasta_file_reader() {
        let fname = "test_data/seq.fasta";
        let mut reader = FastaFileReader::new(fname.to_string());

        assert_eq!(reader.next().unwrap().sequence, "ATCGTACGTACGTACGTAGC");
        assert_eq!(reader.next().unwrap().sequence, "GGGTTTAAACCCGGGTTT");
    }
}
