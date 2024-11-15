use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use rust_htslib::bam::Read;

/// 表示一个 FASTA 记录的结构体
#[derive(Debug, PartialEq)]
pub struct QueryRecord {
    pub qname: String,
    pub sequence: String,
}

impl QueryRecord {
    pub fn from_bam_record(record: &rust_htslib::bam::Record) -> Self {
        let qname = unsafe {
            String::from_utf8_unchecked(record.qname().to_vec())
        };

        let seq = unsafe {
            String::from_utf8_unchecked(record.seq().as_bytes())
        };

        Self { qname: qname, sequence: seq }

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
                });
                current_sequence.clear();
            }
            current_header = Some(line[1..].to_string()); // 去掉 '>'
        } else {
            current_sequence.push_str(&line);
        }
    }

    // 添加最后一条记录
    if let Some(header) = current_header {
        records.push(QueryRecord {
            qname: header,
            sequence: current_sequence,
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
            self.current_header = Some(first_header[1..].to_string());
        }

        if self.current_header.is_none() {
            return None;
        }

        let mut fasta_record = QueryRecord {
            qname: self.current_header.take().unwrap(),
            sequence: String::new(),
        };
        while let Some(line) = self.lines.as_mut().unwrap().next() {
            let line = line.unwrap();
            if line.starts_with(">") {
                self.current_header = Some(line[1..].to_string());
                return Some(fasta_record);
            }

            fasta_record.sequence.push_str(line.trim());
        }

        return Some(fasta_record);
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