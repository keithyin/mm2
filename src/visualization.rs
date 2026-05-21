use minimap2::Strand;

fn revcomp(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            'N' | 'n' => 'N',
            _ => 'N',
        })
        .collect()
}

fn parse_cigar(cigar: &str) -> Vec<(usize, char)> {
    let mut result = Vec::new();
    let mut num = String::new();

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num.push(c);
        } else {
            let len = num.parse::<usize>().unwrap();
            result.push((len, c));
            num.clear();
        }
    }

    result
}

pub fn blast_like_alignment(
    query: &str,
    reference: &str,
    aln: &minimap2::Mapping,
    line_width: Option<usize>,
) {
    let line_width = line_width.unwrap_or(110);
    let mut q_st = aln.query_start as usize;
    let mut q_en = aln.query_end as usize;
    let r_st = aln.target_start as usize;
    let mut strand = aln.strand;
    let cigar_str = aln
        .alignment
        .as_ref()
        .unwrap()
        .cigar
        .as_deref()
        .unwrap();

    let mut query_seq = query.to_string();

    // 负链处理
    if strand == Strand::Reverse {
        query_seq = revcomp(&query_seq);

        let new_q_st = query_seq.len() - q_en;
        let new_q_en = query_seq.len() - q_st;

        q_st = new_q_st;
        q_en = new_q_en;
        strand = Strand::Forward;
    }

    let query_chars: Vec<char> = query_seq.chars().collect();
    let ref_chars: Vec<char> = reference.chars().collect();

    let mut q_pos = q_st;
    let mut r_pos = r_st;

    let mut q_line: Vec<char> = Vec::new();
    let mut r_line: Vec<char> = Vec::new();
    let mut mid_line: Vec<char> = Vec::new();

    let mut q_coords: Vec<Option<usize>> = Vec::new();
    let mut r_coords: Vec<Option<usize>> = Vec::new();

    for (length, op) in cigar_str.to_vec() {
        match op {
            0 | 7 | 8 => {
                for _ in 0..length {
                    let qb = query_chars[q_pos];
                    let rb = ref_chars[r_pos];

                    q_line.push(qb);
                    r_line.push(rb);

                    if qb == rb {
                        mid_line.push('|');
                    } else {
                        mid_line.push(' ');
                    }

                    q_coords.push(Some(q_pos));
                    r_coords.push(Some(r_pos));

                    q_pos += 1;
                    r_pos += 1;
                }
            }

            1 => {
                for _ in 0..length {
                    let qb = query_chars[q_pos];

                    q_line.push(qb);
                    r_line.push('-');
                    mid_line.push(' ');

                    q_coords.push(Some(q_pos));
                    r_coords.push(None);

                    q_pos += 1;
                }
            }

            2 => {
                for _ in 0..length {
                    let rb = ref_chars[r_pos];

                    q_line.push('-');
                    r_line.push(rb);
                    mid_line.push(' ');

                    q_coords.push(None);
                    r_coords.push(Some(r_pos));

                    r_pos += 1;
                }
            }

            4 => {}

            // H / P 等忽略
            _ => {}
        }
    }

    let total_len = q_line.len();

    for i in (0..total_len).step_by(line_width) {
        let end = std::cmp::min(i + line_width, total_len);

        let q_seg = &q_line[i..end];
        let r_seg = &r_line[i..end];
        let m_seg = &mid_line[i..end];

        let qc_seg = &q_coords[i..end];
        let rc_seg = &r_coords[i..end];

        let mut q_start: Option<usize> = None;
        let mut q_end: Option<usize> = None;

        for coord in qc_seg {
            if let Some(c) = coord {
                if q_start.is_none() {
                    q_start = Some(*c);
                }
                q_end = Some(*c);
            }
        }

        if q_start.is_none() {
            continue;
        }

        let mut q_start = q_start.unwrap() + 1;
        let mut q_end = q_end.unwrap() + 1;

        let mut r_start: Option<usize> = None;
        let mut r_end: Option<usize> = None;

        for coord in rc_seg {
            if let Some(c) = coord {
                if r_start.is_none() {
                    r_start = Some(*c);
                }
                r_end = Some(*c);
            }
        }

        if r_start.is_none() {
            continue;
        }

        let r_start = r_start.unwrap() + 1;
        let r_end = r_end.unwrap() + 1;

        if aln.strand == Strand::Reverse {
            let tmp = query_seq.len() - q_end;
            q_end = query_seq.len() - q_start;
            q_start = tmp;

            q_start += 1;
            q_end += 1;

            std::mem::swap(&mut q_start, &mut q_end);
        }

        let r_str: String = r_seg.iter().collect();
        let q_str: String = q_seg.iter().collect();
        let m_str: String = m_seg.iter().collect();

        println!("Sbjct  {:<6} {} {}", r_start, r_str, r_end);
        println!("              {}", m_str);
        println!("Query  {:<6} {} {}", q_start, q_str, q_end);
        println!();
    }
}
