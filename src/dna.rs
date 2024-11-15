static COMPLEMENT_TABLE: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'C' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table
};

pub fn reverse_complement(dna: &str) -> String {
    // 预分配结果字符串容量
    let mut result = dna.as_bytes().to_vec();
    let len = result.len();

    for (i, &base) in dna.as_bytes().iter().enumerate() {
        let complement = COMPLEMENT_TABLE[base as usize];
        
        result[len - i - 1] = complement;
    }

    // 将结果转回字符串
    unsafe { String::from_utf8_unchecked(result) }
}
