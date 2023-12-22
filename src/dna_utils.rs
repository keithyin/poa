
use std::collections::HashMap;

use lazy_static::lazy_static;

lazy_static! {
    static ref BASE_COMPLEMENT: HashMap<u8, u8> = {
        let mut map = HashMap::new();
        map.insert('A' as u8, 'T' as u8);
        map.insert('T' as u8, 'A' as u8);

        map.insert('C' as u8, 'G' as u8);
        map.insert('G' as u8, 'C' as u8);
        map
    };
}

pub fn reverse_complement(read: &[u8]) -> Vec<u8> {
    read.iter()
        .map(|base| *BASE_COMPLEMENT.get(base).unwrap())
        .collect::<_>()
}