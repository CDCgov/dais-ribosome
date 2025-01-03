#![feature(let_chains)]

use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug, PartialEq, Eq, Hash)]
struct CodonKey {
    position: u32,
    codon: [u8; 3],
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Please provide a filename as argument");
        return;
    }

    let file = File::open(&args[1]).unwrap_or_else(|e| {
        eprintln!("Error opening file: {}", e);
        std::process::exit(1);
    });

    let reader = BufReader::new(file);
    let mut vec_of_hashmaps = Vec::new();
    let mut ref_product = Vec::new();
    let mut index = 0;
    let mut vec_sites = Vec::new();

    for line in reader.lines().map_while(Result::ok) {
        if line.contains('|') {
            ref_product.push(line.clone());
            vec_of_hashmaps.push(HashMap::new());
            vec_sites.push(HashSet::new());
            index = vec_of_hashmaps.len() - 1
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() == 3
            && let (Ok(position), Ok(codon), Ok(count)) = (
                parts[0].parse(),
                parts[1].as_bytes()[..3].try_into(),
                parts[2].parse::<u32>(),
            )
        {
            vec_sites[index].insert(position);
            vec_of_hashmaps[index].insert(CodonKey { position, codon }, count);
        }
    }

    let mut sum_codons = 0;
    let mut sum_sites = 0;
    let mut sum_counts = 0;

    for (set, set_name) in ref_product.iter().enumerate() {
        let codons = vec_of_hashmaps[set].len();

        let mut split = set_name.split('|');
        let (ref_id, protein) = (split.next().unwrap(), split.next().unwrap());
        let sites = vec_sites[set].len();
        sum_sites += sites;
        sum_codons += codons;

        for count in vec_of_hashmaps[set].values() {
            sum_counts += count;
        }

        println!("{ref_id}\t{protein}\t{codons}\t{sites}");
    }
    println!("Total\tALL\t{sum_codons}\t{sum_sites}");
    println!("Total counts:\t{sum_counts}");
}
