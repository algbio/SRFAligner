extern crate daachorse;
use daachorse::DoubleArrayAhoCorasick;
use std::env;

use std::{
    fs::File,
    io::{prelude::*, BufReader},
    path::Path,
};

fn lines_from_file(filename: impl AsRef<Path>) -> Vec<String> {
    let file = File::open(filename).expect("no such file");
    let buf = BufReader::new(file);
    buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect()
}

// ---

fn main() {
    let args: Vec<String> = env::args().collect(); // nodes.txt ids.txt reads.txt read_ids.txt

    let nodes : Vec<String> = lines_from_file(args[1].clone());
    let nodeslen : Vec<usize> = nodes.iter().cloned().map(|n| n.len()).collect();
    let ids : Vec<String> = lines_from_file(args[2].clone());
    let reads : Vec<String> = lines_from_file(args[3].clone());
    let rids : Vec<String> = lines_from_file(args[4].clone());
    let automaton : DoubleArrayAhoCorasick<usize> = DoubleArrayAhoCorasick::new(nodes).unwrap();

    let mut printed = false;
    for (read,readid) in reads.into_iter().zip(rids.into_iter()) {
        let readlen = read.len();
        let it = automaton.find_overlapping_iter(read);

        for arg in it {
            if !printed {
                printed = true;
            } else {
                print!("\n");
            }
            print!("{}\t{}\t{}\t{}\t+\t>{}\t{}\t0\t{}\t0\t0\t255", readid, readlen, arg.start(), arg.end(), ids[arg.value()], nodeslen[arg.value()],  nodeslen[arg.value()]);
        }
    }
}
