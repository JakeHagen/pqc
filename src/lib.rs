use std::str;
use ndarray::OwnedRepr;
use phf::phf_map;
use serde::{Serialize, Deserialize};

static TI: phf::Map<&[u8], phf::Map<&[u8], bool>> = phf_map! {
    b"A" => phf_map! {
        b"G" => true,
        b"T" => false,
        b"C" => false,
    },
    b"C" => phf_map! {
        b"T" => true,
        b"A" => false,
        b"G" => false,
    },
    b"G" => phf_map! {
        b"A" => true,
        b"C" => false,
        b"T" => false,
    },
    b"T" => phf_map! {
        b"C" => true,
        b"A" => false,
        b"G" => false,
    },
};

fn decode_genotype(g: i32) -> i32 {
    return (g >> 1) - 1
}

fn is_ti(r: &[u8], a: &[u8]) -> bool {
    return *TI.get(r).unwrap().get(a).unwrap()
}

pub fn get_popfreq(rec: &rust_htslib::bcf::Record, field: &str) -> f32 {
    match rec.info(field.as_bytes()).float() {
        Ok(pf) => match pf {
            Some(pf) => return pf[0],
            None => return 0.0 
        },
        Err(_) => panic!("something went wrong when trying to get the population frequency"),
    }
}

fn is_snv(gt: i32, alleles: &Vec<&[u8]>) -> bool {
    return alleles[0].len() == alleles[gt as usize].len() && alleles[gt as usize].len() == 1
}

fn is_hom(gt: &[i32]) -> bool {
    return gt[0] == gt[1]
}

pub fn _is_variant(gt: &[i32]) -> bool {
    if gt[0] + gt[1] < 4 {
        return false
    }
    if gt[0] + gt[1] > 4 {
        return true
    }
    for g in gt {
        if *g == 4 {
            return true
        }
    }

    return false
}
pub fn is_variant(gt: &[i32]) -> (bool, usize) {
    for (vidx, g) in gt.iter().enumerate() {
        if *g >= 4 {
            return (true, vidx)
        }
    }
    return (false, 0)
}

const NVARS: usize = 0;
const NSNVS: usize = 1;
const NINDELS: usize = 2;
const NHETS: usize = 3;
const NHOMALTS: usize = 4;
const NTIS: usize = 5;
const NTVS: usize = 6;

pub const RAW: usize = 0;
pub const QUAL: usize = 1;
pub const RARE: usize = 2;
pub const QUAL_RARE: usize = 3;

/*
#[derive(Serialize, Deserialize)]
pub struct Metrics {
    gq: i32,
    freq: f32,
    population: String,
    samples: Vec<String>,
    raw: Vec<[i32;7]>,
    qual: Vec<[i32;7]>,
    rare: Vec<[i32;7]>,
    qual_rare: Vec<[i32;7]>,
}

impl Metrics {
    pub fn new(samples: Vec<String>, gq: i32, freq: f32, population: String) -> Metrics {
        let ns = samples.len();
        return Metrics {
            gq,
            freq,
            population,
            samples,
            raw: vec![[0,0,0,0,0,0,0]; ns],
            qual: vec![[0,0,0,0,0,0,0]; ns],
            rare: vec![[0,0,0,0,0,0,0]; ns],
            qual_rare: vec![[0,0,0,0,0,0,0]; ns],
        }
    }
    pub fn add_raw(&mut self, sidx: usize, arr: &[i32;7]) {
        for (i, n) in arr.iter().enumerate() {
            self.raw[sidx][i] += n;
        }
    }
    pub fn add_qual(&mut self, sidx: usize, arr: &[i32;7]) {
        for (i, n) in arr.iter().enumerate() {
            self.qual[sidx][i] += n;
        }
    }
    pub fn add_rare(&mut self, sidx: usize, arr: &[i32;7]) {
        for (i, n) in arr.iter().enumerate() {
            self.rare[sidx][i] += n;
        }
    }
    pub fn add_qual_rare(&mut self, sidx: usize, arr: &[i32;7]) {
        for (i, n) in arr.iter().enumerate() {
            self.qual_rare[sidx][i] += n;
        }
    }
    pub fn cat(&mut self, metrics: Metrics) {
        if self.gq != metrics.gq {
            panic!("GQ thresholds do not match");
        }
        if self.freq != metrics.freq {
            panic!("Rarity thresholds do not match");
        }
        if self.population != metrics.population {
            panic!("The population AF source does not match");
        }
        if self.samples != metrics.samples {
            panic!("Samples do not match");
        }
        for (sidx, arr) in metrics.raw.iter().enumerate() {
            for (i, a) in arr.iter().enumerate() {
                self.raw[sidx][i] += a
            }
        }
        for (sidx, arr) in metrics.qual.iter().enumerate() {
            for (i, a) in arr.iter().enumerate() {
                self.qual[sidx][i] += a
            }
        }
        for (sidx, arr) in metrics.rare.iter().enumerate() {
            for (i, a) in arr.iter().enumerate() {
                self.rare[sidx][i] += a
            }
        }
        for (sidx, arr) in metrics.qual_rare.iter().enumerate() {
            for (i, a) in arr.iter().enumerate() {
                self.qual_rare[sidx][i] += a
            }
        }
    }
    pub fn to_csv(self, iowtr: Box<dyn std::io::Write>) {
        let mut wtr = csv::WriterBuilder::new().has_headers(false).from_writer(iowtr);
        _ = wtr.write_record(&[
            "sample", "n_vars", "n_snvs", "n_indels", "n_hets", "n_hom_alts", "n_tis", "n_tvs",
            &format!("gq{}_n_vars", self.gq), &format!("gq{}_n_snvs", self.gq), &format!("gq{}_n_indels", self.gq), &format!("gq{}_n_hets", self.gq), &format!("gq{}_n_hom_alts", self.gq), &format!("gq{}_n_tis", self.gq), &format!("gq{}_n_tvs", self.gq),
            &format!("{}{}_n_vars", self.population, self.freq), &format!("{}{}_n_snvs", self.population, self.freq), &format!("{}{}_n_indels", self.population, self.freq), &format!("{}{}_n_hets", self.population, self.freq), &format!("{}{}_n_hom_alts", self.population, self.freq), &format!("{}{}_n_tis", self.population, self.freq), &format!("{}{}_n_tvs", self.population, self.freq),
            &format!("gq{}_{}{}_n_vars", self.gq, self.population, self.freq), &format!("gq{}_{}{}_n_snvs", self.gq, self.population, self.freq), &format!("gq{}_{}{}_n_indels", self.gq, self.population, self.freq), &format!("gq{}_{}{}_n_hets", self.gq, self.population, self.freq), &format!("gq{}_{}{}_n_hom_alts", self.gq, self.population, self.freq), &format!("gq{}_{}{}_n_tis", self.gq, self.population, self.freq), &format!("gq{}_{}{}_n_tvs", self.gq, self.population, self.freq),
        ]);
        for sidx in 0..self.samples.len() {
            _ = wtr.write_record(&[
                &self.samples[sidx], &format!("{}", self.raw[sidx][NVARS]), &format!("{}", self.raw[sidx][NSNVS]), &format!("{}", self.raw[sidx][NINDELS]), &format!("{}", self.raw[sidx][NHETS]), &format!("{}", self.raw[sidx][NHOMALTS]), &format!("{}", self.raw[sidx][NTIS]), &format!("{}", self.raw[sidx][NTVS]),
                &format!("{}", self.qual[sidx][NVARS]), &format!("{}", self.qual[sidx][NSNVS]), &format!("{}", self.qual[sidx][NINDELS]), &format!("{}", self.qual[sidx][NHETS]), &format!("{}", self.qual[sidx][NHOMALTS]), &format!("{}", self.qual[sidx][NTIS]), &format!("{}", self.qual[sidx][NTVS]),
                &format!("{}", self.rare[sidx][NVARS]), &format!("{}", self.rare[sidx][NSNVS]), &format!("{}", self.rare[sidx][NINDELS]), &format!("{}", self.rare[sidx][NHETS]), &format!("{}", self.rare[sidx][NHOMALTS]), &format!("{}", self.rare[sidx][NTIS]), &format!("{}", self.rare[sidx][NTVS]),
                &format!("{}", self.qual_rare[sidx][NVARS]), &format!("{}", self.qual_rare[sidx][NSNVS]), &format!("{}", self.qual_rare[sidx][NINDELS]), &format!("{}", self.qual_rare[sidx][NHETS]), &format!("{}", self.qual_rare[sidx][NHOMALTS]), &format!("{}", self.qual_rare[sidx][NTIS]), &format!("{}", self.qual_rare[sidx][NTVS]),
            ]);
        }
    }
}
*/
#[derive(Serialize, Deserialize)]
pub struct Metrics {
    gq: i32,
    freq: f32,
    population: String,
    samples: Vec<String>,
    data: ndarray::ArrayBase<OwnedRepr<i32>, ndarray::Dim<[usize; 3]>>,
}

impl Metrics {
    pub fn new(samples: Vec<String>, gq: i32, freq: f32, population: String) -> Metrics {
        let ns = samples.len();
        return Metrics {
            gq,
            freq,
            population,
            samples,
            data: ndarray::Array3::<i32>::zeros((ns,7,4)),
        }
    }

    pub fn update(&mut self, sidx: usize, m: usize, gt: &[i32], alleles: &Vec<&[u8]>) {
        // this ignores -1/1 GTs, will see how it affects performance
        //if !is_variant(gt) {
        //    return a
        //}
        let (is_var, vidx) = is_variant(gt);
        if !is_var { return }
    
        // add to n_vars
        self.data[[sidx, NVARS, m]] += 1;
    
        
        if is_snv(decode_genotype(gt[vidx]), &alleles) {
            // add to n_snvs
            self.data[[sidx, NSNVS, m]] += 1;
            if is_ti(alleles[0], alleles[decode_genotype(gt[vidx]) as usize]) {
                // add to n_tis
                self.data[[sidx, NTIS, m]] += 1
            } else {
                // add to n_tvs
                self.data[[sidx, NTVS, m]] += 1
            }
        } else {
            // add to n_indels
            self.data[[sidx, NINDELS, m]] += 1;
        }
    
        if is_hom(gt) {
            // add to n_hom_alts
            self.data[[sidx, NHOMALTS, m]] += 1
        } else {
            // add to n_hets
            self.data[[sidx, NHETS, m]] += 1
        }
    }
    
    pub fn cat(&mut self, metrics: Metrics) {
        if self.gq != metrics.gq {
            panic!("GQ thresholds do not match");
        }
        if self.freq != metrics.freq {
            panic!("Rarity thresholds do not match");
        }
        if self.population != metrics.population {
            panic!("The population AF source does not match");
        }
        if self.samples != metrics.samples {
            panic!("Samples do not match");
        }
        self.data += &metrics.data;
    }
    
    pub fn to_csv(self, iowtr: Box<dyn std::io::Write>) {
        let mut wtr = csv::WriterBuilder::new().has_headers(false).from_writer(iowtr);
        _ = wtr.write_record(&[
            "sample", "n_vars", "n_snvs", "n_indels", "n_hets", "n_hom_alts", "n_tis", "n_tvs",
            &format!("gq{}_n_vars", self.gq), &format!("gq{}_n_snvs", self.gq), &format!("gq{}_n_indels", self.gq), &format!("gq{}_n_hets", self.gq), &format!("gq{}_n_hom_alts", self.gq), &format!("gq{}_n_tis", self.gq), &format!("gq{}_n_tvs", self.gq),
            &format!("{}{}_n_vars", self.population, self.freq), &format!("{}{}_n_snvs", self.population, self.freq), &format!("{}{}_n_indels", self.population, self.freq), &format!("{}{}_n_hets", self.population, self.freq), &format!("{}{}_n_hom_alts", self.population, self.freq), &format!("{}{}_n_tis", self.population, self.freq), &format!("{}{}_n_tvs", self.population, self.freq),
            &format!("gq{}_{}{}_n_vars", self.gq, self.population, self.freq), &format!("gq{}_{}{}_n_snvs", self.gq, self.population, self.freq), &format!("gq{}_{}{}_n_indels", self.gq, self.population, self.freq), &format!("gq{}_{}{}_n_hets", self.gq, self.population, self.freq), &format!("gq{}_{}{}_n_hom_alts", self.gq, self.population, self.freq), &format!("gq{}_{}{}_n_tis", self.gq, self.population, self.freq), &format!("gq{}_{}{}_n_tvs", self.gq, self.population, self.freq),
        ]);
        for sidx in 0..self.samples.len() {
            _ = wtr.write_record(&[
                &self.samples[sidx], &format!("{}", self.data[[sidx, NVARS, RAW]]), &format!("{}", self.data[[sidx, NSNVS, RAW]]), &format!("{}", self.data[[sidx, NINDELS, RAW]]), &format!("{}", self.data[[sidx, NHETS, RAW]]), &format!("{}", self.data[[sidx, NHOMALTS, RAW]]), &format!("{}", self.data[[sidx, NTIS, RAW]]), &format!("{}", self.data[[sidx, NTVS, RAW]]),
                &format!("{}", self.data[[sidx, NVARS, QUAL]]), &format!("{}", self.data[[sidx, NSNVS, QUAL]]), &format!("{}", self.data[[sidx, NINDELS, QUAL]]), &format!("{}", self.data[[sidx, NHETS, QUAL]]), &format!("{}", self.data[[sidx, NHOMALTS, QUAL]]), &format!("{}", self.data[[sidx, NTIS, QUAL]]), &format!("{}", self.data[[sidx, NTVS, QUAL]]),
                &format!("{}", self.data[[sidx, NVARS, RARE]]), &format!("{}", self.data[[sidx, NSNVS, RARE]]), &format!("{}", self.data[[sidx, NINDELS, RARE]]), &format!("{}", self.data[[sidx, NHETS, RARE]]), &format!("{}", self.data[[sidx, NHOMALTS, RARE]]), &format!("{}", self.data[[sidx, NTIS, RARE]]), &format!("{}", self.data[[sidx, NTVS, RARE]]),
                &format!("{}", self.data[[sidx, NVARS, QUAL_RARE]]), &format!("{}", self.data[[sidx, NSNVS, QUAL_RARE]]), &format!("{}", self.data[[sidx, NINDELS, QUAL_RARE]]), &format!("{}", self.data[[sidx, NHETS, QUAL_RARE]]), &format!("{}", self.data[[sidx, NHOMALTS, QUAL_RARE]]), &format!("{}", self.data[[sidx, NTIS, QUAL_RARE]]), &format!("{}", self.data[[sidx, NTVS, QUAL_RARE]]),
            ]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_decode_genotype() {
        for (e, d) in std::iter::zip([0,2,4,6,8], [-1,0,1,2,3])  {
            assert_eq!(decode_genotype(e), d)
        }
    }
    #[test]
    fn test_is_snv() {
        let alleles_v: Vec<Vec<&[u8]>> = vec![vec![b"A", b"C"]];
        for alleles in alleles_v {
            assert!(is_snv(1, &alleles));
        }
    }
    #[test]
    fn test_not_snv() {
        let alleles_v: Vec<Vec<&[u8]>> = vec![vec![b"A", b"AC"], vec![b"ACT", b"A"], vec![b"ACC", b"AAA"]];
        for alleles in alleles_v {
            assert!(!is_snv(1, &alleles));
        }
    }
    #[test]
    fn test_is_ti() {
        let alleles_v = vec![vec![b"A", b"G"], vec![b"C", b"T"], vec![b"G", b"A"], vec![b"T", b"C"]];
        for alleles in alleles_v {
            assert!(is_ti(alleles[0], alleles[1]))
        }
    }
    #[test]
    fn test_not_ti() {
        let alleles_v = vec![
            vec![b"A", b"T"], vec![b"A", b"C"],
            vec![b"C", b"A"], vec![b"C", b"G"],
            vec![b"G", b"C"], vec![b"G", b"T"],
            vec![b"T", b"A"], vec![b"T", b"G"]
        ];
        for alleles in alleles_v {
            assert!(!is_ti(alleles[0], alleles[1]))
        }
    }
    #[test]
    fn test_is_hom() {
        let gt_v = vec![vec![0, 0], vec![1, 1], vec![2, 2], vec![3, 3]];
        for gt in gt_v {
            assert!(is_hom(&gt))
        }
    }
    #[test]
    fn test_is_not_hom() {
        let gt_v = vec![vec![0, 1], vec![1, 0], vec![2, 1], vec![-1, 0]];
        for gt in gt_v {
            assert!(!is_hom(&gt))
        }
    }
}