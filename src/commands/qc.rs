use rust_htslib::bcf::Read;
use std::io;
use std::str;
use std::fs::File;
use std::io::BufWriter;


pub fn qc(
    input: Option<&str>,
    output: Option<&str>,
    field: &String,
    rare: &f32,
    gq: &i32,
    threads: &usize,
) {
    let mut bcf = bcflib::get_rdr(input);
    bcf.set_threads(threads.clone())
        .expect("unable to set reader threads");

    let hdrv = bcf.header().clone();
    let n_samples = usize::try_from(hdrv.sample_count()).unwrap();
    let samples = hdrv.samples();
    let mut string_samples:Vec<String> = vec![];
    for s in samples {
        string_samples.push(String::from_utf8(s.to_vec()).unwrap())
    }

    //let mut metrics = pqc::Metrics::new(string_samples, *gq, *rare, field.to_string());
    let mut metrics = pqc::Metrics::new(string_samples, *gq, *rare, field.to_string());

    for record_result in bcf.records() {
        let record = record_result.expect("fail to read record");
        let pf = pqc::get_popfreq(&record, &field);
        let alleles = record.alleles();

        let gqs = record.format(b"GQ").integer().expect("Couldn't retrieve GQ field");
        let gts = record.format(b"GT").integer().expect("Couldn't retrieve GT field");

        for sidx in 0..n_samples {
            metrics.update(sidx, pqc::RAW, gts[sidx], &alleles);
            if pf <= *rare {
                metrics.update(sidx, pqc::RARE, gts[sidx], &alleles);
            }
            if gqs[sidx][0] >= *gq {
                metrics.update(sidx, pqc::QUAL, gts[sidx], &alleles);
            }
            if pf <= *rare && gqs[sidx][0] >= *gq {
                metrics.update(sidx, pqc::QUAL_RARE, gts[sidx], &alleles);
            }
        }
    }
    
    let iowtr: Box<dyn io::Write> = match output {
        None => Box::new(io::stdout()),
        Some("-") => Box::new(io::stdout()),
        Some(o) => Box::new(File::create(o).expect("unable to create output file"))
    };
    let mut wtr = BufWriter::new(iowtr);
    serde_json::to_writer_pretty(&mut wtr, &metrics).expect("could not write json");

}