pub fn cat(input: &Option<Vec<String>>, output: Option<&str>) {
    let iowtr: Box<dyn std::io::Write> = match output {
        None => Box::new(std::io::stdout()),
        Some("-") => Box::new(std::io::stdout()),
        Some(o) => Box::new(std::fs::File::open(o).expect("unable to create output file"))
    };

    let mut iordrs:Vec<Box<dyn std::io::Read>> = vec![];
    match input {
        None => iordrs.push(Box::new(std::io::BufReader::new(std::io::stdin()))),
        Some(rs) => {
            for r in rs {
                if r == "-" {
                    if rs.len() > 1 {
                        panic!("can not read from stdin and line")
                    }
                    iordrs.push(Box::new(std::io::BufReader::new(std::io::stdin())))
                } else {
                    iordrs.push(Box::new(std::io::BufReader::new(std::fs::File::open(r).expect("couldnt open input file"))))
                }
            }
        }
    }

    let mut metrics: pqc::Metrics = serde_json::from_reader(&mut iordrs[0]).expect("");
    for rdr in iordrs[1..].iter_mut() {
        let m: pqc::Metrics = serde_json::from_reader(rdr).expect("");
        metrics.cat(m)
    }

    metrics.to_csv(iowtr);

}