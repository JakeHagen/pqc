pub fn cat(input: &Option<Vec<String>>, list: Option<&str>, output: Option<&str>) {
    let iowtr: Box<dyn std::io::Write> = match output {
        None => Box::new(std::io::stdout()),
        Some("-") => Box::new(std::io::stdout()),
        Some(o) => Box::new(std::fs::File::create(o).expect("unable to create output file"))
    };

    let mut files:Vec<String> = vec![];
    match list {
        None => {
            match input {
                None => files.push("-".to_owned()),
                Some(fv) => {
                    for f in fv {
                        if f == "-" {
                            if fv.len() > 1 {
                                panic!("can not read from stdin and line")
                            }
                            files.push("-".to_owned())
                        }
                        files.push(f.to_owned());
                    }
                }
            }
        },
        Some(fl) => {
            let json_list = std::io::BufReader::new(std::fs::File::open(fl).expect("couldnt open list of jsons"));
            for json_file_r in std::io::BufRead::lines(json_list) {
                let json_file = json_file_r.expect("could not get line from json list file");
                files.push(json_file)
            }
        }
    }

    if files[0] == "-" {
        if files.len() > 1 { panic!("can not read from stdin and line") }
        let rdr = std::io::BufReader::new(std::io::stdin());
        let metrics: pqc::Metrics = serde_json::from_reader(rdr).expect("");
        metrics.to_csv(iowtr);
    } else {
        let mut count: i32 = 0;
        let rdr = std::io::BufReader::new(std::fs::File::open(&files[0]).expect("couldnt open input file"));
        let mut metrics: pqc::Metrics = serde_json::from_reader(rdr).expect("");
        for f in &files[1..] {
            let rdr = std::io::BufReader::new(std::fs::File::open(f).expect("couldnt open input file"));
            let m: pqc::Metrics = serde_json::from_reader(rdr).expect("");
            metrics.cat(m);
            count += 1;
            if count.rem_euclid(100) == 0 {
                eprintln!("{} done", count);
            }
        }
        metrics.to_csv(iowtr);
    }
}