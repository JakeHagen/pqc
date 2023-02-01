mod commands;
use commands::*;
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    QC {
        input: Option<String>,
        #[clap(long, short)]
        output: Option<String>,
        #[clap(long, short)]
        field: String,
        #[clap(long, short)]
        rare: f32,
        #[clap(long, short)]
        gq: i32,
        #[clap(long, value_parser, default_value_t = 1)]
        threads: usize,
    },
    Cat {
        input: Option<Vec<String>>,
        #[clap(long, short)]
        output: Option<String>,
    }
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Commands::QC { input, output, field, rare, gq, threads } => {
            qc::qc(input.as_deref(), output.as_deref(), field, rare, gq, threads)
        },
        Commands::Cat { input, output } => {
            cat::cat(input, output.as_deref())
        },
    }
}
