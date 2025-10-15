use anyhow::{bail, Context, Result};
use bio::alphabets::dna;
use bio::io::fasta;
use clap::{ArgAction, Parser};
use csv::ReaderBuilder;
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::PathBuf;

/// Conventional FASTA wrap width.
const FASTA_WRAP_WIDTH: usize = 60;

/// Re-orient FASTA reads to a constant direction using a TSV with per-read orientation.
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Cli {
    /// Input FASTA (can be .fa/.fasta(.gz)); use '-' for stdin (plain text, not gz)
    #[arg(short = 'f', long)]
    fasta: String,

    /// Tab-delimited table with headers (can be .tsv/.txt(.gz))
    #[arg(short = 't', long)]
    table: PathBuf,

    /// Output FASTA path (default: stdout)
    #[arg(short = 'o', long)]
    out: Option<PathBuf>,

    /// Name of the read ID column in the table
    #[arg(long, default_value = "ReadName")]
    id_col: String,

    /// Name of the orientation column in the table ('+' for cDNA, '-' for rc(cDNA))
    #[arg(long, default_value = "orientation")]
    orientation_col: String,

    /// Target orientation to keep as-is; reads not matching are reverse-complemented. Allowed: '+' or '-'
    #[arg(long, default_value = "+")]
    target_orientation: String,

    /// If true, drop reads missing in the table (instead of passing through unchanged)
    #[arg(long, action = ArgAction::SetTrue)]
    drop_missing: bool,

    /// Append a suffix to headers of flipped reads (e.g., '/rc'); empty = no suffix
    #[arg(long, default_value = "")]
    flipped_suffix: String,
}

fn open_writer(path: &Option<PathBuf>) -> Result<Box<dyn Write>> {
    Ok(match path {
        Some(p) => Box::new(BufWriter::new(
            File::create(p).with_context(|| format!("create {:?}", p))?,
        )),
        None => Box::new(BufWriter::new(io::stdout())),
    })
}

fn open_text(path: &str) -> Result<Box<dyn Read>> {
    if path == "-" {
        // stdin (expect plain text; if gz, pipe through zcat/gunzip externally)
        return Ok(Box::new(io::stdin()));
    }
    let fh = File::open(path).with_context(|| format!("open '{}'", path))?;
    if path.ends_with(".gz") {
        Ok(Box::new(MultiGzDecoder::new(fh)))
    } else {
        Ok(Box::new(fh))
    }
}

fn load_orientation_map(cli: &Cli) -> Result<HashMap<String, u8>> {
    // Support gz TSV by looking at extension.
    let rdr: Box<dyn Read> = if cli.table.to_string_lossy().ends_with(".gz") {
        Box::new(MultiGzDecoder::new(
            File::open(&cli.table).with_context(|| format!("open {:?}", &cli.table))?,
        ))
    } else {
        Box::new(File::open(&cli.table).with_context(|| format!("open {:?}", &cli.table))?)
    };

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(rdr);

    let headers = reader.headers().context("reading TSV headers")?.clone();

    let id_idx = headers
        .iter()
        .position(|h| h == &cli.id_col)
        .with_context(|| format!("column '{}' not found", cli.id_col))?;
    let ori_idx = headers
        .iter()
        .position(|h| h == &cli.orientation_col)
        .with_context(|| format!("column '{}' not found", cli.orientation_col))?;

    let mut map = HashMap::with_capacity(1 << 16);
    for rec in reader.records() {
        let rec = rec?;
        let id = rec.get(id_idx).unwrap().to_owned();
        let ori_field = rec.get(ori_idx).unwrap().trim().as_bytes();
        // Accept '+', '-', or words starting with those.
        let ori = if !ori_field.is_empty() {
            match ori_field[0] as char {
                '+' => b'+',
                '-' => b'-',
                _ => {
                    let s = String::from_utf8_lossy(ori_field).to_ascii_lowercase();
                    if s.starts_with("plus") || s.starts_with("fwd") || s == "1" {
                        b'+'
                    } else if s.starts_with("minus") || s.starts_with("rev") || s == "0" || s == "rc" {
                        b'-'
                    } else {
                        bail!("Unrecognized orientation value '{}' for read '{}'", s, id);
                    }
                }
            }
        } else {
            bail!("Empty orientation for read '{}'", id);
        };
        map.insert(id, ori);
    }
    Ok(map)
}

fn wrap_and_write<W: Write>(w: &mut W, seq: &[u8]) -> Result<()> {
    for chunk in seq.chunks(FASTA_WRAP_WIDTH) {
        w.write_all(chunk)?;
        w.write_all(b"\n")?;
    }
    Ok(())
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let target = match cli.target_orientation.as_str() {
        "+" => b'+',
        "-" => b'-',
        other => bail!("--target-orientation must be '+' or '-', got '{}'", other),
    };

    let ori_map = load_orientation_map(&cli).context("loading orientation table")?;
    let mut out = open_writer(&cli.out)?;

    // Open FASTA (plain or gz). Use '-' to read from stdin (plain).
    let handle = open_text(&cli.fasta)?;
    let reader = fasta::Reader::new(handle);

    let mut n_total: u64 = 0;
    let mut n_flipped: u64 = 0;
    let mut n_missing: u64 = 0;

    for result in reader.records() {
        let record = result.context("parsing FASTA record")?;
        n_total += 1;

        let id = record.id().to_string();
        let desc = record.desc().unwrap_or("");
        let mut header = id.clone();
        if !desc.is_empty() {
            header.push(' ');
            header.push_str(desc);
        }

        // Decide action
        let action = match ori_map.get(&id) {
            Some(&ori) => {
                if ori == target { "keep" } else { "flip" }
            }
            None => {
                if cli.drop_missing {
                    n_missing += 1;
                    continue; // skip this record
                } else {
                    "keep"
                }
            }
        };

        // Sequence handling
        let mut seq = record.seq().to_vec();
        if action == "flip" {
            n_flipped += 1;
            seq = dna::revcomp(&seq);
            if !cli.flipped_suffix.is_empty() {
                header.push_str(&cli.flipped_suffix);
            }
        }

        // Emit FASTA with wrapping
        writeln!(out, ">{}", header)?;
        wrap_and_write(&mut out, &seq)?;
    }

    // Progress to stderr
    eprintln!(
        "processed={} flipped={} missing_in_table={} ({} mode) | wrap={} cols",
        n_total,
        n_flipped,
        n_missing,
        if cli.drop_missing { "dropped" } else { "kept" },
        FASTA_WRAP_WIDTH
    );

    Ok(())
}
