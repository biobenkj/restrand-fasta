use anyhow::{bail, Context, Result};
use bio::alphabets::dna;
use bio::io::{fasta, fastq};
use clap::{ArgAction, Parser};
use csv::ReaderBuilder;
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::PathBuf;

/// Conventional FASTA wrap width.
const FASTA_WRAP_WIDTH: usize = 60;

/// Re-orient FASTA/FASTQ reads to a constant direction using a TSV with per-read orientation or embedded orientation tags.
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Cli {
    /// Input FASTA/FASTQ (can be .fa/.fasta/.fq/.fastq(.gz)); use '-' for stdin (plain text, not gz)
    #[arg(short = 'f', long)]
    fasta: String,

    /// Tab-delimited table with headers (can be .tsv/.txt(.gz)); not required for --fastq mode
    #[arg(short = 't', long)]
    table: Option<PathBuf>,

    /// Output path (default: stdout)
    #[arg(short = 'o', long)]
    out: Option<PathBuf>,

    /// Process as FASTQ and read orientation from header (looks for 'orientation:+' or 'orientation:-')
    #[arg(long, action = ArgAction::SetTrue)]
    fastq: bool,

    /// Name of the read ID column in the table (FASTA mode only)
    #[arg(long, default_value = "ReadName")]
    id_col: String,

    /// Name of the orientation column in the table ('+' for cDNA, '-' for rc(cDNA)) (FASTA mode only)
    #[arg(long, default_value = "orientation")]
    orientation_col: String,

    /// Target orientation to keep as-is; reads not matching are reverse-complemented. Allowed: '+' or '-'
    #[arg(long, default_value = "+")]
    target_orientation: String,

    /// If true, drop reads missing in the table (instead of passing through unchanged) (FASTA mode only)
    #[arg(long, action = ArgAction::SetTrue)]
    drop_missing: bool,

    /// Append a suffix to headers of flipped reads (e.g., '/rc'); empty = no suffix (FASTA mode only)
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

fn load_orientation_map(table_path: &PathBuf, id_col: &str, orientation_col: &str) -> Result<HashMap<String, u8>> {
    // Support gz TSV by looking at extension.
    let rdr: Box<dyn Read> = if table_path.to_string_lossy().ends_with(".gz") {
        Box::new(MultiGzDecoder::new(
            File::open(table_path).with_context(|| format!("open {:?}", table_path))?,
        ))
    } else {
        Box::new(File::open(table_path).with_context(|| format!("open {:?}", table_path))?)
    };

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(rdr);

    let headers = reader.headers().context("reading TSV headers")?.clone();

    let id_idx = headers
        .iter()
        .position(|h| h == id_col)
        .with_context(|| format!("column '{}' not found", id_col))?;
    let ori_idx = headers
        .iter()
        .position(|h| h == orientation_col)
        .with_context(|| format!("column '{}' not found", orientation_col))?;

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

/// Extract orientation from FASTQ header (looks for "orientation:+" or "orientation:-")
fn extract_orientation_from_header(header: &str) -> Option<u8> {
    if let Some(start) = header.find("orientation:") {
        let rest = &header[start + "orientation:".len()..];
        if let Some(first_char) = rest.chars().next() {
            return match first_char {
                '+' => Some(b'+'),
                '-' => Some(b'-'),
                _ => None,
            };
        }
    }
    None
}

/// Update header to change orientation:- to orientation:+
fn update_orientation_in_header(header: &str) -> String {
    header.replace("orientation:-", "orientation:+")
}

fn wrap_and_write<W: Write>(w: &mut W, seq: &[u8]) -> Result<()> {
    for chunk in seq.chunks(FASTA_WRAP_WIDTH) {
        w.write_all(chunk)?;
        w.write_all(b"\n")?;
    }
    Ok(())
}

fn process_fastq(cli: &Cli, target: u8) -> Result<()> {
    let handle = open_text(&cli.fasta)?;
    let reader = fastq::Reader::new(handle);
    let mut out = open_writer(&cli.out)?;

    let mut n_total: u64 = 0;
    let mut n_flipped: u64 = 0;
    let mut n_no_orientation: u64 = 0;

    for result in reader.records() {
        let record = result.context("parsing FASTQ record")?;
        n_total += 1;

        let id = record.id().to_string();
        let desc = record.desc().unwrap_or("");
        let mut header = id.clone();
        if !desc.is_empty() {
            header.push(' ');
            header.push_str(desc);
        }

        // Extract orientation from header
        let full_header = if desc.is_empty() {
            id.as_str()
        } else {
            header.as_str()
        };

        let ori = extract_orientation_from_header(full_header);

        let mut seq = record.seq().to_vec();
        let mut qual = record.qual().to_vec();
        let mut output_header = header.clone();

        match ori {
            Some(o) if o != target => {
                // Need to flip
                n_flipped += 1;
                seq = dna::revcomp(&seq);
                qual.reverse(); // Reverse quality scores to match reversed sequence
                output_header = update_orientation_in_header(&output_header);
            }
            Some(_) => {
                // Already at target orientation, keep as-is
            }
            None => {
                // No orientation tag found, keep as-is
                n_no_orientation += 1;
            }
        }

        // Write FASTQ record
        writeln!(out, "@{}", output_header)?;
        out.write_all(&seq)?;
        out.write_all(b"\n")?;
        writeln!(out, "+")?;
        out.write_all(&qual)?;
        out.write_all(b"\n")?;
    }

    eprintln!(
        "FASTQ mode: processed={} flipped={} no_orientation_tag={}",
        n_total, n_flipped, n_no_orientation
    );

    Ok(())
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let target = match cli.target_orientation.as_str() {
        "+" => b'+',
        "-" => b'-',
        other => bail!("--target-orientation must be '+' or '-', got '{}'", other),
    };

    // If FASTQ mode, process as FASTQ
    if cli.fastq {
        return process_fastq(&cli, target);
    }

    // FASTA mode requires a table
    let table = cli.table.as_ref()
        .context("--table is required for FASTA mode (or use --fastq for FASTQ mode)")?;

    let ori_map = load_orientation_map(table, &cli.id_col, &cli.orientation_col)
        .context("loading orientation table")?;
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
        "FASTA mode: processed={} flipped={} missing_in_table={} ({} mode) | wrap={} cols",
        n_total,
        n_flipped,
        n_missing,
        if cli.drop_missing { "dropped" } else { "kept" },
        FASTA_WRAP_WIDTH
    );

    Ok(())
}
