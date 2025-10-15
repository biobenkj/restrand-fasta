# restrand-fasta

[![CI](https://github.com/biobenkj/restrand-fasta/actions/workflows/ci.yml/badge.svg)](https://github.com/biobenkj/restrand-fasta/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/biobenkj/restrand-fasta?display_name=tag&sort=semver)](https://github.com/biobenkj/restrand-fasta/releases)

### Re-orient ONT reads based on molecule architecture table

`restrand-fasta` flips FASTA sequences to a constant molecular orientation using a TSV with `ReadName` and `orientation` (`+` = cDNA, `-` = rc(cDNA)). Output FASTA is always wrapped to 60 characters as a common convention/compatibility.

## Install

### From source
```
cargo install --git https://github.com/biobenkj/restrand-fasta
```

### Download binaries
Grab macOS (universal) and Linux builds from [Releases](https://github.com/biobenkj/restrand-fasta/releases).

## Usage

```
restrand-fasta \
  -f reads.fa.gz \
  -t annotations.tsv \
  --target-orientation + \
  --flipped-suffix "/rc" \
  > reoriented.fa
```

- Input FASTA can be gzipped; TSV can be gzipped.  
- Reads absent from the table pass through unchanged; add `--drop-missing` to drop them.  
- Header is preserved; if flipped, optional suffix is appended.

## TSV columns expected

- `ReadName` (string, must match FASTA IDs)  
- `orientation` (`+` for cDNA, `-` for rc(cDNA); also accepts `plus/fwd/1` and `minus/rev/0/rc`)

### Example TSV

```
ReadName	orientation
09232f67-71ed-46be-9518-da3561e0351e	-
ea1e3068-7185-4dea-8473-b5e1c3624b19	+
eab28f58-f4db-463f-a4a9-4d2c354ea54f	-
```

## MSRV

- Minimum Supported Rust Version: **1.74**

## License

Dual-licensed under **MIT** or **Apache-2.0**.  
You may choose either license.  
See [`LICENSE-MIT`](LICENSE-MIT) and [`LICENSE-APACHE`](LICENSE-APACHE) for details.

## Security

If you discover a security vulnerability in `restrand-fasta`, please **do not** open a public issue.  
Instead, email **Ben Johnson** at the address listed on the GitHub profile for [@biobenkj](https://github.com/biobenkj).  
We will coordinate a fix and disclosure.
