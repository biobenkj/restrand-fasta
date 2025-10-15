# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] - 2025-10-15
### Added
- Cross-architecture compile compatibility
- Basic tests added

## [0.2.0] - 2025-10-14
### Added
- Fixed **60-column** FASTA wrapping (no CLI override)
- GitHub Actions **CI**: fmt, clippy, tests on Ubuntu/macOS/Windows
- GitHub Actions **Release**: macOS **universal2** binary; Linux (x86_64, aarch64) builds
- `.cargo/config.toml` with ThinLTO + symbol stripping for release
- Makefile targets (`universal`, `build`, `fmt`, `clippy`, `test`)

### Changed
- Improved TSV orientation parsing (accepts `plus/fwd/1` and `minus/rev/0/rc`).

## [0.1.0] - 2025-10-14
### Added
- Initial CLI to re-orient FASTA sequences using per-read orientation table.
