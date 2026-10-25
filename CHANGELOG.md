# Changelog

All notable changes to this project will be documented in this file.


## [Unreleased]

### Added

### Changed

### Fixed

## [1.0.0] - 2026-07-17

### Added
- Initial public release
- Reading and writing EPIC V2 files

### Changed
- Improved documentation

### Fixed
- Packaging metadata

## [1.0.1] - 2026-07-17

### Added
- Nothing

### Changed
- Change package name in toml

### Fixed
- Packaging metadata

## [1.0.2] - 2026-09-13

### Changed
- Delimiter is taken from the file extension (`.csv` → comma, `.txt` → tab); content sniffing and semicolon support removed
- `BetasLoader` loads betas and the Peters manifest lazily (only when needed)

## [1.0.3] - 2026-10-25

### Changed
- Split `BetasLoader` into I/O helpers (`io.py`) and Peters QC exclusions (`exclusions.py`); public loader API is unchanged
- README now leads with loading (`example/usage_epicv2-io.py`); manifest usage is documented last
- Releases are tag-driven (`v*.*.*` → `release.yml`); the package version is taken from the git tag (setuptools-scm) 

## [1.0.4] - 2026-10-25

### Changed
-  Hotfix: Fix PyPI publish workflow and README

