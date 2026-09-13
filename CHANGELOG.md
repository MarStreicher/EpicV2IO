# Changelog

All notable changes to this project will be documented in this file.


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

