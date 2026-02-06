# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-02-06

### Added
- **Search parameters**: `window_size`, `threshold`, and `stranded` options on `search()`, `search_cells()`, `submit_search()`, and all methods that delegate to them
- **Batch sequence search**: `search_sequences()` method for querying multiple DNA sequences at once
- **Coverage analysis**: `get_coverage()`, `get_sequence_coverage()`, `get_coverage_data()`, `download_coverage_wig()`, and `get_coverage_filter_options()` methods on `MalvaClient`
- **`CoverageResult` model**: with `to_dataframe()`, `plot()`, `download_wig()`, and `get_filter_options()` methods
- **Dataset discovery**: `get_datasets_hierarchy()`, `get_dataset_studies()`, `get_study_samples()`, `get_sample_details()`, `get_filter_values()`, and `get_overview_stats()` methods
- **CLI flags**: `--window-size`, `--threshold`, `--stranded` for the `search` command
- **Query parameters guide**: new documentation page with recommended settings per use case
- **Tutorials section** in documentation with notebook integration
- **CHANGELOG.md** for tracking releases
- `sphinx-design` and `nbsphinx` extensions for documentation

### Fixed
- **Job management methods** (`submit_search`, `get_job_status`, `get_job_results`, `list_jobs`, `wait_for_job`, `cancel_job`) were defined as standalone functions instead of `MalvaClient` methods, causing `NameError` at runtime
- Removed non-existent `get_coverage_for_region` from `__init__.py` exports

### Changed
- Version is now derived from `importlib.metadata` with `pyproject.toml` as single source of truth
- Documentation overhauled with grid cards layout, notebook integration, and correct URLs
- Updated all repository and documentation URLs to `malva-bio` organization
- Documentation: removed separate examples section; all content consolidated into quickstart guide

## [0.1.0] - 2025-01-01

### Added
- Initial release
- `MalvaClient` with search, sample browsing, and download functionality
- `SearchResult` and `SingleCellResult` models
- CLI with `search`, `config`, `login`, and `quota` commands
- Sphinx documentation with Furo theme
