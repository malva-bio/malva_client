# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- Opt-in Leiden grouping: `search()`, `search_genes()` and `search_federated()`
  accept `group_by="leiden"` (default `None`, which omits the field and keeps
  the legacy request unchanged). `SearchResult.grouping` reports which
  dimension the server actually returned (`grouping["applied"]` is the
  authoritative value, including the `cell_type` fallback and its
  `fallback_reason`), and `SearchResult.cluster_annotations` exposes the
  per-(sample, cluster) annotation table. Cluster-keyed frames gain a
  `cluster_id` column whose `cell_type` is resolved from those annotations.
  Federated searches expose the aggregate block in `df.attrs["grouping"]` and
  the per-provider maps in `df.attrs["cluster_annotations"]`. `retrieve_cells`
  is unchanged.

## [0.7.0] - 2026-11-11
### Added
- `find_metacells(method=...)` supports the optional `"hellinger"` and
  `"panel_lda"` server scoring spaces while keeping `"cosine"` as the exact
  legacy default.
- `MetacellMatchResult.assign_labels()` now uses server per-cell-type top-k
  scores for the new spaces (`weight="auto"`), and still supports the legacy
  similarity/rank neighbour aggregation via `weight="similarity"`/`"rank"`.

## [0.5.0] - 2026-11-09
### Added
- `stream_cell_expression()` returns a `CellExpressionStream` that yields
  sparse per-cell batches and converts them with `to_anndata()`,
  `iter_anndata()` or `write_h5ad()`, optionally enriched with sample metadata.
- `SearchResult.to_anndata()` builds a cells × genes AnnData object (obs are
  `sample_id` × `cell_type` pairs) with provenance in `uns["malva"]`.
- `find_metacells()` and the `MetacellMatchResult` wrapper map local clusters to
  Malva metacells, with `assign_labels()` writing cell-type labels back to the
  input AnnData. `compare_clusters_to_metacells()` remains as a deprecated
  alias.
- `differential_expression()`, `differential_expression_from_cells()` and
  `DifferentialExpressionResult` expose server-computed DE tables.
- `project_cells_to_metacells()` projects an AnnData object onto metacells.
- Workspace/provider context on `MalvaClient` (`workspace`, `dataset_id`,
  `provider_id`), `set_workspace()`, `get_workspace()` and `get_datasets()`.
  Config file/env support uses `MALVA_WORKSPACE`, `MALVA_DATASET_ID` and
  `MALVA_PROVIDER_ID`.
- `provider=` scoping on sample/metadata methods, sent as the `?provider=`
  query parameter. Explicit arguments override the client default.
- `search_federated()` with per-provider provenance in `df.attrs`, visible
  failed-shard errors, an optional `strict=True`, and a stable empty-result
  schema.
- Study discovery falls back to the dataset hierarchy when
  `/samples/api/studies` or the per-dataset studies endpoint is empty.
- Canonical Explorer endpoints are preferred, with automatic one-time fallback
  to the historical paths. `MALVA_ENDPOINT_STYLE=auto|canonical|legacy` and
  `client.endpoint_style` expose the resolved style.
- `retrieve_cells()` accepts the DataFrame returned by `search_federated()` and
  retrieves cells from a completed shard, optionally selected with `provider=`.
- `scripts/live_check.py` runs the whole public client surface against a
  deployed server, including a Scanpy reference query, and writes a PASS/FAIL
  report. The `live-check` CI job uploads that report as a build artifact.

### Fixed
- `retrieve_cells()` now forwards an explicitly selected federated shard's
  provider id as `?provider=` to every search-service retrieval call and skips
  the local values export, so remote shard jobs resolve on their own service.
- Federated searches without explicit `provider_ids` now resolve the default
  provider list from `get_workspace()` before submitting, so canonical
  `POST /api/search` is dispatched as a federated job. Failed or empty
  workspace lookups omit the field and log a warning.
- `retrieve_cells()` now treats cell-expression values-export timeouts as
  recoverable and falls back through the full-search, cell-ID, global-ID and
  Explorer fetchers while preserving the final aggregated error.
- Unknown federation provider ids are rejected before submission instead of
  producing a zero-shard job that stays `running`.
- A zero-shard federated job that never leaves `running` now raises a clear
  error after the poll budget.
- `get_sample_details()` surfaces the server error body and provider context.
- `get_study_samples()` raises instead of silently returning an empty list when
  the hierarchy reports samples that the provider route cannot return.

## [0.4.0] - 2026-10-26
### Added
- Cell Explorer client methods: build, list, preview, delete, download and load
  per-user metacell cohorts stored as Zarr. The cohort schema is
  `malva.cell_explorer.v1` (`obs/*`, `var/gene`, `X`, `layers/raw`, optional
  `obsm/X_umap`).
- `malva_client cell-explorer` CLI commands for listing, building and
  downloading cohorts.
- Large Explorer responses above the API Gateway 10 MB cap are followed
  transparently when the Explorer returns a `303` redirect to a presigned S3
  URL.

## [0.3.4] - 2026-06-05
### Added
- Coverage APIs accept `stage_callback` and log separate search, server-result generation, download, and parse stages for large coverage responses.
- `get_cells_by_metadata()` for search-independent all-cell denominator retrieval by sample and/or cell type. Database-wide retrieval is available with explicit `include_all_database_cells=True`.
- Aggregate result columns now use a single Expression Explorer-style schema: `rel`, `exp`, `pct`, `raw_kmers`, and `cell_count`. Legacy expression names are still accepted by helper methods such as `aggregate_by(..., expr_column=...)` when they can be mapped unambiguously.

### Changed
- Quickstart and tutorials now avoid placeholder sample IDs and document the meaning of aggregate and per-cell expression values.
- `retrieve_cells()` documentation clarifies that denominator cells should be fetched independently with `get_cells_by_metadata()`, while fallback positive-cell-only endpoints may report `value=1`.
- All-cell metadata decoding uses categorical cell types and rounded/clipped `uint16` total counts for lower memory use and faster DataFrame construction.

### Fixed
- `retrieve_cells()` now prefers the direct msgpack per-cell expression values endpoint, which uses the same raw k-mer value resolver as the Explorer matrix export without generating a ZIP. Aggregate searches therefore return raw per-cell values instead of all `1.0` positive indicators when values are available.
- Removed stale docs examples that produced errors with placeholder IDs or non-existent tutorial links.

## [0.3.2] - 2026-05-07
### Added
- `retrieve_cells()` for direct per-cell retrieval from an existing search job,
  including cell IDs, feature metadata, sample metadata, and per-cell
  normalization factors when the metadata endpoint is available.
- `project_cells()` and sample-scoped `get_coexpression(..., filter_sample_ids=...)`
  support for projecting selected cells or selected samples onto coexpression
  indices.

### Changed
- Cell-level workflows now use the two-step search-then-retrieve path
  consistently across docs, examples, and tests.

## [0.3.1] - 2026-02-27
### Bugfix
- [Internal] An issue where the availability of lazy-loaded data was not checked before loading the dataframe

## [0.3.0] - 2026-02-26

### Added
- **`min_kmer_presence` / `max_kmer_presence`** parameters on `search()`,
  `search_sequences()`, and `search_genes()` - filter k-mers
  by database prevalence to reduce noise from rare errors or ubiquitous
  repeats
- **Unified `search_sequences()`** now accepts either a single `str` or a
  `List[str]`.  A single string (up to 500 kb) is sent directly; a list
  (total ≤ 100 kb) is submitted as a FASTA batch.  `search_sequence()` is
  removed; pass a string to `search_sequences()` instead

### Changed
- **All search methods submit jobs asynchronously** (`wait_for_completion`
  is now handled by client-side polling only).  This prevents 502 Bad
  Gateway errors from nginx proxy read timeouts on long-running searches
- `stranded` parameter now correctly maps to server's `unstranded_mode`
  field (previously the key was sent verbatim and ignored by the server)
- `window_size` and `threshold` parameters removed - these were server-
  internal and never had any effect on results

### Fixed
- **502 Bad Gateway on long searches**: POST `/search` no longer holds
  the connection open; the server starts the job and returns a `job_id`
  immediately
- **"still running" after Redis write failure**: `GET /search/{job_id}`
  now falls back to the on-disk result file when Redis holds a stale
  `running`/`pending` entry (happens when the Redis connection resets
  after a long job)
- **Empty results for all searches**: `SearchResult._convert_to_dataframe`
  now handles the `_fmt: col` columnar expression format returned by the
  server (previously every gene was silently skipped)
- **Non-JSON error responses**: `_request` now catches
  `requests.exceptions.JSONDecodeError` (simplejson variant) and reports
  the HTTP status code and body instead of a cryptic "Request failed" error

## [0.2.0] - 2026-02-06

### Added
- **Search parameters**: `window_size`, `threshold`, and `stranded` options on `search()`, `submit_search()`, and all methods that delegate to them
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
