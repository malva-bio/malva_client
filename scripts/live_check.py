#!/usr/bin/env python3
"""Run an end-to-end live check of the Malva client against a deployed server.

The script exercises the public client surface in order, including a reference
query: a small Scanpy example dataset is clustered locally and only the cluster
pseudobulk profiles are sent to the server. Every step is recorded as PASS,
FAIL or SKIP and the result is written as Markdown and JSON so CI can publish it
as a build artifact.

Example:
    export MALVA_API_TOKEN="..."
    python scripts/live_check.py --out reports
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import time
import traceback
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

import requests

from malva_client import MalvaClient
from malva_client.exceptions import MalvaAPIError

DEFAULT_URL = "https://malva.mdc-berlin.de"
DEFAULT_GENE = "CD3D"


@dataclass
class Check:
    name: str
    status: str
    seconds: float = 0.0
    details: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None
    traceback: Optional[str] = None


def _trim(value: Any, limit: int = 240) -> str:
    text = str(value)
    return text if len(text) <= limit else text[: limit - 3] + "..."


def describe(value: Any) -> Dict[str, Any]:
    """Summarise a client return value without dumping large payloads."""
    if value is None:
        return {"value": None}
    try:
        import pandas as pd
    except ImportError:
        pd = None
    if pd is not None and isinstance(value, pd.DataFrame):
        return {
            "type": "DataFrame",
            "rows": int(value.shape[0]),
            "columns": [str(c) for c in value.columns][:20],
            "head": json.loads(value.head(3).to_json(orient="records", date_format="iso")),
        }
    if isinstance(value, dict):
        summary: Dict[str, Any] = {"type": "dict", "keys": [str(k) for k in list(value)[:20]]}
        preview = {
            str(k): _trim(v)
            for k, v in list(value.items())[:8]
            if isinstance(v, (str, int, float, bool, type(None)))
        }
        if preview:
            summary["preview"] = preview
        for key in ("datasets", "studies", "samples", "results", "clusters", "cells"):
            item = value.get(key)
            if isinstance(item, (list, tuple)):
                summary[f"n_{key}"] = len(item)
        return summary
    if isinstance(value, (list, tuple)):
        summary = {"type": type(value).__name__, "n": len(value)}
        if value and isinstance(value[0], dict):
            summary["first_keys"] = [str(k) for k in list(value[0])[:20]]
        return summary
    frame = getattr(value, "df", None)
    if pd is None or not isinstance(frame, pd.DataFrame):
        dataframe = getattr(value, "to_dataframe", None)
        if callable(dataframe):
            try:
                frame = dataframe()
            except Exception:  # pragma: no cover - defensive
                frame = None
    if pd is not None and isinstance(frame, pd.DataFrame):
        summary = {
            "type": type(value).__name__,
            "rows": int(frame.shape[0]),
            "columns": [str(c) for c in frame.columns][:20],
            "head": json.loads(frame.head(3).to_json(orient="records", date_format="iso")),
        }
        if "rank" in frame.columns:
            best = frame[frame["rank"] == 1]
            summary["top_match_per_cluster"] = json.loads(
                best.head(25).to_json(orient="records", date_format="iso")
            )
        if "cell_type" in frame.columns:
            counts = frame["cell_type"].astype(object).fillna("(none)").value_counts().head(10)
            summary["cell_type_counts"] = {str(k): int(v) for k, v in counts.items()}
        return summary
    attributes = {}
    for attr in ("ncells", "n_cells", "n_genes", "n_obs", "n_vars", "shape", "job_id"):
        if hasattr(value, attr):
            try:
                candidate = getattr(value, attr)
                attributes[attr] = (
                    candidate if isinstance(candidate, (int, float, str)) else _trim(candidate)
                )
            except Exception:  # pragma: no cover - defensive
                pass
    if attributes:
        return {"type": type(value).__name__, **attributes}
    return {"type": type(value).__name__, "repr": _trim(value)}


class SkipCheck(Exception):
    """Raise from a check to record it as skipped instead of failed."""


_TRANSIENT_ERRORS = (
    "503",
    "504",
    "no available server",
    "service unavailable",
    "gateway time",
    "connection",
    "timeout",
)


def retrying(fn: Callable[[], Any], attempts: int = 3, delay: float = 8.0) -> Any:
    """Retry a call whose failure looks like a temporary server outage."""
    for attempt in range(1, attempts + 1):
        try:
            return fn()
        except Exception as exc:
            transient = isinstance(exc, MalvaAPIError) and any(
                token in str(exc).lower() for token in _TRANSIENT_ERRORS
            )
            if attempt == attempts or not transient:
                raise
            print(f"  retrying after {type(exc).__name__}: {_trim(exc)}")
            time.sleep(delay * attempt)
    return None


class Recorder:
    def __init__(self) -> None:
        self.checks: List[Check] = []

    def record(self, name: str, fn: Callable[[], Any], skip: Optional[str] = None) -> Any:
        if skip is not None:
            self.checks.append(Check(name=name, status="SKIP", details={"reason": skip}))
            print(f"SKIP {name}: {skip}")
            return None
        started = time.monotonic()
        caught: List[warnings.WarningMessage] = []
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                value = fn()
        except SkipCheck as exc:
            self.checks.append(Check(name=name, status="SKIP", details={"reason": str(exc)}))
            print(f"SKIP {name}: {exc}")
            return None
        except Exception as exc:
            self.checks.append(
                Check(
                    name=name,
                    status="FAIL",
                    seconds=time.monotonic() - started,
                    error=f"{type(exc).__name__}: {exc}",
                    traceback=traceback.format_exc(limit=6),
                )
            )
            print(f"FAIL {name}: {type(exc).__name__}: {exc}")
            return None
        details = describe(value)
        messages = [f"{type(item.message).__name__}: {item.message}" for item in caught]
        if messages:
            details["warnings"] = [_trim(message, 400) for message in messages[:10]]
        check = Check(
            name=name,
            status="PASS",
            seconds=time.monotonic() - started,
            details=details,
        )
        self.checks.append(check)
        print(f"PASS {name} ({check.seconds:.2f}s)")
        return value

    @property
    def result(self) -> str:
        if any(check.status == "FAIL" for check in self.checks):
            return "NO PASS"
        if not any(check.status == "PASS" for check in self.checks):
            return "SKIPPED"
        return "PASS"


def select_dataset(client: MalvaClient) -> Optional[str]:
    try:
        hierarchy = client.get_datasets_hierarchy()
    except MalvaAPIError:
        hierarchy = {}
    candidates: List[str] = []
    for dataset in hierarchy.get("datasets") or []:
        if not isinstance(dataset, dict):
            continue
        dataset_id = dataset.get("dataset_id") or dataset.get("id")
        if dataset_id and str(dataset_id) not in candidates:
            candidates.append(str(dataset_id))
    if candidates:
        return candidates[0]
    try:
        for dataset in client.get_datasets():
            dataset_id = dataset.get("dataset_id") or dataset.get("id")
            if dataset_id:
                return str(dataset_id)
    except MalvaAPIError:
        return None
    return None


def scanpy_query_dataset() -> Any:
    """Load a small, bundled Scanpy example dataset for the reference query."""
    import scanpy as sc

    adata = sc.datasets.pbmc68k_reduced()
    if adata.raw is not None:
        adata = adata.raw.to_adata()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    groupby = "bulk_labels" if "bulk_labels" in adata.obs else adata.obs.columns[0]
    adata.obs[groupby] = adata.obs[groupby].astype("category")
    return adata, groupby


def write_report(
    out_dir: Path,
    server: str,
    token_set: bool,
    recorder: Recorder,
    dataset_id: Optional[str],
    result: Optional[str] = None,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    headline = result or recorder.result
    payload = {
        "server": server,
        "token_set": token_set,
        "result": headline,
        "dataset_id": dataset_id,
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "platform": {"python": platform.python_version(), "machine": platform.platform()},
        "checks": [
            {
                "name": check.name,
                "status": check.status,
                "seconds": round(check.seconds, 3),
                "details": check.details,
                "error": check.error,
            }
            for check in recorder.checks
        ],
    }
    (out_dir / "live_check.json").write_text(json.dumps(payload, indent=2, default=str))

    lines = [
        "# Malva client live check",
        "",
        f"- **Result**: {headline}",
        f"- **Server**: {server}",
        f"- **Token configured**: {'yes' if token_set else 'no'}",
        f"- **Dataset**: {dataset_id or '(none discovered)'}",
        f"- **Generated**: {payload['generated_at']}",
        f"- **Platform**: Python {payload['platform']['python']} on {payload['platform']['machine']}",
        "",
        "| Check | Status | Time (s) | Note |",
        "| --- | --- | --- | --- |",
    ]
    for check in recorder.checks:
        note = check.error or check.details.get("reason") or ""
        lines.append(
            f"| {check.name} | {check.status} | {check.seconds:.2f} | {_trim(note, 160)} |"
        )
    lines += ["", "## Details", ""]
    for check in recorder.checks:
        lines.append(f"### {check.name} — {check.status}")
        lines.append("")
        if check.error:
            lines += ["```", check.error, "```", ""]
        if check.traceback and check.status == "FAIL":
            lines += [
                "<details><summary>traceback</summary>",
                "",
                "```",
                check.traceback.strip(),
                "```",
                "",
                "</details>",
                "",
            ]
        if check.details:
            lines += ["```json", json.dumps(check.details, indent=2, default=str)[:4000], "```", ""]
    report = out_dir / "live_check.md"
    report.write_text("\n".join(lines) + "\n")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default=os.environ.get("MALVA_API_URL") or DEFAULT_URL)
    parser.add_argument("--token", default=os.environ.get("MALVA_API_TOKEN"))
    parser.add_argument("--out", default="reports", help="Directory for the Markdown/JSON report.")
    parser.add_argument("--dataset", default=None, help="Dataset id; discovered when omitted.")
    parser.add_argument("--gene", default=DEFAULT_GENE, help="Gene symbol used for search checks.")
    parser.add_argument("--top-n", type=int, default=25, help="Metacell neighbours per cluster.")
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument(
        "--skip-reference-query",
        action="store_true",
        help="Skip the Scanpy reference query (useful when scanpy is unavailable).",
    )
    args = parser.parse_args()

    token_set = bool(args.token)
    recorder = Recorder()

    recorder.record(
        "health",
        lambda: requests.get(
            args.url.rstrip("/") + "/health", timeout=min(args.timeout, 30)
        ).json(),
    )

    if not token_set:
        reason = "MALVA_API_TOKEN is not set"
        for name in (
            "authenticate",
            "datasets_hierarchy",
            "studies",
            "overview_stats",
            "available_filters",
            "samples",
            "gene_search",
            "retrieve_cells",
            "reference_query",
        ):
            recorder.record(name, lambda: None, skip=reason)
        write_report(Path(args.out), args.url, token_set, recorder, None, result="SKIPPED")
        return 0

    client = recorder.record(
        "client_init",
        lambda: MalvaClient(base_url=args.url, api_token=args.token, timeout=args.timeout),
    )
    if client is None:
        report = write_report(Path(args.out), args.url, token_set, recorder, args.dataset)
        print(f"report: {report}")
        print(f"result: {recorder.result}")
        return 1

    recorder.record("authenticate", client.get_quota_status)
    hierarchy = recorder.record("datasets_hierarchy", client.get_datasets_hierarchy)
    recorder.record("studies", client.get_studies)
    recorder.record("overview_stats", client.get_overview_stats)
    recorder.record("available_filters", client.get_available_filters)
    recorder.record("samples", client.get_samples)

    dataset_id = args.dataset or select_dataset(client)
    if dataset_id is None and isinstance(hierarchy, dict):
        for dataset in hierarchy.get("datasets") or []:
            if isinstance(dataset, dict) and (dataset.get("dataset_id") or dataset.get("id")):
                dataset_id = str(dataset.get("dataset_id") or dataset.get("id"))
                break

    search_result = recorder.record("gene_search", lambda: client.search(args.gene))
    recorder.record(
        "retrieve_cells",
        lambda: client.retrieve_cells(search_result, features=args.gene),
        skip=None if search_result is not None else "gene search did not complete",
    )

    if dataset_id is None:
        for name in ("umap_coordinates", "search_to_anndata", "coexpression"):
            recorder.record(name, lambda: None, skip="no dataset id discovered")
    else:

        def _umap() -> Any:
            try:
                return client.get_umap_coordinates(dataset_id)
            except MalvaAPIError as exc:
                if "not found" in str(exc).lower():
                    raise SkipCheck("the dataset has no UMAP coordinates") from exc
                raise

        recorder.record("umap_coordinates", _umap)
        recorder.record(
            "search_to_anndata",
            lambda: search_result.to_anndata() if search_result is not None else None,
            skip=None if search_result is not None else "gene search did not complete",
        )

        def _coexpression() -> Any:
            frame = search_result.df
            sample_ids: List[int] = []
            if "rel" in frame.columns and not frame.empty:
                top = frame.sort_values("rel", ascending=False)["sample_id"].drop_duplicates()
                sample_ids = [int(value) for value in top.head(3)]
            return retrying(
                lambda: client.get_coexpression(
                    search_result.job_id,
                    dataset_id,
                    filter_sample_ids=sample_ids or None,
                    max_wait=240,
                )
            )

        recorder.record(
            "coexpression",
            _coexpression,
            skip=None if search_result is not None else "gene search did not complete",
        )

    if args.skip_reference_query:
        recorder.record("reference_query", lambda: None, skip="--skip-reference-query")
    elif dataset_id is None:
        recorder.record("reference_query", lambda: None, skip="no dataset id discovered")
    else:
        adata, groupby = scanpy_query_dataset()
        recorder.record(
            "reference_query",
            lambda: retrying(
                lambda: client.find_metacells(
                    adata,
                    groupby=groupby,
                    dataset_id=dataset_id,
                    top_n=args.top_n,
                )
            ),
        )

    report = write_report(Path(args.out), args.url, token_set, recorder, dataset_id)
    print(f"report: {report}")
    print(f"result: {recorder.result}")
    return 0 if recorder.result in {"PASS", "SKIPPED"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
