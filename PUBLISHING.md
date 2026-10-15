# Publishing to PyPI

## Versioning

The version lives in exactly two places and they must match:

- `pyproject.toml` → `[project].version` (the source of truth),
- `malva_client/__init__.py` → the `__version__ = "x.y.z"` fallback (used only
  when the package metadata is unavailable, e.g. running from a plain checkout).

The project follows [Semantic Versioning](https://semver.org/) and
[Keep a Changelog](https://keepachangelog.com/) (`CHANGELOG.md`).

To cut a release:

1. On a branch, bump `version` in `pyproject.toml` and the fallback in
   `malva_client/__init__.py`, and add a `CHANGELOG.md` entry.
2. Merge to `main`.
3. Tag and push: `git tag vX.Y.Z && git push origin vX.Y.Z` (or run the
   `release-pypi` workflow manually from the Actions tab).

`release-pypi.yml` fails if the pushed tag does not match the `pyproject.toml`
version, so a release can never publish a version other than the one tagged.

## Releasing from CI (recommended)

`.github/workflows/release-pypi.yml` builds the sdist and the wheel and uploads
everything to PyPI. The build itself lives in
`.github/workflows/build-wheels.yml`, which also runs on every pull request
(without publishing) so a PR shows whether the package builds before it is
merged.

Trigger it either way:

- GitHub → **Actions** → **release-pypi** → **Run workflow** (leave *Publish to
  PyPI* checked), or
- push a tag: `git tag v0.3.5 && git push origin v0.3.5`.

Run it from a green `main`. PyPI rejects a version that already exists, so a
release is one-shot per version bump.

### One-time authentication setup

Pick one of the two options below.

**Trusted publishing (recommended — nothing to rotate).** On
<https://pypi.org/manage/account/publishing/> add a *pending publisher*:

| Field | Value |
|---|---|
| PyPI project name | `malva-client` |
| Owner | `malva-bio` |
| Repository name | `malva_client` |
| Workflow name | `release-pypi.yml` |
| Environment name | `pypi` |

The workflow already requests `id-token: write` and targets the `pypi`
environment, so that is all.

**API token.** Create a project token at
<https://pypi.org/manage/account/token/>, add it as the repository secret
`PYPI_API_TOKEN`, then in the `publish` job swap the trusted-publishing step for:

```yaml
    steps:
      - uses: actions/download-artifact@v4
        with:
          name: dist
          path: dist
      - uses: pypa/gh-action-pypi-publish@release/v1
        with:
          password: ${{ secrets.PYPI_API_TOKEN }}
```

## Integration tests (real platform)

`tests/test_integration.py` hits the real Malva platform to check that the
configured endpoints work. It is skipped unless `MALVA_API_TOKEN` is set. It
discovers a dataset id from the metadata hierarchy for the UMAP/coexpression
checks (those need a dataset whose coexpression index is loaded).

- **Locally**: export the token (the client already reads `MALVA_API_TOKEN` and
  `MALVA_API_URL`), e.g. from the macOS keychain or a file outside the repo, then
  run `pytest -m integration`.
- **In CI**: add a repository **secret** `MALVA_API_TOKEN` and a repository
  **variable** `MALVA_API_URL` under *Settings → Secrets and variables →
  Actions*. The `integration` job passes them as environment variables and skips
  automatically when the token is absent (e.g. fork pull requests).

## Publishing locally

Keep this for a dry run or an offline release:

```bash
pip install build twine
python -m build
python -m twine upload dist/*
```
