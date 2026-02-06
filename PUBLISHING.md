## Publishing to PyPI

```bash
pip install build twine
python -m build
python -m twine upload dist/*
```