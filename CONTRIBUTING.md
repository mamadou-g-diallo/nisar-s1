# Contributing

Issues and pull requests are welcome.

1. Install for development: `pip install -e ".[dev]"`.
2. Keep functions generic (no site-specific paths); put data and credentials outside the repository.
3. Add or update a test in `tests/` (synthetic data, no network).
4. Run `ruff check src tests` and `pytest` before opening a pull request.
5. Describe in the pull request what changes for users, and add a line to `CHANGELOG.md`.

By contributing you agree that your contribution is licensed under the Apache License 2.0.
