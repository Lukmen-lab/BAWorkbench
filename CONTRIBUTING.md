# Contributing

Contributions are welcome through GitHub issues and pull requests.

## Development setup

1. Fork and clone the repository.
2. Create the environment with `setup_env.bat` (Windows) or `./setup_env.sh` (macOS/Linux).
3. Create a focused branch.
4. Add or update tests with each behavior change.
5. Run `python -m pytest -q` before opening a pull request.

## Chemistry changes

Changes to bile-acid cores, stereochemistry, attachment sites, reaction rules, or partner structures should include:

- the structure source or literature reference;
- an explanation of stereochemical assumptions;
- representative structure and exact-mass tests;
- updated limitations or warning flags where uncertainty remains.

Do not commit runtime `outputs/`, conformer caches, confidential data, or unpublished experimental results.

By submitting a contribution, you agree that it may be distributed under the MIT License.
