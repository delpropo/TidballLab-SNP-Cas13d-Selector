# Copilot Instructions for TidballLab-SNP-Cas13d-Selector

## 1. Getting Started

- Ensure you have Python (≥3.8), Conda, and Snakemake installed.
- Clone the repository:
  ```bash
  git clone https://github.com/delpropo/TidballLab-SNP-Cas13d-Selector.git
  cd TidballLab-SNP-Cas13d-Selector
  ```
- Set up the environment:
  ```bash
  conda env create -f environment.yml
  conda activate tidball_snp_selector
  ```
- Install Python dependencies (if needed):
  ```bash
  pip install -r requirements.txt
  ```

## 2. Running the Workflow

- Configure the workflow in `config/config.yaml`.
- Perform a dry run:
  ```bash
  snakemake --dry-run
  ```
- Run the workflow:
  ```bash
  snakemake --cores 4
  ```
- For containerized execution, see the README for Apptainer/Singularity instructions.

## 3. Project Structure

- `src/tidball_snp_selector/`: Main Python package
- `workflow/`: Snakemake files (Snakefile, rules, envs, scripts)
- `config/`: Configuration and schemas
- `data/`: Data storage
- `tests/`: Test suite
- `docs/`: Documentation (MkDocs)

## 4. Testing

- Run all tests:
  ```bash
  pytest
  ```
- Add new tests in the `tests/` directory.

## 5. Documentation

- User and developer docs are in `docs/`.
- To view locally (if using MkDocs):
  ```bash
  mkdocs serve
  ```

## 6. Contributing

- Follow PEP8 and project coding standards.
- Open issues or pull requests for bugs, features, or improvements.
- See `CONTRIBUTING.md` if available.

## 7. Troubleshooting

- Check environment setup and config files.
- Review error messages and logs.
- For help, open an issue on GitHub.

---

For more details, see the main `README.md` and documentation in `docs/`.




