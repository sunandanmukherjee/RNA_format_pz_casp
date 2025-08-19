# rna-submission-prep

<p align="center">
  <a href="https://github.com/sunandanmukherjee/RNA_format_pz_casp"><img src="https://img.shields.io/badge/GitHub-Repo-blue.svg" alt="GitHub Repo"></a>
  <a href="https://github.com/sunandanmukherjee/RNA_format_pz_casp/blob/main/LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
  <img src="https://img.shields.io/badge/python-3.7+-blue.svg" alt="Python 3.7+">
  <img src="https://img.shields.io/badge/dependency-Biopython-green.svg" alt="Dependency: Biopython">
</p>

A command-line pipeline for formatting 3D RNA models to meet the strict submission guidelines of the [RNA-Puzzles](https://www.rnabase.org/puzzles/) and [CASP](https://predictioncenter.org/) experiments.

## 🧬 Overview

Participating in community-wide experiments like RNA-Puzzles and CASP requires submitted 3D models to adhere to a specific format. This often involves matching atom order, residue numbering, and chain identifiers to an official template PDB file. Manually editing PDB files to meet these requirements is tedious, time-consuming, and prone to errors.

**RNA_format_pz_casp** automates this entire process with a single command. It provides a simple workflow to split, reformat, and combine your 3D models according to an official template, ensuring your submissions are valid and correctly processed.

## ✨ Key Features

-   **One-Command Pipeline**: Automates splitting, formatting, and combining models into a final submission file.
-   **Rank-Based Model Handling**: Processes multiple input models in your specified ranked order.
-   **Template-Based Formatting**: Automatically reorders atoms and adjusts residue/chain information to match an official template.
-   **Clean Workflow**: Intermediate files are automatically removed, keeping your workspace tidy.
-   **Flexible**: Supports both multi-model and single-model PDB files as input.

## 🚀 Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/sunandanmukherjee/rna-submission-prep.git
    cd rna-submission-prep
    ```

2.  **Install dependencies:**
    This project requires Biopython. Install it using the provided `requirements.txt` file.
    ```bash
    pip install -r requirements.txt
    ```

## 📖 Recommended Usage: The Automated Pipeline

The `formatPDB_pipeline.py` script is the easiest and recommended way to format your models for submission. It handles all steps automatically.

### Command

```bash
python formatPDB_pipeline.py -t <template.pdb> -m <model1.pdb> <model2.pdb> ... -o <output.pdb>
```

### Arguments

-   `-t`, `--template` (Required): Path to the official template PDB file provided by RNA-Puzzles/CASP.
-   `-m`, `--models` (Required): A space-separated list of your model PDB files, ordered by rank (best model first).
-   `-o`, `--output` (Optional): The name for the final, combined submission PDB file. If omitted, a name will be generated based on the template file (e.g., `formatted_submission_template.pdb`).
-   `-k`, `--keep-intermediate` (Optional): A flag to prevent the script from deleting the intermediate split and formatted files. Useful for debugging.

---

### 📝 Walkthrough Example

Let's say you are participating in an RNA-Puzzles challenge.

1.  **Get Files**:
    -   You download the official template: `PZ_TARGET_template.pdb`.
    -   Your prediction program generates your top 3 models: `rank1.pdb`, `rank2.pdb`, and `rank3.pdb`.

2.  **Run the Pipeline**:
    Execute a single command to process all models and create a final submission file named `MySubmission.pdb`.

    ```bash
    python formatPDB_pipeline.py \
      --template PZ_TARGET_template.pdb \
      --models rank1.pdb rank2.pdb rank3.pdb \
      --output MySubmission.pdb
    ```

3.  **Submit**:
    The script will create `MySubmission.pdb`, which contains your three models, correctly formatted and combined with `MODEL`/`ENDMDL` records. This file is now ready for upload to the submission server.

### 📁 Examples Directory

This repository includes an `examples/` directory containing a sample template and model files so you can test the pipeline immediately.

To run the example:
```bash
# Navigate to the root directory of the repository
python formatPDB_pipeline.py \
  --template examples/template.pdb \
  --models examples/model_set.pdb \
  --output examples/formatted_example_submission.pdb
```
This will generate a properly formatted submission file inside the `examples` directory.

---

## 🔬 Manual Workflow (Step-by-Step)

If you need more granular control or want to debug a specific step, you can use the individual scripts.

#### Step 1: Split a Multi-Model PDB File

If your prediction tool outputs a single PDB file containing multiple models, use `split_PDB.py` to separate them.

```bash
python split_PDB.py <your_multi_model_file.pdb> --model
```
This will create `your_multi_model_file_model1.pdb`, `_model2.pdb`, etc.

#### Step 2: Format a Single Model

Use `format_PDB.py` to adjust a single model file based on the official template.

```bash
python format_PDB.py --template <template.pdb> --model <your_single_model.pdb>
```
This creates a new file named `formatted_your_single_model.pdb`.

## 📜 License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## 🙌 Contributing

Contributions are welcome! If you have suggestions for improving this tool, feel free to open an issue or submit a pull request.
