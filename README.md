# AMULET

AMULET is a gene-first workflow for inferring gene modules and module-level
relationships from single-cell RNA-seq data. The method constructs gene-gene
relation signals from mutual-exclusion profiles, trains soft gene-to-hyperedge
assignments, and uses the learned gene modules to support cell identity
interpretation.

This repository is a compact release package for review and reproduction. It
contains the core implementation, a minimal simulation dataset, and one runnable
example script. Generated results, notebooks, exploratory scripts, and temporary
files are intentionally excluded.

## Repository Structure

```text
.
├── run_306.py                  # Minimal end-to-end example script
├── CME_CPU.py                  # Numba helper for CME computation
├── datasets/
│   └── adata_306.h5ad          # Small simulation example
├── src/
│   ├── cme_supervision.py      # CME, p-value, and supervision mask construction
│   ├── gene_ambiguity.py       # Garbage-hyperedge ambiguity score
│   ├── train.py                # Training loss and optimization loop
│   ├── supervision_pipeline.py # Gene-to-hyperedge assignment pipeline
│   ├── metagene_tree.py        # Module-level directed inclusion aggregation
│   ├── cell_type_evaluation.py # Prototype-style cell-type evaluation
│   ├── cme_visualization.py    # Diagnostic heatmaps
│   └── loss.py                 # Relation and hierarchy loss utilities
├── environment.yml             # Minimal conda environment
└── requirements-minimal.txt    # Minimal pip dependencies
```

## Installation

Create a Python 3.11 environment and install the required packages:

```bash
conda create -n amulet python=3.11
conda activate amulet
pip install -r requirements-minimal.txt
```

Alternatively:

```bash
conda env create -f environment.yml
conda activate amulet
```

## Input Data

The workflow expects an `.h5ad` file with:

- `adata.X`: cell-by-gene expression matrix.
- `adata.var_names`: gene names.
- optional `adata.obs["cell_type"]` or a similar label column for evaluation.
- optional gene module information encoded in gene names or `adata.var` for
  simulation diagnostics.

The included `datasets/adata_306.h5ad` is only a minimal demonstration dataset.
For new data, place the `.h5ad` file under `datasets/` and update `ADATA_PATH`
near the top of `run_306.py`.

## Running the Workflow

Run the minimal demo:

```bash
python run_306.py
```

By default, outputs are written to:

```text
results/<dataset_name>_results/
```

For a new dataset, copy `run_306.py` to a new script and adjust:

- `ADATA_PATH`
- `t_CME`, `t_p`, `t_Jaccard`, `t_inclusion`
- `num_biological_modules`
- `num_garbage_modules`
- `training_epochs`

## Method Steps

1. Compute the conditional mutual exclusion (CME) matrix from gene expression.
2. Estimate empirical p-values for CME scores by permutation.
3. Construct relation supervision masks:
   - positive gene pairs from Jaccard similarity of mutual-exclusion profiles;
   - negative gene pairs from significant CME;
   - weak-positive inclusion pairs from directed inclusion profiles.
4. Compute gene ambiguity from relation contradictions and bidirectional
   inclusion conflicts.
5. Train soft gene-to-hyperedge assignments. The garbage hyperedge receives
   ambiguity-weighted genes and is excluded from ordinary relation loss.
6. Convert soft assignments to discrete gene modules by `argmax`.
7. Aggregate gene-level directed inclusion into module-level inclusion edges.
8. If cell-type labels are available, compute module activity prototypes and
   evaluate cell assignment.

## Main Outputs

- `summary.json`: run configuration and key metrics.
- `gene_assignment_diagnostics.csv`: per-gene module probabilities and ambiguity.
- `gene_modules.csv`: genes grouped by learned hyperedge.
- `garbage_genes_only.csv`: genes assigned to the garbage hyperedge.
- `module_inclusion_selected_edges.csv`: directed module-level inclusion edges.
- `module_cell_type_table.csv`: dominant cell type for each module.
- `cell_assignment_from_modules.csv`: prototype-based cell assignment.
- `training_loss_history.csv/png`: optimization trajectory.
- `combined_supervision_heatmaps.png`: CME and supervision diagnostics.
- `hyperedge_run_summary.png`: learned gene-to-hyperedge assignment summary.

## Minimal Demo

The included 306 simulation file is a small demo used to verify that the workflow
can be executed end-to-end. It is not intended to be the main benchmark dataset.
Running `python run_306.py` regenerates all demo outputs under
`results/adata_306_results/`.

## Reproducibility Notes

- The example script fixes `seed = 0`.
- Results are regenerated locally and are not included in the release package.
- The package contains no notebooks or intermediate exploratory outputs.
