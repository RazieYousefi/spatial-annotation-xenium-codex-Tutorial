# Project: Cell-type annotation in spatial omics, verified by Gaussian KDE

## Part 1 — Project description (instructions for Claude)

### Goal
Demonstrate a complete spatial omics analysis on two public datasets, written as Jupyter notebooks:

1. **Annotate cell types myself** with the standard workflow: QC → normalization → PCA →
   neighbor graph → Leiden clustering → marker genes/proteins → manual naming of clusters.
2. **Validate spatially:** my labels on the tissue, plus neighborhood enrichment.
3. **Verify every cell with a Gaussian KDE** (scikit-learn `KernelDensity`): a density per cell
   type in PCA space; each cell is checked against the density of its own type and against
   competing types.
4. **Compare my annotation with the original publications:** cell-level agreement with the
   authors' labels, cell-type composition, and whether I recover their main findings.

The deliverable is a set of clean, readable notebooks that show command of spatial omics
analysis, with a density-based check of how much to trust each cell's label.

### Current status
Only the data is in place. No code has been written yet.

### Data (in `data/`, do not re-download)
| | Xenium (transcriptomics) | CODEX (proteomics) |
|---|---|---|
| Original work | Janesick et al., Nat Commun 2023 | Schürch et al., Cell 2020 |
| Files | `data/xenium_breast/`: `xe1_cell_feature_matrix.h5` (cell × gene counts), `xe1_cells.csv.gz` (centroids, QC), `Cell_Barcode_Type_Matrices.xlsx` (authors' labels) | `data/codex_crc/CRC_clusters_neighborhoods_markers.csv` (one table: marker intensities, coordinates, patient, authors' labels) |
| Tissue | FFPE breast cancer, replicate 1 | CRC invasive front, 35 patients, 140 regions |
| Features | 313-gene panel, counts | 56 protein markers, intensities |
| Authors' labels | sheet "Xenium R1 Fig1-5 (supervised)" (`Barcode`, `Cluster`) | column `ClusterName` |
| Sources | Zenodo 10076046 | Mendeley 10.17632/mpjzbtfgfr.1 |

Both licensed CC BY 4.0. Notes for loading:
- Xenium: keep only features of type "Gene Expression" (drop negative controls / blanks);
  join matrix barcodes to `cell_id` in `xe1_cells.csv.gz` and to `Barcode` in the label sheet.
- CODEX: marker columns contain `:Cyc_`; drop nuclear stains and empty channels
  (HOECHST, DRAQ5, empty, blank); coordinates `X:X`, `Y:Y`; patient column `patients`.

**Authors' labels are never used to build my annotation or the KDE.** They enter only in
notebook 04.

### Notebooks
```
helpers.py                          small shared functions only (loading, KDE wrapper, plotting)
01_xenium_annotation.ipynb          load → QC → normalize → PCA → Leiden → markers → my labels → spatial check
02_codex_annotation.ipynb           same workflow for CODEX
03_kde_verification.ipynb           KDE check of my labels, both datasets
04_comparison_with_original.ipynb   my labels vs. the authors', both datasets
results/                            saved AnnData (.h5ad), tables, figures
```
Notebooks pass results through `results/*.h5ad` so each one can be rerun on its own.

### What each notebook does
**01 / 02 — My annotation**
- Load into AnnData (counts or intensities in `.X`, coordinates in `.obsm["spatial"]`).
- QC: Xenium — drop cells with very few transcripts; CODEX — drop obvious debris.
- Normalize: Xenium — library size + log1p; CODEX — arcsinh. Then scale.
- PCA → kNN graph → Leiden. Resolution chosen by marker coherence, recorded in the notebook.
- Markers per cluster (rank test) and a dot plot of canonical markers.
- **My annotation:** a cluster → cell-type dictionary written by hand in the notebook, with the
  marker evidence for each name in a comment. Unclear clusters are `unknown`.
- Spatial check: cell-type maps (CODEX: one region per patient group) and neighborhood
  enrichment (squidpy).
- Save `results/<dataset>_annotated.h5ad`.

**03 — KDE verification** (both datasets)
- Space: the first k PCs used for clustering.
- KDE: `sklearn.neighbors.KernelDensity`, Gaussian kernel. The PCs have different spreads and
  `KernelDensity` uses one bandwidth for all directions, so each type's points are **whitened**
  first (centered and scaled by that type's covariance), which is equivalent to a bandwidth
  matrix H = h² Σ̂. The whitening term is added back to the log-density so types stay comparable.
- Bandwidth per type: Scott's rule as reference; final h by held-out log-likelihood on a small grid.
- At most N cells per type as KDE support (random subsample) to keep it laptop-fast.
- Cross-fitting so no cell scores itself: split cells into two halves by spatial block (Xenium,
  ~500 µm tiles) or by patient (CODEX); fit on one half, score the other, swap.
- Per cell with my label c (all in log space, `scipy.special.logsumexp`):
  - log p(x | c) below the q-th percentile of typical type-c cells → **outlier**
  - posterior p(c | x) ∝ π_c p(x | c) clearly favors another type → **possibly mislabeled**
  - otherwise → **confirmed**
- A short toy example first (two 2-D classes): densities at three bandwidths, the decision
  boundary and the outlier region, so the reader sees what the check does before real data.
- Outputs: flag rates per type, flag map on the tissue, examples of flagged cells' marker
  profiles, sensitivity to k, bandwidth, N, q.

**04 — Comparison with the original work** (both datasets)
- Name mapping table between my types and the authors' (granularity differs; written explicitly).
- Cell-level agreement: confusion matrix, adjusted Rand index, per-type agreement.
- Composition: my vs. the authors' cell-type proportions.
- Main findings: do my labels reproduce the papers' key observations?
  - Xenium: tumor compartments (DCIS vs. invasive regions) and immune populations shown in the paper.
  - CODEX: immune composition differences between patient groups (CLR vs. DII) and
    lymphoid aggregates in CLR patients.
- KDE link: are KDE-flagged cells enriched where my labels and the authors' disagree?
- Short written conclusion at the end.

### Coding rules
- Notebooks read top to bottom as a narrative: a markdown cell explaining each step, then the code.
- Simple and minimal: scanpy, squidpy, numpy, pandas, scipy, scikit-learn, matplotlib.
  Nothing else unless the goal requires it.
- `helpers.py` holds only functions used in more than one notebook. No classes, no frameworks.
- Parameters in one cell at the top of each notebook; fixed random seeds.
- Short docstrings; math in markdown where it helps understanding.
- The manual annotation is my decision: code produces evidence, I write the cluster names.

### Out of scope
Cell segmentation, transcript-level data, images, segmentation-free spatial KDE (possible later
extension), cross-platform reference mapping, deep-learning or automated annotation tools,
hyperparameter searches beyond Leiden resolution, k, bandwidth, N and q.

### Environment
```
pip install scanpy squidpy igraph leidenalg h5py openpyxl scikit-learn matplotlib jupyter
```

---

## Part 2 — Project understanding

**What this project shows.** That I can take spatial omics data from raw matrices to a defensible
cell-type map on two very different platforms, check it against the published work, and
quantify how much to trust each label instead of presenting the map as certain.

**Why manual annotation.** Clustering plus marker-based naming is how most labs annotate spatial
data. Doing it myself, with the marker evidence recorded, demonstrates the core skill. Its
weakness is that it labels clusters, not cells: every cell in a cluster inherits the name,
including doublets, poorly segmented cells and transitional states.

**Why scikit-learn for the KDE.** The estimator is standard; a library implementation gives the
same numbers as a hand-written one and is faster. The project's contribution is how the KDE is
used (per-cell verification of a manual annotation, cross-fitting, linking flags to
disagreements), not the estimator itself. The one decision the library does not make for me,
the bandwidth shape, is handled explicitly by whitening.

**What the KDE adds.** It moves the check from the cluster to the cell. Each cell type becomes a
smooth density in expression space. A cell far from the bulk of its own type's density, or one
another type explains better, is flagged. Every flag has a number behind it (a density and a
posterior), so the check is interpretable.

**Why cross-fitting.** A cell that helps build its own type's density will always look typical.
Fitting on one half and scoring the other removes that bias; splitting by tissue block or
patient also stops near-identical neighbors from vouching for each other.

**Why compare with the original work last.** The authors' labels are an independent expert
annotation. Agreement shows my annotation is sound; reproducing their main findings shows it
is biologically meaningful; disagreements are where the interesting cells are. If KDE flags
concentrate in those disagreements, the verification is doing real work.

**Why two platforms.** Xenium gives sparse counts on hundreds of genes; CODEX gives dense
intensities on tens of proteins. The same workflow on both shows the approach is not tied to
one data type.

**What can go wrong and how we watch for it.**
- Over- or under-clustering → resolution chosen from marker coherence, with one sensitivity check.
- Names driven by expectation → every name backed by listed markers; unclear = `unknown`.
- KDE in too many dimensions → only the first k PCs; sensitivity to k reported.
- Bandwidth decides everything in a KDE → chosen by held-out likelihood, sensitivity shown.
- Different label granularity from the authors → explicit mapping table before comparing.

**Expected outputs per dataset.** Marker dot plot, annotated UMAP, tissue maps, neighborhood
enrichment heatmap, KDE flag map and flag rates per type, confusion matrix and composition
comparison against the authors, reproduction of the papers' main findings, and a short table
linking KDE flags to disagreements.
