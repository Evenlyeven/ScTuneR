<p align="center">
  <img src="ScTuneR_hex.png" width="180"/>
</p>

# ScTuneR: SCTransform, Integration, and PC × Resolution Tuning for Seurat

**ScTuneR** is a command-line R pipeline that takes single-cell samples from raw counts to a clustering-parameter survey. It normalizes each sample with `SCTransform`, merges them, corrects batch effects, and then sweeps a grid of **principal-component counts × clustering resolutions**, rendering every combination as a labelled UMAP panel so the right parameters can be chosen by eye.

It supports two integration methods:

- **`harmony`** — Harmony correction of the PCA embedding, written to reduction `harmony`
- **`rpca`** — Seurat reciprocal-PCA anchor integration of the same embedding, written to reduction `integrated.rpca`

Integration can also be skipped entirely (`--integration FALSE`), in which case the sweep runs on the uncorrected `pca` reduction.

The pipeline runs in four stages, each writing a checkpoint that `--resume_from` can re-enter:

	1.	**Load, split, normalize** — SCTransform per sample → `SCTed_seu_obj_ls.rds`
	2.	**Merge + PCA** — shared variable features, `RunPCA` → `merged_seu_obj.rds`
	3.	**Integration** — Harmony or RPCA → `integrated_seu_obj.rds`
	4.	**Exploration** — the PC × resolution UMAP grid → `UMAP_PCres.png`

## 📦 Features

	•	Accepts a single Seurat object (split by a metadata column) **or** a named list of Seurat objects, in `.rds` or `.RData` format
	•	Per-sample `SCTransform` (`vst.flavor = "v2"`, 3000 variable features)
	•	Optional cell-cycle regression (`S.Score`, `G2M.Score`) during normalization
	•	Two integration methods: `harmony` and `rpca` — or no integration at all
	•	Variable features chosen by cross-sample recurrence rather than union, so the PCA is driven by genes variable *across* samples
	•	`PrepSCTFindMarkers()` run on the merged object, ready for downstream DE
	•	Grid sweep over any comma-separated `--pc_range` × `--res_range`
	•	Every UMAP panel captioned with its PC, resolution, and source reduction
	•	Elbow plot of the merged PCA for sanity-checking the PC range
	•	Grid validated up front, so bad parameters fail in seconds — not after hours of normalization
	•	`--resume_from` re-enters the pipeline at any checkpoint after a failure
	•	Logs the reduction name the follow-up **ScHerdeR** run needs
	•	Timestamped output directory (prevents overwriting)

## 🚀 Quick Start

### 🔧 Requirements

Install required R packages:

```r
install.packages(c(
  "tidyverse",
  "Seurat",
  "optparse",
  "magrittr",
  "patchwork"
))
```

Integration additionally needs the backend for the chosen method:

```r
install.packages("harmony")   # for --integration_method harmony
```

### 🖥️ Basic Usage

**Default run** — split one object by sample, SCTransform, Harmony, then sweep:
```bash
Rscript ScTuneR.R \
  --seurat_obj path/to/seurat_object.rds \
  --output_dir results/ \
  --pc_range "20,30,40" \
  --res_range "0.2,0.6"
```

**Pre-split input** — a named list of Seurat objects, one per sample:
```bash
Rscript ScTuneR.R \
  --seurat_obj path/to/seurat_list.rds \
  --output_dir results/ \
  --sub_analysis FALSE \
  --pc_range "30,40,50"
```

**RPCA integration** instead of Harmony:
```bash
Rscript ScTuneR.R \
  --seurat_obj path/to/seurat_object.rds \
  --output_dir results/ \
  --integration_method rpca \
  --pc_range "20,30"
```

**Wider PC range** — compute more PCs so the sweep can reach them:
```bash
Rscript ScTuneR.R \
  --seurat_obj path/to/seurat_object.rds \
  --output_dir results/ \
  --npcs 100 \
  --pc_range "50,75,100" \
  --res_range "0.4,0.8"
```

**Regress out cell cycle** (human only):
```bash
Rscript ScTuneR.R \
  --seurat_obj path/to/seurat_object.rds \
  --output_dir results/ \
  --regress_cc TRUE
```

**Resume from a checkpoint** — re-run only the sweep on an already-integrated object:
```bash
Rscript ScTuneR.R \
  --seurat_obj results/ScTuneR_20260807_141500/integrated_seu_obj.rds \
  --output_dir results/ \
  --resume_from integrated \
  --integration_method harmony \
  --pc_range "30,40" \
  --res_range "0.2,0.4,0.6"
```

### 📝 Parameters

**Required**
| Parameter | Description |
|-----------|-------------|
| `--seurat_obj` | Path to a Seurat object or a named list of Seurat objects (`.rds`, `.RData`). When `--resume_from` is set, this must point at the matching checkpoint file |

**Input Handling**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--output_dir` | `.` | Directory where the timestamped output folder is created |
| `--sub_analysis` | `TRUE` | `TRUE`: input is one Seurat object, split by `--splitby`. `FALSE`: input is already a named list of Seurat objects |
| `--splitby` | `orig.ident` | Metadata column used to split the object into per-sample pieces |

**Normalization**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--regress_cc` | `FALSE` | Regress out `S.Score` and `G2M.Score` during SCTransform. **Human only** — uses Seurat's built-in `cc.genes` |

**Integration**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--integration` | `TRUE` | Run batch integration. `FALSE` sweeps the uncorrected `pca` reduction |
| `--integration_method` | `harmony` | `harmony` (writes reduction `harmony`) or `rpca` (Seurat reciprocal PCA, writes `integrated.rpca`) |

**Exploration Grid**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--npcs` | `50` | Number of PCs computed by `RunPCA`. Must be ≥ the largest value in `--pc_range` |
| `--pc_range` | `30,40,50` | Comma-separated PC counts to explore (each ≥ 2) |
| `--res_range` | `0.2,0.4,0.6` | Comma-separated clustering resolutions to explore (each > 0) |

**Resume**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--resume_from` | `NULL` | Re-enter the pipeline at a checkpoint instead of starting from raw input. One of `sct` (`--seurat_obj` = `SCTed_seu_obj_ls.rds`), `merged` (= `merged_seu_obj.rds`), or `integrated` (= `integrated_seu_obj.rds`) |

## 📂 Output Structure

Each run creates a timestamped directory: `ScTuneR_20260807_141500/`
<br>
**Stage 1 — normalization**
	•	SCTed_seu_obj_ls.rds — named list of per-sample SCTransformed objects
  <br>
**Stage 2 — merge + PCA**
	•	merged_seu_obj.rds — merged object, `PrepSCTFindMarkers()` and `RunPCA()` done
	•	merged_PCA.png — elbow plot of the computed PCs
  <br>
**Stage 3 — integration**
	•	integrated_seu_obj.rds — object carrying the corrected reduction (skipped when `--integration FALSE`)
  <br>
**Stage 4 — exploration**
	•	UMAP_PCres.png — UMAP grid, one panel per PC × resolution pair, laid out with resolutions across columns and PCs down rows

## 📌 Notes

 - Both integration methods take their batch grouping from the separately SCTransformed samples — i.e. the split set by `--splitby`, or the names of the input list. **Neither reads a metadata column**, and the Harmony batch variable is not settable. Each run logs the groups it actually used.
 - Harmony corrects *every* dimension of the PCA, so `--npcs` sets how many dims it sees. Raise it only when a larger `--pc_range` genuinely needs it — larger values cost time and memory.
 - RPCA finds anchors in dims 1:30 (`RPCAIntegration`'s default) but corrects the full embedding, so `--pc_range` is not capped at 30.
 - The variable features passed to RPCA are handed over explicitly. Setting `VariableFeatures()` alone is not enough: for an `SCTAssay`, `IntegrateLayers` would otherwise overwrite them via `SelectSCTIntegrationFeatures()`.
 - `--pc_range`, `--res_range` and `--npcs` are validated before any heavy computation, so an inconsistent grid fails immediately rather than after normalization and multi-GB writes.
 - When resuming from `integrated`, pass the same `--integration_method` that produced the object — otherwise the expected reduction will not be found.
 - The list form of the input must be **named** (e.g. by sample ID). If several Seurat objects or lists are present in an `.RData` file, the first one is used.
 - Sweep runtime scales as `length(pc_range) × length(res_range)`, each combination running a full `FindNeighbors` → `FindClusters` → `RunUMAP` pass. Clustering uses a fixed seed (`10086`), so panels are reproducible.
 - Feed the chosen PC, resolution, and the logged reduction name into **ScHerdeR** (`--npc`, `--resolution`, `--reduction`) to produce the final clustering.
 - Each run creates a new timestamped output folder to prevent overwriting.
