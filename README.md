# Real-data analyses: model-based network clustering with binary covariates

**Manuscript:** Model-Based Clustering of Network Data with Binary Covariate Information  
**Authors:** Paul H. Atandoh and Kevin H. Lee

This repository contains four real-data analyses comparing network-only and covariate-informed latent-group models. Each analysis uses K=2, three initializations per model, selection of the finite outer-converged fit with minimum VBIC, initialization-stability summaries and between-model partition agreement.

## Analysis entry points

| Script | Dataset | Nodes | Edges | Observed binary covariate (0 / 1) |
|---|---|---:|---:|---|
| `01_technology.R` | GitHub technology | 300 | 1458 | Web / ML developer; 223 / 77 |
| `02_health.R` | Add Health | 247 | 999 | Male / female; 122 / 125 |
| `03_insurance.R` | Dukou insurance | 250 | 667 | Uninsured / insured; 136 / 114 |
| `04_technology_750.R` | Label-aware GitHub technology | 750 | 15335 | Web / ML developer; 500 / 250 |

Required prepared input CSVs are included under `data/`. `R/legacy/` contains the eight original supporting base fitters invoked by these four entry points; run the entry scripts rather than choosing between supporting files.

## Model and settings

For the covariate-informed model:

```text
logit P(Y_ij = 1 | z, x) = theta[z_i] + theta[z_j] + beta * I(x_i == x_j)
```

Edges are independent conditional on latent memberships and the observed covariate. The network-only fit omits beta. Both models analyze the same observed network.

Initialization uses `lda` MMSB with 1000 iterations and burn-in 20. The E-step uses the supplied quadratic MM surrogate and `quadprog`. The M-step uses the original full Newton updates with a full Hessian. Inner/outer tolerances are 1e-4 and each iteration limit is 500. Inner-convergence diagnostics are recorded separately from outer convergence.

| Dataset | Network-only seeds | Covariate seeds |
|---|---|---|
| Smaller technology | 1001, 1002, 1003 | 1001, 1002, 1003 |
| Health | 2001, 2002, 2003 | 1001, 1002, 1003 |
| Insurance | 2001, 2002, 2003 | 1001, 1002, 1003 |
| Technology N750 | 2001, 2002, 2003 | 1001, 1002, 1003 |

## Run in RStudio or Posit Cloud

### 1. Set the working directory before running any script

All paths in this project are relative to the repository root. **Your working directory must be the folder containing `Four_Real_Data_Analyses.Rproj`, the four numbered scripts, `R/` and `data/`.**

In RStudio or Posit Cloud, open `Four_Real_Data_Analyses.Rproj`. Opening this project normally sets the working directory automatically. Check it in the Console:

```r
getwd()
file.exists("R/real_data.R")  # must be TRUE
dir.exists("data")          # must be TRUE
```

If either check is FALSE, set the directory explicitly. For a ZIP extracted into the folder `Four_Real_Data_Analyses` in Posit Cloud, use:

```r
setwd("/cloud/project/Four_Real_Data_Analyses")
```

If you cloned the suggested GitHub repository into Posit Cloud, use its actual folder name, for example:

```r
setwd("/cloud/project/covariate-network-clustering-real-data")
```

If the repository files are directly under `/cloud/project`, use `setwd("/cloud/project")`. On your own computer, replace the path with the full path to your extracted or cloned repository. **Choose the path that matches your files; do not run all these examples in sequence.** The leading `/` in a Posit Cloud path is required.

You can also select the folder in RStudio's Files pane and choose **More → Set As Working Directory**. Before proceeding, run:

```r
stopifnot(
  file.exists("R/real_data.R"),
  dir.exists("data"),
  file.exists("01_technology.R"),
  file.exists("04_technology_750.R")
)
```

If this check fails, correct the folder path or upload the missing folders. The four numbered scripts alone are insufficient; keep the repository's supporting folders together.

### 2. Install packages and run an analysis

Install fitting packages once:

```r
install.packages(c("quadprog", "igraph", "lda", "fossil"))
```

Then source one entry script at a time:

```r
Sys.setenv(REAL_DATA_MODE = "fit")
source("01_technology.R")
```

After it finishes, run the remaining scripts individually:

```r
source("02_health.R")
source("03_insurance.R")
source("04_technology_750.R")
```

The equivalent direct call is `run_real_analysis("health", mode="fit")` after sourcing `R/real_data.R`. Default execution fits both models and saves complete fit objects. The N750 analysis can take several hours across its six starts.

## Outputs

Selected results are written to `outputs/fits/<dataset>/reports/`:

- `model_comparison.csv`: selected-fit parameters, criteria, stability and between-model agreement;
- `network_only_best_clusters.csv`, `covariate_best_clusters.csv`: IDs, hard labels, covariates and probabilities;
- `*_all_starts.csv`: seeds, objectives, iterations, numerical failures and inner/outer convergence;
- `*_pairwise_stability.csv`: pairwise ARIs across stored start partitions;
- `between_model_agreement.csv`: ARI and RI between the selected model partitions.

Each start saves its full parameter/probability object, ELBO trace, log and restart checkpoints. Run configurations and session information are recorded. `outputs/` is excluded from Git commits; retain it locally for restart and save the chosen final result files with each reported experiment.

Stability measures sensitivity to initialization. Between-model ARI measures agreement between two selected partitions, not accuracy against known ground truth. Lower VBIC/VICL indicates stronger support under the specified penalties; modularity is a complementary structural diagnostic.

## Resume after restart

Reopen the same project, confirm or reset the working directory as described above, and source the same entry script in fit mode. Keep `outputs/fits/`.

Completed starts reuse their saved fit. Interrupted starts resume at the latest committed outer EM iteration; unfinished E/inner steps are rerun from the preceding completed iteration. Initial MMSB setup is rerun before restoring the numerical state. Resumed timing excludes downtime and startup recomputation.

Checkpoint identities include input and fitting/shared-code hashes, seeds/settings, R version/platform/RNG, BLAS/LAPACK and fitting package versions. A change stops reuse. For a changed configuration use a new output root, for example:

```r
source("R/real_data.R")
run_real_analysis("health", mode="fit", output_root="outputs_new")
```

This publication update changes documentation/reference packaging only; the four main scripts and fitting code are unchanged from the package used for the author's completed reruns.

## Saved-reference reports

To reconstruct the supplied historical summaries quickly:

```r
Sys.setenv(REAL_DATA_MODE = "report")
source("01_technology.R")
```

This reads saved reference outputs and writes `outputs/historical/<dataset>/`; it does not rerun models. Unknown historical selected-run parameters remain NA rather than being filled from a different start. Later fresh rerun values are documented in `docs/REPRODUCTION_STATUS.md`.

## Reproduction status

The author completed fresh Posit Cloud reruns for the smaller technology, health and insurance networks. Their selected fits report outer and inner convergence with no numerical failure. The smaller technology selected beta is **0.8434849**. The insurance network-only mean stability ARI in the latest rerun is **0.1135449**; its between-model ARI remains **0.1206003**.

The N750 saved selected-fit criteria and both three-start stability calculations have been independently checked against the supplied labels. Its final fresh rerun is pending. Detailed recorded results and timing conventions are in `docs/REPRODUCTION_STATUS.md`.

## Preparation, checks and upload

Exact prepared graphs are the default inputs. Optional preparation code and raw source files are included; see `docs/PREPARATION.md`. The insurance raw-file reconstruction agrees with the supplied 250-node, 667-edge graph.

Run the base-R checks from the repository root:

```r
source("tests/test_model_math.R")
source("tests/test_pipeline.R")
source("tests/test_saved_labels.R")
```

GitHub Actions runs these checks plus the exact saved-node smaller-technology preparation check. These automated checks complement the observed-network reruns and do not run all six expensive fits for each dataset.

See `docs/GITHUB_UPLOAD.md` for publishing steps.
