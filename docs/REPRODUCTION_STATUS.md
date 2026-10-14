# Recorded reproduction status

## Fresh reruns reported by the author

The following values were transcribed from the author's completed Posit Cloud Console output on October 13, 2026. They have the precision printed in the Console. Complete new start-level CSVs, fit objects and session information remain in the author's Posit project and should be retained with final reported results. These records are separate from the older saved baseline files in `reference/original_outputs/`.

| Dataset | Model | Selected run | Seed | VBIC | VICL | Beta | Stability ARI | Between-model ARI |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Smaller technology | Network-only | 1 | 1001 | 9280.577 | 9284.591 | NA | 1.000000 | 0.906599 |
| Smaller technology | Covariate | 2 | 1002 | 9175.154 | 9179.427 | 0.8434849 | 1.000000 | 0.906599 |
| Health | Network-only | 2 | 2002 | 8779.094 | 8972.496 | NA | 0.6130951 | 0.8137504 |
| Health | Covariate | 2 | 1002 | 8776.264 | 8977.489 | 0.2322542 | 0.6374190 | 0.8137504 |
| Insurance | Network-only | 2 | 2002 | 6488.702 | 6811.997 | NA | 0.1135449 | 0.1206003 |
| Insurance | Covariate | 1 | 1001 | 6489.560 | 6815.479 | 0.2369646 | 0.2370809 | 0.1206003 |

Both selected models in each of these three reruns report `EM_converged=TRUE`, `All_inner_converged=TRUE` and `Numerical_failure=FALSE`. These selected-fit flags do not by themselves establish convergence for every unselected start.

## Distinguish stability from agreement

For insurance, the latest network-only mean stability ARI is 0.1135449 (SD 0.1133742), compared with 0.120380276874710 in the archived baseline. The between-model ARI is still 0.1206003. These are different statistics. The observed change in across-start stability has not been attributed to a specific start; the latest all-start diagnostics/pairwise CSVs are needed for that comparison.

For an updated table based on the latest rerun, insurance network-only stability rounds to **0.114**, whereas between-model agreement rounds to **0.121**. Smaller technology beta rounds to **0.843**.

## Timing

The displayed selected-fit runtimes (seconds) are:

| Dataset | Network-only | Covariate |
|---|---:|---:|
| Smaller technology | 521.0987 | 409.2329 |
| Health | 332.7880 | 460.9968 |
| Insurance | 380.6297 | 501.5628 |

Runtimes depend on the execution environment and interruptions. Keep the timing values and conventions from a documented execution; do not combine values from different runs as if measured together.

## Technology N750

Historical evidence supports selected network-only run 3 (seed 2003), selected covariate run 2 (seed 1002), between-model ARI 0.766487175960483 and RI 0.908485981308411. Independently recomputed stability ARIs are 0.973602259951642 (network-only) and 0.995302409035804 (covariate). The final fresh N750 rerun has not yet been supplied and is not labeled complete.

## Original parameter and output records

Older small-technology covariate diagnostics omitted beta for selected run 2; historical report mode therefore leaves it NA. The fresh run above resolves the parameter value. The older baseline files are deliberately unchanged so their provenance and automated saved-result checks remain valid.
