# Manuscript result conventions

Use the selected start with minimum VBIC among finite outer-converged fits. Keep selected-fit criteria, parameters and partition results from the same run. The completed newer reruns are documented in `REPRODUCTION_STATUS.md`; historical baseline files are documented in `reference/README.md`.

The observed dataset counts are smaller technology N=300/E=1458, health N=247/E=999, insurance N=250/E=667, and label-aware technology N=750/E=15335. The insurance edge set reconstructed from the supplied Stata file exactly matches its prepared adjacency; the original preparation expectation of 671 was outdated and was corrected to 667 without altering graph construction rules.

For the three fresh reruns supplied by the author, smaller-technology beta is 0.8434849; insurance network-only mean stability ARI is 0.1135449 (SD 0.1133742), while insurance between-model ARI is 0.1206003. Use distinct labels for those two ARI measures. Health between-model ARI is 0.8137504 and the smaller technology value is 0.906599.

The historical N750 selected partitions reproduce between-model ARI 0.766487175960483 and RI 0.908485981308411. The corresponding selected covariate fit has VBIC 105108.435174673, VICL 105221.392896658 and beta 1.30906897919145. A final fresh N750 rerun remains pending.

The historical quick-report mode uses original saved outputs and may leave absent selected-run parameters NA. Fresh fit mode saves a full fit object for every start. Ordinary/run-1 summaries must not supply parameters for a different selected start. Runtime measurements should come from a consistently documented execution.
