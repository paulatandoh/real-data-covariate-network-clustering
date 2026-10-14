# Publish the real-data analysis repository

Suggested repository name: `covariate-network-clustering-real-data`.

Suggested description: `R code and prepared datasets for four real-data analyses in Model-Based Clustering of Network Data with Binary Covariate Information.`

1. Extract the GitHub-ready ZIP on your computer.
2. Create a repository at https://github.com/new. Choose the desired visibility. Do not generate another README or .gitignore; both are included.
3. Open the empty repository's upload-files link, or choose `Add file` > `Upload files`.
4. Open the extracted `Four_Real_Data_Analyses` folder. Drag its CONTENTS into the upload box, preserving folders. The four numbered R scripts, README.md and Rproj belong at the repository root, beside R/, data/, prepare/, reference/, docs/ and tests/.
5. Include the hidden .gitignore file and .github/ folder. Confirm .github/workflows/r-checks.yml is present. The package has fewer than 100 files and can be uploaded in one batch; folders must be dragged into the upload box rather than selected as ordinary files in a file-only picker.
6. Commit with message `Add four real-data analyses and reproducibility documentation`.
7. Open Actions to inspect the R verification workflow. It checks mathematical components, checkpoint recovery, saved partitions and smaller-technology input reconstruction. It does not fit all expensive observed-network models.
8. Keep generated outputs/fits/ locally. To publish final new result tables, add the relevant report CSVs, all-start diagnostics and session information in a clearly labeled results folder or versioned release. Preserve the distinction between historical reference files and new reruns.

Expected top-level files: 01_technology.R, 02_health.R, 03_insurance.R, 04_technology_750.R, README.md, Four_Real_Data_Analyses.Rproj and .gitignore. Expected folders: R, data, prepare, reference, docs, tests and .github.

Official GitHub instructions:
- https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository
- https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository

This project contains the four real-data analyses. The simulation code can remain in its separately prepared repository.
