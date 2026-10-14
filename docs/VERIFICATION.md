# Verification record

Completed under base R 4.3.3 on Linux:

- Parsed all four entry scripts, eight original supporting fitters, shared runtime/config and preparation files.
- Verified all four input matrices are binary, square, symmetric, loop-free, ID-aligned, and match documented node/edge/covariate counts.
- Reconstructed the smaller GitHub induced graph from exact saved original node IDs and supplied raw SNAP data; all edges and labels agree.
- Reconstructed insurance in pandas from the raw Stata file; its exact edge set and attribute lookup match the supplied 250-node, 667-edge prepared inputs.
- Reconstructed four historical selected-fit reports; selected VBICs and between-model ARIs agree with the saved diagnostics/partitions.
- Independently recomputed both large-network three-start stability calculations and best-partition between-model agreement from the saved label vectors.
- Tested all eight original base-model objectives, E-step coefficients and Newton gradient/Hessian assembly against independent direct/finite-difference calculations.
- Tested interruption after a committed outer iteration, recovery to the same final numerical state as an uninterrupted test-engine execution, completed-fit reuse, and rejection of changed checkpoint identities/seeds.

The checkpoint integration test uses a small deterministic test engine; it does not initialize MMSB or solve a real network fit. The author has supplied completed fresh Posit Cloud outputs for smaller technology, health and insurance. Both selected fits in each dataset report outer and inner convergence without numerical failure. The local automated checks and supplied historical-output reconstruction are separately verified. The final fresh N750 rerun is pending; see REPRODUCTION_STATUS.md. Original package/R versions were not supplied, so a fresh numerical reproduction must be checked against the historical reports after fitting in the chosen environment.
