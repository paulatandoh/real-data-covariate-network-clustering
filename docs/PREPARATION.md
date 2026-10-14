# Optional preparation and provenance

The default analysis needs no preparation: it uses `data/<dataset>/*csv` supplied with the complete export. Preserve their hashes. Optional preparation should be done in a separate working directory, then compare the resulting graph and ordered covariate against the supplied target before substituting inputs.

## Smaller GitHub technology

`prepare/technology_from_saved_nodes.R` reconstructs the induced graph on the exact 300 original GitHub IDs in the saved covariate file. It reads the supplied raw SNAP ZIP and verifies every entry of the adjacency matrix. It does not invent or reconstruct the unrecorded historical rule used to choose those 300 IDs. Run it from the package root; it only verifies/recreates in a temporary workspace.

## Health

`prepare/prepare_addhealth_data_v2.R` is the latest supplied preparation code. It loads `amen::addhealthc9`, removes its single all-missing attribute row, removes missing sex, binarizes and OR-symmetrizes nominations, removes loops/isolates and produces local sequential IDs. Install `amen` before running in a separate workspace. The prepared target is N=247, E=999, male/female counts 122/125. Validate node alignment against the precise original packaged data; equal dimensions alone are insufficient.

## Insurance

`prepare/prepare_insurance_dukou_data.R` uses `haven` and `igraph` and expects `0422allinforawnet.dta` in its working directory. Copy the raw file from `data/raw/` to a separate preparation workspace. The packaged copy corrects the stale expected edge count from 671 to 667. Both endpoints are restricted to Dukou; placeholder contact 99 and self-ties are removed; observed take-up is mapped by contact/network ID after consistency checks; only verified binary-covariate nodes and their observed edges are retained. The resulting 250-node, 667-edge set was independently reconstructed in Python and equals the supplied adjacency. This verifies implementation consistency; attribute ownership should still be justified from the source/codebook in the paper.

## Label-aware GitHub N750

`prepare/prepare_large_github_balanced_N750.R` is the original preparation code. Copy `data/raw/git_web_ml.zip` to a separate working directory first. It chooses a high-degree cross-type starting edge, grows a connected degree-prioritized frontier with label quotas 500 Web/250 ML, and retains all original induced edges. Tie seed is 20260809. Graph ordering and tie behavior can depend on software version; the supplied prepared input and original-ID mapping are the reproducibility target. The file's `balanced` label refers to a development filename, not 50/50 covariate proportions. The other BFS 750-node development graph (1084 edges) is not used for this manuscript analysis.
