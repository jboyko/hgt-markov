# 0001 — HGT modelled as per-lineage Poisson donor→recipient state copy

Date: 2026-05-26
Status: Accepted

## Context

The project tests whether standard Mk inference is robust to violations of the vertical-inheritance assumption introduced by horizontal gene transfer. This requires a simulator that *introduces* HGT as a tunable violation. Several mechanisms are possible:

- **State copy between contemporaneous lineages** (true HGT): pick a donor and recipient alive at the same time, overwrite recipient's state with donor's.
- **Regime shift on a subclade**: change Q (the rate matrix) on a chosen subclade. The existing `code/sim-mk.r` carries commented-out scaffolding (`Q2`, `NoI`, `NewQDesc`) suggesting this was an earlier direction.
- **Recipient-only "jump" event**: recipient resets to a random state from some distribution; no donor lineage involved.

Within "state copy between contemporaneous lineages," the rate can be parameterised per-pair (λ·C(k,2)), per-lineage (λ·k), or per-tree (λ·tree_length).

## Decision

HGT is modelled as a **per-lineage Poisson process of donor→recipient state-copy events**. Each lineage alive at time *t* independently experiences "receive an HGT event" at rate **λ**. When such an event fires, the donor is drawn uniformly from the other co-existing lineages and the donor's current state is copied into the recipient (overwriting). "Silent" transfers (donor and recipient already in the same state) are allowed and counted as events, not redrawn.

## Consequences

- **λ is directly commensurable with the off-diagonals of Q** (both in events per lineage per unit time), so results can be reported as λ/q̄ ratios — the natural axis for "how much HGT relative to vertical substitution breaks Mk."
- λ = 0 reduces the model exactly to standard Mk, giving a clean null.
- The simulator cannot be written edge-by-edge in postorder (as `sim-mk.r` does); HGT couples co-existing lineages and requires a forward time-sweep over epochs between node events.
- **Regime-shift dynamics are out of scope** for the primary sweep. They were considered as an "innocent violator" comparator cell but deferred.
- **Recipient-only jumps were rejected** because they confound HGT with rate inflation: without a donor lineage, the event is just extra noise rather than the specific violation of vertical inheritance under study.
- Per-pair (λ·k²) and per-tree (λ·tree_length) parameterisations were rejected because their effective HGT pressure depends on tree size/shape, complicating cross-replicate comparison.
