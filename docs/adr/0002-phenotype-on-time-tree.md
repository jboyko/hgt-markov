# 0002 — Simulate phenotype directly on the time tree (single clock)

Date: 2026-05-26
Status: Accepted

## Context

In a typical phylogenetic comparative methods (PCM) study, branch lengths come in *substitution units* (estimated from sequence data with possibly heterogeneous rates across lineages), and the discrete-trait model Q is then defined per substitution unit. There are two clocks in play: a real-time clock for the underlying lineage history and a substitution-unit clock for what inference sees.

The empirical use case here is bacterial **phenotype** evolution — antibiotic resistance carried on plasmids, for example — not sequence evolution. The phenotype changes in real time, not per substitution. And HGT is fundamentally a real-time event: two lineages must be alive at the same calendar moment to exchange a plasmid.

We therefore needed to choose: keep the conventional two-clock setup (time tree for HGT exposure, substitution tree for inference), or collapse to a single time clock.

## Decision

**A single time clock is used throughout.** Trees are simulated under a birth-death-sampling process (TreeSim) directly in time units. Branch lengths handed to corHMM for inference are those same time-unit branch lengths. Q is parameterised in events per unit time. λ is in events per lineage per unit time. Non-ultrametricity arises from **heterochronous tip sampling** (bacterial isolates from different collection dates), not from rate heterogeneity across lineages.

## Consequences

- **λ, Q's off-diagonals, and tree branch lengths are all in the same units**, making λ/q̄ a meaningful, dimensionless dial. Results can be reported as "HGT rate equal to N× the mean substitution rate."
- **"Co-existing at time t" is unambiguous**: a branch is alive between its parent node's time and its own child's time (or sampling time). Defining the donor-eligibility set requires no extra convention.
- **Two confounds are deliberately separated.** The "HGT broke Mk" signal is not entangled with "rate heterogeneity across lineages broke Mk." Rate heterogeneity is only present in the HMM-null generator (ADR 0003), where it is the *truth* being recovered.
- The simulator does not need a second pass to map time to substitutions, and inference does not need a time-to-substitution rescaling.
- **Trade-off accepted**: results are not directly comparable to PCM studies that operate on substitution-length trees with possibly heterogeneous molecular clocks. A follow-up could add a substitution-clock cell, but conflating molecular-rate variation with HGT damage was judged to muddy the headline question.
