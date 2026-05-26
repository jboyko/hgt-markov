# 0003 — Inference suite spans HMM models to test whether HGT masquerades as hidden rate heterogeneity

Date: 2026-05-26
Status: Accepted

## Context

A standard robustness study fits the *matched* model (here, Mk-ARD) to data simulated with a violation, and reports parameter bias. That answers "how wrong are the estimates?" but stops there.

corHMM users routinely fit hidden-rate variants (HMM with 2, 3, ... rate categories) and use AIC/AICc to decide whether genuine hidden biological heterogeneity is present. If HGT-perturbed data systematically *looks* like hidden rate heterogeneity to model selection — i.e., AICc spuriously prefers an HMM when truth is Mk-ARD + HGT — then practitioners are being fooled into "discovering" biology that isn't there. The project lead develops corHMM and is uniquely positioned to say this.

The alternative was to keep the inference side simple (matched ARD only), saving compute and keeping the writeup narrow.

## Decision

Every simulated dataset is fit with a **fixed inference suite** of corHMM models spanning vertical complexity: ER, ARD, HMM-ARD with 2 hidden rate categories, and HMM-ARD with 3 hidden rate categories. The suite is identical across all three generators (Mk-null, Mk+HGT, HMM-null), so each generator's behaviour can be plotted on the same axes.

The headline diagnostic is **P(HMM-preferred by AICc) as a function of λ** under the Mk+HGT generator, compared against (i) Mk-null as the false-positive floor and (ii) HMM-null as the genuine-heterogeneity ceiling.

## Consequences

- **The paper's framing shifts** from "HGT biases Mk rate estimates" (known, mild) to "HGT looks like hidden rate heterogeneity to model selection" (sharper, actionable for corHMM users).
- Compute cost roughly quadruples vs. matched-only fitting (4 models × every dataset). Tractable at the planned scale (~3000 datasets, ~12k fits) given parallelisation and `use_RTMB`.
- The suite must be fixed *before* the sweep starts. Adding or changing models mid-sweep invalidates cross-cell comparisons of P(model-X-wins).
- **Explicitly out of scope**: an HGT-aware fitter (e.g. a mixture model or instantaneous-jump term in Q). It would be the "correct" model but requires bespoke implementation outside corHMM and was judged too costly for this study.
- **Tied to ADR 0002**: because no substitution clock is in play, any rate heterogeneity detected by an HMM fit cannot be explained away as molecular-rate variation across lineages — it really is "the discrete-character process looking heterogeneous."
