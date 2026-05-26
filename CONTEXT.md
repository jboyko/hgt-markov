# CONTEXT

Glossary of domain terms used in this project. Implementation details live elsewhere.

## Project framing

This project tests the **robustness of the standard Mk model to violations of the vertical-inheritance assumption** introduced by **HGT**. The empirical use case is bacterial phenotypes carried on mobile elements (e.g. plasmid-borne antibiotic resistance). No genetic sequence is simulated; characters are simulated *directly as phenotypes* on a time-calibrated phylogeny.

## Terms

### HGT (horizontal gene transfer)
A character-state copy event from one co-existing lineage (the **donor**) to another (the **recipient**), overwriting the recipient's current state. Distinct from a **regime shift**, which changes the *rate matrix* on a subclade rather than transferring state between contemporaneous lineages.

### Donor / Recipient
The two lineages involved in an HGT event. The donor's *current* state at the event time is copied into the recipient. Donor's own state is unchanged.

### Co-existing lineages
The set of branches alive at a given time *t* on the time-calibrated tree. A branch is alive between its parent node's time and its own child node's (or sampling) time. Eligibility for HGT pairing is defined strictly by this real-time co-existence, not by substitution-length proximity.

### λ (HGT rate)
The per-lineage rate of *receiving* an HGT event, in units of events per lineage per unit time. Directly commensurable with off-diagonal entries of Q. λ = 0 reduces the model to standard Mk.

### Silent transfer
An HGT event in which donor and recipient happen to share the same state, producing no change. Allowed and counted; do not reject and re-draw.

### Q
The instantaneous rate matrix of the Mk process governing vertical character change. Entries are in events per unit time. For the binary use case Q is 2×2 with two free rates (q01, q10) under ARD.

### Generator
One of the three data-generating models used in the sweep: **Mk-null** (Mk-ARD, λ=0), **Mk+HGT** (Mk-ARD, λ swept), and **HMM-null** (true hidden-rate model, λ=0). The generator is the "truth" that inference is judged against.

### Inference suite
The set of corHMM models fit to every simulated dataset, spanning vertical complexity from ER through ARD to HMM with 2 and 3 hidden rate categories. The same suite is fit regardless of which generator produced the data.

### Heterochronous tip
A tip sampled before the most recent time in the tree — i.e. a non-extant lineage. Bacterial isolates from different collection dates produce heterochronous tips and are the source of non-ultrametricity here.

### Cell
One configuration in the sweep: a (generator, λ) pair (plus any other swept axes added later). Each cell holds many independent replicates (each its own tree, character history, and fits).

### Replicate
One end-to-end run: simulate a fresh tree under the birth-death-sampling process, simulate a character history under the generator, fit the full inference suite, record metrics.

### Brier score
Proper scoring rule used to summarise ancestral-state reconstruction (ASR) error: mean squared difference between corHMM's marginal posterior probability over states at each internal node and the indicator of the simulated true state.

### Model-selection misbehavior
The event that the inference suite's AICc-winner is *not* the structurally matched model — most notably, an HMM model winning when the generator is Mk+HGT. The central diagnostic of the study: how often does HGT make a vertical-only fitter "discover" hidden rate heterogeneity that isn't biologically there?
