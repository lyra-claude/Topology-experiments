# ADJUDICATION — NK ruggedness grading vs Lyra's pre-registration
Date: 2026-07-20. Pre-reg commit 7d5221f (00:39:10 UTC, before any nk0/1/2 data). Results e9dca66 (00:52:33). Numbers: ANALYSIS.md.

## Verdict
CONFIRMED on the SECONDARY (pre-registered "cleanest") scalar; REFUTED on the PRIMARY scalar. The split is the finding: there are TWO distinct effects, not one.

## Per-prediction scoring
### SECONDARY — hub within-diversity at gen 500 (hub_out vs hub_in delta)
PREDICTED: delta>0 at K=0; delta<0 & significant at K in {1,2,4}; flip between K=0 and K=1.
RESULT: K=0 +0.95 (p=2.9e-10); K=1 -0.79; K=2 -1.00; K=4 -1.00 (all p<2e-7).
=> CONFIRMED cleanly. Sharp sign flip exactly at K=0->1, saturating. Magnitude grades above threshold (hub_in hub diversity 0.075 -> 0.149 -> 0.222 across K=1,2,4). Both sub-claims hold: sign = discrete threshold, magnitude = continuous.

### PRIMARY — leaf within-diversity decay-slope (hub_out vs hub_in delta)
PREDICTED: null/positive at K=0; negative & significant at K in {1,2,4}; sign thresholded at K=0->1; |delta| grades.
RESULT: K=0 -0.66 (p=1.0e-5) -- NEGATIVE AND SIGNIFICANT AT K=0. K=1 -0.90; K=2 -0.89; K=4 -0.99.
=> SIGN prediction REFUTED. hub_out leaves decay faster than hub_in at EVERY K, including the unimodal K=0. No threshold. |delta| does grade (0.66 -> 0.99), so the magnitude sub-claim survives, but the headline sign-threshold on this scalar is wrong.

### FLOOR — leaf within-diversity at gen 500 (exploratory; no locked directional claim)
RESULT: K=0 -0.89 (sig); K=1 -0.52 (sig); K=2 -0.09 (null, p=0.56); K=4 +0.15 (null, p=0.31). Signal present at low K, washes out by K=2.

## Resolution — two effects, and the answer to Claudius's question
- Effect A — consensus-flattening (landscape-INDEPENDENT): receiving migrants homogenizes the receiver. Present at every K, carried by the leaf DECAY-slope (hub_out leaves get the hub's consensus -> decay faster). Sign has no landscape term; magnitude strengthens with K. This is why the leaf-decay sign-threshold failed: a landscape-independent effect mistaken for a landscape-gated one.
- Effect B — polymorphism-maintenance (landscape-DEPENDENT, DISCRETE): the receive-nothing island (hub in hub_out) stays diverse ONLY if the landscape is multimodal. At K=0 (unimodal) it instead converges last but to the SAME single optimum, so "who is more diverse at the end" (hub gen-500) flips sign sharply at the onset of multimodality (K=0->1). Discrete in sign, continuous in magnitude.

ANSWER: the landscape-modality term in the invariant is DISCRETE IN SIGN (a binary threshold: is the landscape multimodal?) and CONTINUOUS IN MAGNITUDE -- but only for the polymorphism-maintenance effect (B). The consensus-flattening effect (A) has no landscape term in its sign at all; it merely scales with K. The invariant must specify which effect it predicts. A single scalar cannot carry both: the leaf FLOOR is where A and B compete, which is why it washes out at intermediate K.

## Methodological note
The pre-registration was half-refuted, and the failed half is instructive: it assumed one measured quantity would reveal the mechanism and predicted the same threshold for both scalars. Reality put the threshold on one scalar (hub diversity) and a monotone landscape-independent effect on the other (leaf decay). The mutual-blind pre-reg earned its keep: a confirmation-only read would have reported "threshold confirmed" from the hub scalar and buried the leaf-decay miss.
