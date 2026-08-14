# FedIRT Gradient-Space Anomaly Screening (Proof of Concept)

A proof-of-concept for detecting contaminated student response rows in a 2PL IRT
model within FedIRT's federated EM algorithm. Contaminated rows are flagged using unsupervised anomaly
detection on **per-student gradient vectors**, rather than on raw responses.

## Motivation

FedIRT (Zhou, Luo & Ji) aggregates per-student gradients across sites into a shared
item-parameter updates. Its FedIRT-DP extension bounds the influence of extreme rows by
clipping every student's gradient, which absorbs contamination without identifying
faulty rows. This repo asks a complementary question: at a single site, before
aggregation, can we **detect** which rows are contaminated using only the gradient
vectors already computed in the E-step with no labels, no raw-response rules?

## Approach

- Generate clean 2PL data (single school, following the Study 1 / Study 3 DGP).
- Inject contamination with known ground-truth labels.
- Compute each student's contribution to the item-parameter gradient (a `2J`-length
  vector), evaluated at the true item parameters — the same per-student quantity
  FedIRT-DP clips, kept instead of summed away.
- Score rows with LOF / ABOD on those gradient features.
- Evaluate detection with AUC against the known labels, benchmarked against two
  baselines: LOF on raw responses, and a trivial row-sum rule.

## Key result

Detection difficulty depends entirely on the contamination type:

| Scenario | Rowsum baseline | Gradient-space LOF | Notes |
|---|---|---|---|
| All-ones (toy) | ~0.97 | ~.04 | Any naive rule catches it, but inverse relation on LOF |
| Partial straightlining | ~0.75 | ~0.42 | Mild/ambiguous, no method is clearly stronger|
| **Random responding** | **~0.38 (chance)** | **~0.76** |  |

On **random responding**, where the naive row-sum rule collapses to chance (a random responder's total looks normal),
label-free detection on gradient features recovers contaminated rows at ~0.76 AUC (mean over 20 seeds, SD ~0.04).
LOF and ABOD perform equivalently, suggesting features could be more significant than detector choice.

## Future Directions
- Conducting formal sweep across range of 0.1-0.5 of contaminated row proportions
- Examining relation of effective detection and number of E/M stages
- Other methods of anomaly detection (One-class SVM, Isolation Forest)