# Example data

A synthetic dataset small enough to ship, so you can confirm your install works
and see real output before you have your own mass-spec export.

Everything here is **fabricated**. The proteins do not exist and the spectra were
never measured. The point is that the right answer is known in advance, which
makes this a positive control as well as a smoke test.

| file | what it is |
|---|---|
| `example.mzid` | an mzIdentML export with three runs — `sample_01`, `sample_02` (treatment) and `sample_03` (control) |
| `structures/TOYA.pdb`, `structures/TOYB.pdb` | matching structures, so the run needs no internet |
| `peg8.json` | a PEG8 scaffold definition, the per-bond form filled in |

## Run it

From the repository root:

```bash
covalent-pocket-mapper \
  --mzid examples/example.mzid \
  --polymer examples/peg8.json \
  --structures examples/structures --no-download \
  --treatment sample_01 sample_02 \
  --control sample_03 \
  --output /tmp/cpm-example \
  --permutations 200
```

It finishes in a few seconds and needs no network.

## What you should see

Two proteins, built to give opposite answers:

```
 gene accession  n_signal_sites  credible_region_A   z_score  p_value  control_fraction_reachable
TOYA1      TOYA               5               15.7  3.665824 0.004975                         0.0
TOYB1      TOYB               5                4.5 -1.493984 0.970149                         NaN
```

**TOYA is the positive control.** Its five capture sites were placed around one
surface pocket, within reach of a PEG8 scaffold. The tool finds a pocket at
`p = 0.005`, and none of its three control-only sites are reachable from that
pocket — `control_fraction_reachable` is 0.0, which is what specificity looks
like.

The pocket was planted at `(30, 0, 0)` and recovered at about `(22, -1, -1)`,
roughly 8 Å away. That is not an error to worry about: the credible region is
16 Å wide, so the true pocket sits well inside it. Five sites on one protein
localise a pocket to a neighbourhood, not to a point, and the credible region
is the column that tells you so.

**TOYB is the negative control.** Its five sites were spread as far apart on the
surface as possible, so no single anchor point explains them. It comes back at
`p = 0.97` — correctly unconvincing. Note its credible region is *narrower* than
TOYA's, which is the trap the main README warns about: a tight region means few
locations are feasible, not that the evidence is good. Read `p_value` first.

## Regenerating

These files were generated, not hand-written. If you want to change the
geometry, the site layout, or the number of proteins, the parameters that matter
are the globule radius, the number of residues, and how far the capture sites
sit from the planted pocket. A protein much larger than the scaffold's reach is
what makes the test discriminating — when the scaffold can span the whole
protein, a random scatter of sites fits about as well as a real cluster, and
`z_score` collapses toward zero.
