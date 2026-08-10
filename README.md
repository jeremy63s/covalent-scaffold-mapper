# covalent-pocket-mapper

Finds candidate covalent binding pockets from scaffold pulldown mass spectrometry.

When a chemical scaffold with a reactive group is pulled down with its bound
proteins and the result is digested and sequenced, the peptides that come back
tell you *where* on each protein the scaffold attached. If a protein is a real
binder, those attachment points should not be scattered at random over its
surface — every one of them should be reachable from the single pocket the
scaffold was anchored to.

The scaffold is a polymer, so how far it reaches is not one number but a
distribution: P(r), the chance its two ends are a distance r apart, which
follows from its chemistry. Each capture site therefore casts a cloud of
probability into the space around it, and the tool multiplies those clouds
together. Where they agree is the pocket. Where any one of them is zero — past
the point the chain physically cannot span — the answer is zero no matter how
well the rest agree.

## What you need

one mzIdentML export from your search — `.mzid` or `.mzid.gz`.
In Scaffold this is `Export > mzIdentML`. That single file carries everything
the analysis uses: peptide sequences, the protein each maps to, the peptide's
start and stop position within that protein, and which run each spectrum came
from. No spreadsheet exports, peptide CSVs, or FASTA files are needed.

**Required:** your scaffold's chemistry, as a short form you fill in once. See
[Describing your scaffold](#describing-your-scaffold) — it is table lookups, not
a structure file.

Neither is needed just to try the tool: a synthetic export, matching structures,
and a filled-in PEG8 form ship in [examples/](examples/). See [Check it
works](#check-it-works).

**Optional:** a negative control run in the same export — a
scrambled scaffold, beads-only, or DMSO condition. See
[Running without a control](#running-without-a-control).

**Optional:** a directory of PDB structures named by accession (`P08670.pdb`).
Without one, structures are downloaded from AlphaFold and cached, so the second
run is fast. Long proteins are split across several AlphaFold fragments; these
are stitched into a single model, keeping the highest-confidence copy of each
residue. Alternative isoforms are ignored — only the canonical sequence is
numbered the way your peptide positions are.

## Install

Needs Python 3.9 or newer and nothing else — no compilers, no system libraries,
no external tools.

```bash
git clone https://github.com/jeremy63s/covalent-pocket-mapper
cd covalent-pocket-mapper

python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

pip install .
```

That installs the dependencies and puts a `covalent-pocket-mapper` command on your
path, which works from any directory. Check it came through:

```bash
covalent-pocket-mapper --help
```

If you plan to edit the code, use `pip install -e .` instead so your changes
take effect without reinstalling.

If you would rather not install anything, `pip install -r requirements.txt` and
then run `python -m covalent_pocket_mapper` from the repository root (the folder that
contains `covalent_pocket_mapper/`). Every command below works the same way, just with
that longer form.

### Check it works

A small synthetic dataset ships with the repository, so you can confirm the
install before you have your own data. From the repository root:

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

A few seconds, no network. You should get two proteins with deliberately
opposite answers — `TOYA1` at `p = 0.005`, a real pocket planted where the tool
finds it, and `TOYB1` at `p = 0.97`, sites scattered so that no single pocket
explains them. If you see those two lines, everything works. [examples/](examples/)
explains how the data were built and what the numbers mean.

## Run

The workflow is five steps: list the samples, do a small trial run, run the full
set, read the results table, then open a plot to inspect a candidate. The table
is the answer; the plots are for checking it.

With no arguments, `covalent-pocket-mapper` asks for everything one question at a
time, including the scaffold form, which it explains as it goes — a fine way to
start. The steps below are the same thing, run explicitly.

### 1. List the samples

Your treatment and control names live *inside* the export, not in the filename.
List them before choosing:

```bash
covalent-pocket-mapper --mzid ~/data/experiment.mzid.gz --list-samples
```

Pick which runs are treatment and, if you have one, which is the negative
control (a scrambled scaffold, beads-only, or DMSO run).

### 2. Trial run

Structures download on demand, so the first full run on a large export is slow.
Confirm the data parses and your scaffold settings behave on a handful of
proteins first:

```bash
covalent-pocket-mapper --mzid ~/data/experiment.mzid.gz \
  --treatment run_01 run_02 \
  --control run_03 \
  --output ~/results/trial \
  --polymer peg8.json \
  --max-proteins 5 --permutations 50
```

`--max-proteins` scores only the first few proteins that clear the site
threshold — enough to prove the pipeline runs, not a real result.

### 3. Full run

Drop the trial caps and point at a fresh output directory. Leaving
`--max-proteins` off scores everything:

```bash
covalent-pocket-mapper --mzid ~/data/experiment.mzid.gz \
  --treatment run_01 run_02 \
  --control run_03 \
  --output ~/results/experiment \
  --polymer peg8.json
```

`--polymer` is the saved per-bond form describing your scaffold's chemistry; the
guided run offers to write one for you, and [Describing your
scaffold](#describing-your-scaffold) explains what goes in it. If you took the
persistence-length shortcut instead, swap it for
`--persistence-length 3.8 --scaffold-size 28.9`.

### 4. Read the results

`binding_site_scores.csv` is the answer, ranked most convincing first. Read
`p_value` first: it is the evidence, and a real binder shows a low one. Then
check `credible_region_A` to see how tightly the pocket is actually pinned down,
and `null_feasible_fraction` — when few random scatters admit any pocket at all
but yours does, that is the strongest signal the tool produces.

Proteins with an empty `binding_site_x` had no surface point reachable from
every capture site. That is a real answer, not a failure: no single pocket
explains them.

This is where you find candidates — not in the plots, which are drawn only for
the top-ranked proteins and are there to *inspect* a hit you already found in
the table. A tight-looking cloud with a `p_value` near 1 is still noise.

### 5. View a plot

List the plots and open the protein you want to look at:

```bash
ls ~/results/experiment/plots/
```

On macOS or Linux, open the file directly — `open FILE.html`, `xdg-open
FILE.html`, or just double-click it. On WSL (a Windows machine running Linux),
hand it to your Windows browser instead:

```bash
explorer.exe "$(wslpath -w ~/results/experiment/plots/PROTEIN_ACCESSION_pocket.html)"
```

`wslpath -w` rewrites the Linux path into the form Windows needs; the file then
opens in your default browser. What each colour means is in [Output](#output)
below.

Downloaded structures cache under the output directory by default, so a repeat
run on the same export starts immediately. To share one cache across several
experiments, point them all at the same `--cache ~/some/dir`.

## Describing your scaffold

The scaffold's reach is derived from its chemistry rather than guessed, so this
is the one thing you have to supply beyond the mzIdentML file. There are two
routes. Run `covalent-pocket-mapper` with no arguments and it walks you through
either one.

### Route 1: the per-bond form (most faithful)

Look up three numbers for **each backbone bond in one repeat unit**:

| what | where it comes from |
|---|---|
| **bond length** | a standard chemistry table: C–C 1.53 Å, C–O 1.43 Å, C–N 1.47 Å, amide C–N 1.33 Å, C=C 1.34 Å, aromatic 1.39 Å |
| **bond angle** | the same table: a single-bonded carbon or oxygen sits near 109–112°, an sp² centre near 120° |
| **twist** | the bond's chemistry: an ordinary single bond is **RIS**, a double bond, amide, or aromatic bond is **rigid** |

For each RIS bond you also need its **rotamer energies** — how much a kink costs
versus staying straight, in kcal/mol. This is the only number not on a geometry
table; it comes from a polymer reference such as Flory's *Statistical Mechanics
of Chain Molecules* or the Polymer Handbook. A typical alkane gauche penalty is
0.5–0.9 kcal/mol. PEG is the interesting exception: its O–C–C–O gauche effect
runs the other way and *favours* the kink.

Finally, say how many repeat units the scaffold has — 8 for PEG8.

That is the whole definition. No structure file, no drawing, no simulation
input deck. Saved as JSON it looks like this (`--polymer peg8.json`):

```json
{
  "name": "PEG",
  "temperature_K": 298,
  "n_units": 8,
  "bonds": [
    {"length_A": 1.53, "angle_deg": 109.5,
     "torsion": {"type": "ris", "states_deg": [180, 60, -60], "energies_kcal": [0.0, 0.5, 0.5]}},
    {"length_A": 1.43, "angle_deg": 111.5,
     "torsion": {"type": "ris", "states_deg": [180, 60, -60], "energies_kcal": [0.3, 0.0, 0.0]}},
    {"length_A": 1.43, "angle_deg": 111.5,
     "torsion": {"type": "ris", "states_deg": [180, 60, -60], "energies_kcal": [0.0, 0.5, 0.5]}}
  ]
}
```

A `rigid` bond takes `{"type": "rigid", "fixed_deg": 180}` and needs no
energies; `{"type": "free"}` is an unhindered bond with no barrier at all.

### Route 2: the persistence-length shortcut

If you cannot find rotamer energies for your polymer, skip the per-bond detail
entirely. Look up one number — the **persistence length**, tabulated for
essentially everything common — and give the scaffold's total contour length:

```bash
covalent-pocket-mapper --persistence-length 3.8 --scaffold-size 28.9 ...
```

You lose a little accuracy at the very-short-chain extreme, where the smooth
worm-like-chain approximation is least apt, but it costs one lookup instead of
a form.

### What this replaces

Earlier versions took a scaffold *size* and a hand-tuned *rigidity* slider. Both
are gone. Size is now the contour length the bonds imply, and rigidity is no
longer a number you pick — it emerges from the rotamer energies, which is where
it actually comes from. Worth knowing what the physics says: **PEG8 has a
28.9 Å contour length but a most-likely end-to-end distance of only 10.6 Å.**
Using the stretched-out length as a reach radius, as the old model invited, is
off by nearly a factor of three.

## How well is the pocket actually located?

Multiplying the reach clouds gives a posterior over locations, not a single
point, and the tool reports the spread honestly as `credible_region_A` — the
width of the smallest set of surface points holding half the posterior mass.

This is where the polymer model earns its keep. A hard reach radius is flat
inside its own shell: if every capture site already sits within reach, every
candidate centre scores identically and the pocket is genuinely unlocatable. On
a real run the old model left RPS14 indeterminate across **56 Å**. A P(r)
likelihood is never flat — it peaks and falls away — so the same data now
resolve to a credible region, and proteins with more sites tighten sharply:

| protein | sites | credible region |
|---|---|---|
| ENO1 | 5 | 3.7 Å |
| ANXA2 | 7 | 5.0 Å |
| ACTB | 5 | 8.1 Å |
| RPS14 | 4 | 19.2 Å |

Four sites on a small protein still leave real ambiguity, and the column says
so rather than printing a coordinate that implies otherwise.

**One trap worth knowing.** `site_posterior` is the posterior mass sitting on
the single best point, and it is conditional on a pocket existing at all. A
protein where only one surface point is even reachable gets a posterior near
1.0 — that means the feasible set is tiny, not that the evidence is strong.
`p_value` is what measures evidence. Read the two together.


## Running without a control

A negative control is optional. The pipeline changes shape depending on whether
you have one:

**With a control.** A peptide counts as signal only if it appears in the
treatment runs and never in the control. Control peptides are kept and scored
separately: they are measured against the same fitted pocket, so you can see
whether non-specific capture lands in the same place. This is reported as
`control_fraction_reachable` and drawn in red in the plots. It matters more here
than it used to: because a single unreachable site zeroes out a protein, the
control is also what keeps a spurious site from vetoing a real binder.

**Without a control.** Every confidently detected peptide counts as signal, and
the geometry alone decides. This still works — the significance test is internal
to each protein, asking whether *these* sites are explained by a single pocket
better than a random scatter of the same size would be.

It is noticeably harsher, though, and it is worth knowing why. Without a control
a protein carries every peptide it ever yielded rather than only the specific
ones: on a real run ENO1 went from 5 sites to 20. Since one unreachable site
vetoes the protein, more sites means more chances to be vetoed, and the share of
proteins returning a pocket fell from roughly half to about a fifth.

`--min-psms` is the lever. Requiring more spectra per peptide drops the
one-off identifications that are most likely spurious, and the hit rate climbs
back: on the same data, 6 of 28 proteins at `--min-psms 1` versus 7 of 17 at
`--min-psms 3`. If you have a control, use it — it removes bad sites on evidence
rather than on abundance.

Either way, `z_score` and `p_value` are the primary ranking. A control adds a
second, independent line of evidence rather than being load-bearing.

## Output

Written to your chosen output directory:

| file | contents |
|---|---|
| `binding_site_scores.csv` | one row per scored protein, ranked most convincing first |
| `peptide_sites.csv` | every peptide used, with its position and counts |
| `skipped_proteins.csv` | proteins that could not be scored, and why |
| `plots/*.html` | interactive 3D views, one per top-ranked protein |
| `structures/` | cached PDB downloads, reused on later runs |

Columns worth knowing in `binding_site_scores.csv`:

- **`binding_site_x/y/z`** — the most likely pocket, in the coordinate frame of
  that protein's PDB file, so you can open the structure and look at it.
- **`p_value`** — share of random site scatters that explain themselves at least
  as well as the real ones do. **This is the evidence, and the main ranking.**
- **`credible_region_A`** — width of the smallest region holding half the
  posterior. Small means the pocket is pinned down; large means the coordinate
  above is the best point in a broad region.
- **`site_probability`** — the reach density at the pocket, as a per-site
  geometric mean (1/Å³), so proteins with different numbers of sites stay
  comparable.
- **`site_posterior`** — posterior mass on the single best point. Conditional on
  a pocket existing at all: see the trap noted above.
- **`log_likelihood`** — the summed log reach density; `-inf` means no surface
  point is reachable from every capture site.
- **`null_feasible_fraction`** — how often a random scatter of the same size
  admits any pocket at all. Low values mean the real pattern is doing something
  chance rarely manages, and it is often more telling than `z_score`.
- **`control_fraction_reachable`** — share of control sites also within reach of
  the fitted pocket, when you supply a control. Low is good: it means the pocket
  explains the real sites and not the background.
- **`mean_site_distance_A`** — average distance from the pocket to its capture
  sites, worth sanity-checking against your scaffold's most-likely span.

The plots are self-contained HTML — no server, no internet. Green markers are
capture sites, each tethered to the pocket by a dotted line standing in for the
scaffold; red markers are control sites; and the coloured cloud is the credible
region, shaded by likelihood, with the best point marked by a cross.

## Tuning

| flag | what it changes |
|---|---|
| `--min-psms` | spectra a peptide needs before it counts as signal |
| `--min-sites` | attachment points a protein needs to be scored at all (default 4) |
| `--grid-spacing` | resolution of the search over surface locations, in Å (default 1.5; finer is slower and rarely changes the answer) |
| `--plddt-floor` | lowest model confidence a residue may have (default 50; AlphaFold stores this in the B-factor column) |
| `--permutations` | random draws in the significance test (default 200; raise for finer p-values) |
| `--max-proteins` | stop after N proteins, for a quick trial run |
| `--no-download` | use only local structures |
| `--esmfold` | fold proteins AlphaFold lacks, one 400-residue window each |
| `--no-plots` | skip the 3D output |
| `--plot-top` | how many top-ranked proteins get a plot (default 5) |
| `--structures` | directory of local PDB files named by accession, e.g. `P08670.pdb` |
| `--cache` | where downloaded structures are kept; share one across runs to avoid re-downloading |
| `--seed` | fix the random seed so the permutation test is reproducible |
| `--no-control` | run without a negative control and don't prompt for one |

## How it works

1. **Parse** the mzIdentML export into peptide-to-protein matches with positions
   and sample of origin. Decoys and common contaminants are dropped.
2. **Label** each peptide as signal or control, per the section above.
3. **Fetch** a structure per protein — local directory, then cache, then
   AlphaFold, then optionally ESMFold. Every network step degrades to "no
   structure for this protein" rather than failing the run.
4. **Locate** each peptide's N-terminal residue in 3D. That is where the
   scaffold's reactive group attached, so it is the point the geometry is about.
   Residues below the confidence floor are skipped.
5. **Build** a reach cloud around each capture site from the scaffold's P(r),
   converted from a radial distribution to a density in space by dividing out
   the 4πr² shell factor. Without that division, a point right next to a capture
   site would look implausible purely because a thin shell has little area.
6. **Multiply** the clouds — summed in log space — over a grid of candidate
   locations in the solvent shell just outside the protein, since a covalent
   partner has to be on the surface. Interior cavities are excluded by a
   directional test: a buried point has backbone in every direction, so the mean
   direction to its neighbours cancels, while a surface point sees backbone only
   on one side. The threshold is deliberately permissive, because binding
   pockets are concave and a strict convexity test discards them.
7. **Test** against a null: draw the same number of sites from anywhere on the
   protein's surface and score them the same way, including the common and
   informative outcome that no single pocket can reach a random scatter at all.

   The null draws from surface residues rather than all residues, which matters
   more than it sounds. A scaffold can only reach the outside of a protein, and
   buried residues sit far closer together than surface ones — in a typical
   structure, roughly 13 Å from the centre versus 24 Å. A null that included
   them would be scoring against a tight interior cloud, making chance look
   better than it is and pushing every real score down. Burial is estimated by
   counting neighbouring alpha-carbons within 12 Å; the most crowded 40% are
   treated as interior and excluded.

## Limitations

- The N-terminus of a detected protein fragment is a proxy for the attachment site. It is
  the right proxy for scaffolds that react at a digestion-adjacent position, but
  it is a proxy. 
- A protein needs a structure and at least four attachment points inside the
  modelled region. Long disordered proteins often fail both.
- The null assumes attachment could have happened anywhere on the modelled
  surface. It does not model residue-level chemical preference, so a protein
  whose reactive residues are themselves clustered can score well without
  binding anything.
- Predicted structures are single conformations. A pocket that only exists in a
  bound conformation will not be found.
- With few capture sites the pocket stays broadly located even though the
  posterior is proper — four sites on a small protein can leave a credible
  region tens of angstroms wide. `credible_region_A` reports it.
- **One unreachable site vetoes the whole protein.** That is the intended
  behaviour — it is what makes the multiplication strict — but it means a single
  spurious capture site can zero out a real binder. A negative control is the
  defence, since it removes spurious sites before the geometry runs. On a real
  run about half the proteins came back with no reachable pocket, and the ones
  eliminated were largely the usual sticky non-specific binders (vimentin,
  filamin, talin, tropomyosin).
- All sites are weighted equally, which assumes one pocket per protein. A
  protein with two genuine binding sites will usually return no pocket at all
  rather than either one.
- The reach model treats the scaffold as a free chain. It does not know about
  the protein getting in the way, so a path that would have to pass through the
  protein is still counted as reachable.
