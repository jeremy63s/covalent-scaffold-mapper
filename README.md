# scaffold-binding

Finds candidate covalent binding pockets from scaffold pulldown mass spectrometry.

When a chemical scaffold with a reactive group is pulled down with its bound
proteins and the result is digested and sequenced, the peptides that come back
tell you *where* on each protein the scaffold attached. If a protein is a real
binder, those attachment points should not be scattered at random over its
surface — they should sit within one scaffold's reach of a single pocket. This
tool tests that, protein by protein, against a null of randomly placed
attachment points on the same structure.

## What you need

one mzIdentML export from your search — `.mzid` or `.mzid.gz`.
In Scaffold this is `Export > mzIdentML`. That single file carries everything
the analysis uses: peptide sequences, the protein each maps to, the peptide's
start and stop position within that protein, and which run each spectrum came
from. No spreadsheet exports, peptide CSVs, or FASTA files are needed.

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
git clone https://github.com/jeremy63s/covalent-scaffold-mapper
cd ~/covalent-scaffold-mapper

python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

pip install .
```

That installs the dependencies and puts a `scaffold-binding` command on your
path, which works from any directory. Check it came through:

```bash
scaffold-binding --help
```

If you plan to edit the code, use `pip install -e .` instead so your changes
take effect without reinstalling.

If you would rather not install anything, `pip install -r requirements.txt` and
then run `python -m scaffold_binding` from inside this directory. Every command
below works the same way, just with that longer form.

## Run

Start here. With no arguments it asks for everything it needs, one question at a
time, and explains the two scaffold parameters as it goes:

```bash
scaffold-binding
```

Once you know what you want, give it everything up front:

```bash
scaffold-binding \
  --mzid ~/data/experiment.mzid.gz \
  --treatment run_01 run_02 \
  --control run_03 \
  --output ~/results/experiment \
  --scaffold-size 25 \
  --rigidity 0.7
```

Your sample names come from the export itself, so check them before choosing:

```bash
scaffold-binding --mzid ~/data/experiment.mzid.gz --list-samples
```

### A first run worth doing

Structures are downloaded on demand, so a full run on a large export takes a
while the first time. Try a handful of proteins first to confirm everything
works and to see whether your scaffold settings behave sensibly:

```bash
scaffold-binding --mzid ~/data/experiment.mzid.gz \
  --treatment run_01 --no-control \
  --output ~/results/trial \
  --scaffold-size 25 --rigidity 0.7 \
  --max-proteins 5 --permutations 50
```

Downloaded structures are cached under the output directory, so point later runs
at the same `--cache` and they start immediately.

Results, the structure cache, and
plots all go under the output directory you select.

## The two scaffold parameters

These are the knobs that encode your chemistry. Both are yours to set; neither
is fitted from the data.

### `--scaffold-size` (angstroms)

How far the scaffold reaches from its anchor point. This becomes the radius of
the shell fitted around each candidate pocket. Because it is fixed by the
molecule rather than fitted, only the pocket *centre* is optimised — which is
what makes the resulting number interpretable rather than a best-fit artifact.

### `--rigidity` (0 to 1)

How much the scaffold can fold back on itself, which decides how far *inside*
the shell an attachment point may sit and still be consistent with the pocket.

A flexible scaffold with a PEG8 linker can bend, so a site anywhere closer than
its full length is perfectly plausible and its rigidity would be 0. A rigid scaffold holds its length, so a
site well inside the shell means the pocket is probably somewhere else.

With `R` = scaffold size and `d` = distance from a site to the pocket centre:

| rigidity | unpenalised range | meaning |
|---|---|---|
| `0.0` | `0 ≤ d ≤ R` | fully flexible; nothing inside the shell is penalised |
| `0.5` | `0.5R ≤ d ≤ R` | moderate |
| `0.7` | `0.7R ≤ d ≤ R` | sites closer in than `0.7R` are penalised by how far past that inner ring they sit |
| `1.0` | `d = R` | fully rigid; anything off the shell surface is penalised |

Worked example at `--scaffold-size 20 --rigidity 0.7`: the inner ring is at
14 Å. A site 16 Å from the pocket centre is free. One at 10 Å is penalised by
4 Å (how far inside 14 Å it fell). One at 5 Å is penalised by 9 Å.

**Sites outside the shell are penalised by their distance beyond it, at every
rigidity.** No amount of flexibility lets a scaffold reach further than its own
length, so that side of the penalty never changes.

## Running without a control

A negative control is optional. The pipeline changes shape depending on whether
you have one:

**With a control.** A peptide counts as signal only if it appears in the
treatment runs and never in the control. Control peptides are kept and scored
separately: they are measured against the same fitted pocket, so you can see
whether non-specific capture lands in the same place. This is reported as
`signal_vs_control` and drawn in red in the plots.

**Without a control.** Every confidently detected peptide counts as signal, and
the geometry alone decides. This still works because the significance test is
internal to each protein: it asks whether *these* attachment points cluster
better than the same number of points scattered anywhere on the same protein's
surface. Use `--min-psms` to require more spectra per peptide, since it becomes
your main abundance filter.

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

- **`fraction_reachable`** — share of this protein's attachment points the
  scaffold could reach from the fitted pocket.
- **`mean_penalty`** — average angstroms by which points miss the reachable
  zone. Zero is a perfect fit.
- **`z_score`** — how much better the real points fit than random points on the
  same protein. Higher is better; this is the main ranking.
- **`p_value`** — share of random draws that fit at least as well.
- **`signal_vs_control`** — present only when you supply a control.
- **`center_x/y/z`** — coordinates of the inferred pocket, in the frame of that
  protein's PDB file, so you can open the structure and look at it.

The plots are self-contained HTML — no server, no internet. Green points are in
reach of the fitted pocket, amber are not, red are control sites. The blue shell
is the scaffold's reach and the purple inner surface (drawn when rigidity is
between 0 and 1) is the inner limit.

## Tuning

| flag | what it changes |
|---|---|
| `--min-psms` | spectra a peptide needs before it counts as signal |
| `--min-sites` | attachment points a protein needs to be scored at all (default 4 — three points define a sphere, so fewer cannot be tested) |
| `--plddt-floor` | lowest model confidence a residue may have (default 50; AlphaFold stores this in the B-factor column) |
| `--permutations` | random draws in the significance test (default 200; raise for finer p-values) |
| `--max-proteins` | stop after N proteins, for a quick trial run |
| `--no-download` | use only local structures |
| `--esmfold` | fold proteins AlphaFold lacks, one 400-residue window each |
| `--no-plots` | skip the 3D output |

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
5. **Fit** a pocket centre minimising the mean penalty over those points, with
   the radius held at your scaffold size. Several starting positions are tried
   so the optimiser does not settle into a local minimum.
6. **Test** against a null: draw the same number of points from residues
   anywhere on that protein's surface, refit a pocket from scratch with the same
   objective, and repeat. `z_score` compares the real fit to that distribution.

   The null draws from surface residues rather than all residues, which matters
   more than it sounds. A scaffold can only reach the outside of a protein, and
   buried residues sit far closer together than surface ones — in a typical
   structure, roughly 13 Å from the centre versus 24 Å. A null that included
   them would be fitting spheres to a tight interior cloud, making chance look
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
