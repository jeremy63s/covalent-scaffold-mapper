"""Get per-residue 3D coordinates for each protein.

Structures are looked for in this order: a local directory the user points at,
then the download cache, then the AlphaFold database. Very long proteins are
stored in AlphaFold as overlapping fragments, which are stitched into one
coordinate table keeping the highest-confidence copy of each residue. ESMFold is
available as an opt-in last resort for proteins AlphaFold does not cover.

Every network step degrades to "no structure for this protein" rather than
raising, so one unreachable service cannot end the run.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import pandas as pd

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"
ESMFOLD_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"

# ESMFold refuses sequences longer than this, so long proteins are folded as a
# single window chosen to cover as many detected peptides as possible.
ESMFOLD_MAX_RESIDUES = 400

# AlphaFold writes its per-residue confidence (pLDDT, 0-100) into the B-factor
# column. Residues below this are too disordered to trust a coordinate for.
DEFAULT_PLDDT_FLOOR = 50.0

COORD_COLUMNS = ["residue_num", "residue_name", "x", "y", "z", "plddt"]

# AlphaFold names its models AF-<accession>-F<n> for the canonical sequence and
# AF-<accession>-<isoform>-F<n> for alternative isoforms. Only fragments of one
# isoform share a residue numbering, so mixing them would build a chimera that
# does not exist. The canonical sequence is what mzIdentML positions refer to.
_CANONICAL_ENTRY = r"AF-{acc}-F\d+"
_ANY_AF_ENTRY = r"AF-.+?-(?:\d+-)?F\d+"


def extract_ca_coords(pdb_path: str | Path) -> pd.DataFrame:
    """Read alpha-carbon coordinates and confidence from a PDB file."""
    from Bio import PDB  # imported lazily so the module loads without biopython

    parser = PDB.PDBParser(QUIET=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = parser.get_structure("protein", str(pdb_path))

    residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue) and "CA" in residue:
                    ca = residue["CA"]
                    x, y, z = (float(v) for v in ca.get_coord())
                    residues.append(
                        {
                            "residue_num": residue.get_id()[1],
                            "residue_name": residue.get_resname(),
                            "x": x,
                            "y": y,
                            "z": z,
                            "plddt": float(ca.get_bfactor()),
                        }
                    )
        break  # a predicted structure has one model; take the first

    return pd.DataFrame(residues, columns=COORD_COLUMNS)


def _merge_fragments(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Combine fragment coordinate tables, keeping the best copy per residue."""
    merged = pd.concat(frames, ignore_index=True)
    merged = merged.sort_values(["residue_num", "plddt"], ascending=[True, False])
    merged = merged.drop_duplicates(subset="residue_num", keep="first")
    return merged.sort_values("residue_num").reset_index(drop=True)


def fetch_uniprot_sequence(accession: str, *, timeout: int = 15) -> str | None:
    """Return the protein's amino-acid sequence, or None if it cannot be read."""
    import requests

    try:
        response = requests.get(UNIPROT_FASTA.format(acc=accession), timeout=timeout)
        if not response.ok:
            return None
        lines = response.text.strip().split("\n")
        return "".join(lines[1:]) or None
    except requests.RequestException:
        return None


class StructureProvider:
    """Resolves an accession to a coordinate table, with caching."""

    def __init__(
        self,
        cache_dir: Path,
        *,
        local_dir: Path | None = None,
        allow_download: bool = True,
        allow_esmfold: bool = False,
        timeout: int = 30,
    ):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.local_dir = Path(local_dir) if local_dir else None
        self.allow_download = allow_download
        self.allow_esmfold = allow_esmfold
        self.timeout = timeout
        self.notes: dict[str, str] = {}

    # -- lookup steps ----------------------------------------------------

    @staticmethod
    def _keep(path: Path, accession: str) -> bool:
        """Drop alternative-isoform models; keep canonical and user-named files."""
        stem = path.stem
        if re.fullmatch(_ANY_AF_ENTRY, stem):
            return bool(re.fullmatch(_CANONICAL_ENTRY.format(acc=re.escape(accession)), stem))
        return True

    def _find(self, directory: Path, accession: str) -> list[Path]:
        """Structure files for this accession, whatever naming scheme they use."""
        matches = sorted(directory.glob(f"*{accession}*.pdb"))
        matches += sorted(directory.glob(f"*{accession}*.cif"))
        return [p for p in matches if p.is_file() and self._keep(p, accession)]

    def _local_files(self, accession: str) -> list[Path]:
        """Structure files sitting in the directory the user pointed at."""
        if not self.local_dir:
            return []
        return self._find(self.local_dir, accession)

    def _cached_files(self, accession: str) -> list[Path]:
        return self._find(self.cache_dir, accession)

    def _download_alphafold(self, accession: str) -> list[Path]:
        """Download every AlphaFold fragment for this accession into the cache."""
        import requests

        try:
            response = requests.get(
                ALPHAFOLD_API.format(acc=accession), timeout=self.timeout
            )
            if not response.ok:
                return []
            entries = response.json()
        except (requests.RequestException, ValueError):
            return []

        saved = []
        canonical = re.compile(_CANONICAL_ENTRY.format(acc=re.escape(accession)))
        for entry in entries:
            pdb_url = entry.get("pdbUrl")
            entry_id = entry.get("entryId") or accession
            if not pdb_url:
                continue
            # The API returns every isoform; only the canonical one is numbered
            # the way the mzIdentML peptide positions are.
            if not canonical.fullmatch(entry_id):
                continue
            out_path = self.cache_dir / f"{entry_id}.pdb"
            if not out_path.exists():
                try:
                    pdb_response = requests.get(pdb_url, timeout=self.timeout)
                    if not pdb_response.ok:
                        continue
                    out_path.write_text(pdb_response.text)
                except requests.RequestException:
                    continue
            saved.append(out_path)
        return saved

    def _fold_with_esmfold(
        self, accession: str, focus_positions: list[int] | None
    ) -> pd.DataFrame | None:
        """Fold a single window of the sequence, centred on the detected sites."""
        import requests

        sequence = fetch_uniprot_sequence(accession, timeout=self.timeout)
        if not sequence:
            return None

        offset = 0
        if len(sequence) > ESMFOLD_MAX_RESIDUES:
            offset = _densest_window_start(
                len(sequence), focus_positions or [], ESMFOLD_MAX_RESIDUES
            )
            sequence = sequence[offset : offset + ESMFOLD_MAX_RESIDUES]

        try:
            response = requests.post(ESMFOLD_API, data=sequence, timeout=300)
        except requests.RequestException as exc:
            self.notes[accession] = f"ESMFold unreachable ({exc.__class__.__name__})"
            return None
        if not response.ok:
            self.notes[accession] = f"ESMFold declined the sequence ({response.status_code})"
            return None

        out_path = self.cache_dir / f"{accession}_ESM.pdb"
        out_path.write_text(response.text)
        coords = extract_ca_coords(out_path)
        if coords.empty:
            return None
        # ESMFold numbers the window from 1; shift back onto the full protein.
        coords["residue_num"] = coords["residue_num"] + offset
        self.notes[accession] = (
            f"ESMFold model of residues {offset + 1}-{offset + len(sequence)}"
        )
        return coords

    # -- public API ------------------------------------------------------

    def get(
        self, accession: str, focus_positions: list[int] | None = None
    ) -> pd.DataFrame | None:
        """Coordinates for one protein, or None if no structure could be had."""
        files = self._local_files(accession) or self._cached_files(accession)
        source = "local file" if files else None

        if not files and self.allow_download:
            files = self._download_alphafold(accession)
            source = "AlphaFold" if files else None

        frames = []
        for path in files:
            try:
                frame = extract_ca_coords(path)
            except Exception as exc:  # a malformed file should not stop the run
                self.notes[accession] = f"Could not parse {path.name}: {exc}"
                continue
            if not frame.empty:
                frames.append(frame)

        if frames:
            coords = _merge_fragments(frames)
            if source:
                note = f"{source} ({len(frames)} model{'s' if len(frames) > 1 else ''})"
                self.notes.setdefault(accession, note)
            return coords

        if self.allow_esmfold:
            return self._fold_with_esmfold(accession, focus_positions)

        self.notes.setdefault(accession, "no structure available")
        return None


def _densest_window_start(
    sequence_length: int, positions: list[int], window: int
) -> int:
    """Zero-based start of the window covering the most detected positions."""
    if sequence_length <= window or not positions:
        return 0
    best_start, best_count = 0, -1
    for start in range(0, sequence_length - window + 1, 10):
        end = start + window
        count = sum(1 for p in positions if start < p <= end)
        if count > best_count:
            best_start, best_count = start, count
    return best_start


def coverage_fraction(coords: pd.DataFrame, sequence_length: int | None) -> float | None:
    """How much of the protein the structure actually spans."""
    if not sequence_length or coords is None or coords.empty:
        return None
    return len(coords) / sequence_length
