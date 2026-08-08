"""Read a Scaffold / mzIdentML export into a peptide-spectrum-match table.

The only required input to the whole pipeline is this file. Everything the
scoring needs -- peptide sequence, parent protein accession, the peptide's
start/stop position in that protein, and which sample it came from -- is read
from here.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd

# Scaffold writes the peptide's position into PeptideEvidence; a peptide shared
# between proteins appears once per parent, which is why one PSM can yield
# several rows.

_CONTAMINANT_PREFIXES = ("Cont_", "CONT_", "contam_")


def _open(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _namespace(root: ET.Element) -> dict[str, str]:
    match = re.match(r"\{(.+)\}", root.tag)
    return {"mz": match.group(1)} if match else {"mz": ""}


def clean_accession(accession: str | None) -> str | None:
    """Strip Scaffold's group-size suffix, e.g. 'P08195 (+1)' -> 'P08195'."""
    if not accession:
        return None
    return re.sub(r"\s*\(\+\d+\)", "", accession).strip()


def parse_mzid(path: str | Path, *, passing_only: bool = True) -> pd.DataFrame:
    """Parse an mzIdentML file into one row per peptide-to-protein match.

    Columns: peptide_sequence, accession, start, stop, charge, score,
    pass_threshold, is_decoy, sample.
    """
    path = Path(path)
    with _open(path) as handle:
        tree = ET.parse(handle)
    root = tree.getroot()
    ns = _namespace(root)

    peptide_seqs: dict[str, str] = {}
    for pep in root.findall(".//mz:Peptide", ns):
        seq_el = pep.find("mz:PeptideSequence", ns)
        if seq_el is not None and seq_el.text:
            peptide_seqs[pep.get("id")] = seq_el.text

    protein_accs: dict[str, str] = {}
    for dbseq in root.findall(".//mz:DBSequence", ns):
        protein_accs[dbseq.get("id")] = dbseq.get("accession") or dbseq.get("id")

    # SpectraData names the run each spectrum came from; the id is usually the
    # short run name, with the source file name as a last resort.
    sample_names: dict[str, str] = {}
    for spec in root.findall(".//mz:SpectraData", ns):
        sid = spec.get("id")
        location = spec.get("location")
        name = spec.get("name") or sid
        if not name and location:
            name = Path(location).stem
        sample_names[sid] = name or sid

    evidence: dict[str, dict] = {}
    for pev in root.findall(".//mz:PeptideEvidence", ns):
        db_ref = pev.get("dBSequence_ref")
        if db_ref is None:
            continue
        start, end = pev.get("start"), pev.get("end")
        evidence[pev.get("id")] = {
            "accession": protein_accs.get(db_ref, db_ref),
            "start": int(start) if start else None,
            "stop": int(end) if end else None,
            "is_decoy": pev.get("isDecoy", "false") == "true",
        }

    records = []
    for result in root.findall(".//mz:SpectrumIdentificationResult", ns):
        spectra_ref = result.get("spectraData_ref", "")
        sample = sample_names.get(spectra_ref, spectra_ref)
        for item in result.findall("mz:SpectrumIdentificationItem", ns):
            passes = item.get("passThreshold") == "true"
            if passing_only and not passes:
                continue
            score = None
            for cv in item.findall("mz:cvParam", ns):
                if "score" in (cv.get("name") or "").lower():
                    score = cv.get("value")
            sequence = peptide_seqs.get(item.get("peptide_ref"), "")
            for ref in item.findall("mz:PeptideEvidenceRef", ns):
                info = evidence.get(ref.get("peptideEvidence_ref"))
                if not info:
                    continue
                records.append(
                    {
                        "peptide_sequence": sequence,
                        "accession": info["accession"],
                        "start": info["start"],
                        "stop": info["stop"],
                        "charge": item.get("chargeState"),
                        "score": score,
                        "pass_threshold": passes,
                        "is_decoy": info["is_decoy"],
                        "sample": sample,
                    }
                )

    if not records:
        raise ValueError(
            f"No peptide-spectrum matches found in {path}. "
            "Is this an mzIdentML export?"
        )

    frame = pd.DataFrame(records)
    frame["accession"] = frame["accession"].map(clean_accession)
    return frame


def drop_contaminants_and_decoys(frame: pd.DataFrame) -> pd.DataFrame:
    """Remove decoy hits and the common-contaminant entries Scaffold prefixes."""
    keep = ~frame["is_decoy"]
    keep &= ~frame["accession"].fillna("").str.startswith(_CONTAMINANT_PREFIXES)
    return frame[keep].reset_index(drop=True)


def list_samples(frame: pd.DataFrame) -> list[str]:
    """Sample names present in the export, in a stable order."""
    return sorted(frame["sample"].dropna().unique().tolist())


def parse_proteins(path: str | Path) -> pd.DataFrame:
    """Per-protein metadata from the same file: description, gene, length.

    Scaffold carries the FASTA header through as a 'protein description' term,
    so gene symbols come free with the export and need no lookup.
    """
    path = Path(path)
    with _open(path) as handle:
        tree = ET.parse(handle)
    root = tree.getroot()
    ns = _namespace(root)

    records = []
    for dbseq in root.findall(".//mz:DBSequence", ns):
        description = None
        for cv in dbseq.findall("mz:cvParam", ns):
            if cv.get("name") == "protein description":
                description = cv.get("value")
        length = dbseq.get("length")
        gene_match = re.search(r"\bGN=(\S+)", description or "")
        records.append(
            {
                "accession": clean_accession(dbseq.get("accession") or dbseq.get("id")),
                "gene": gene_match.group(1) if gene_match else None,
                "description": description,
                "sequence_length": int(length) if length else None,
            }
        )

    frame = pd.DataFrame(records).drop_duplicates(subset="accession")
    return frame.reset_index(drop=True)
