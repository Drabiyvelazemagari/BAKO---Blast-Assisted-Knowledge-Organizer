#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import csv
import json
import hashlib
import shutil
import subprocess
import time
import random
import urllib.parse
import urllib.request
import urllib.error
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from PyQt6 import QtWidgets, QtCore

# -----------------------------
# App metadata
# -----------------------------
__app_name__ = "BAKO"
__version__ = "2.0.0"
__build__ = "pyqt6-crossplatform"
__repo__ = ""  # optional: add repo URL




# -----------------------------
# Utilities
# -----------------------------

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

def safe_write_text(path: Path, s: str) -> None:
    """Write text atomically (best-effort).

    Writes to a temporary file in the same directory and then replaces the target.
    This prevents partial/corrupt files if the process is interrupted.
    """
    ensure_dir(path.parent)
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    with tmp.open("w", encoding="utf-8", newline="\n") as f:
        f.write(s)
    os.replace(tmp, path)

def which_or_raise(bin_name: str) -> str:
    """
    Return full path to a required external executable on PATH, or raise with OS-specific guidance.
    On Windows, we also try '<name>.exe' for convenience.
    """
    p = shutil.which(bin_name)

    if not p and sys.platform.startswith("win") and not bin_name.lower().endswith(".exe"):
        p = shutil.which(bin_name + ".exe")

    if not p:
        if sys.platform.startswith("win"):
            raise RuntimeError(
                f"'{bin_name}' not found on PATH.\n\n"
                "Windows install options:\n"
                "  • Download the binary (e.g., prodigal.exe) and add its folder to PATH\n"
                "    Example: put it in C:\\tools\\prodigal\\ and add C:\\tools\\prodigal\\ to PATH\n"
                "  • Or install via Conda: conda install -c bioconda prodigal\n\n"
                "After installing, open a NEW terminal and confirm with:\n"
                f"  {bin_name} -v\n"
            )
        raise RuntimeError(
            f"'{bin_name}' not found on PATH.\n"
            f"macOS: brew install {bin_name}\n"
            "Linux: install via your package manager\n"
        )
    return p

def run_cmd(cmd: List[str]) -> str:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out = p.stdout or ""
    if p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {' '.join(cmd)}\n\n{out}")
    return out

def default_cache_root() -> Path:
    if os.name == "nt":
        appdata = os.environ.get("APPDATA")
        if appdata:
            return Path(appdata) / "BAKO" / "cds_cache"
    return Path.home() / ".bako" / "cds_cache"

def default_meta_cache_path() -> Path:
    """Stable location for NCBI metadata cache (JSON)."""
    if os.name == "nt":
        appdata = os.environ.get("APPDATA")
        if appdata:
            return Path(appdata) / "BAKO" / "ncbi_meta_cache.json"
    return Path.home() / ".bako" / "ncbi_meta_cache.json"


# -----------------------------
# NCBI metadata enrichment (replicon + taxonomy)
# -----------------------------

NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

NCBI_FAIL_COOLDOWN_S = 6 * 60 * 60  # retry failed NCBI lookups after 6 hours

def _http_get_json(url: str, timeout_s: int = 30) -> dict:
    """HTTP GET JSON with retries, backoff, and timeouts.

    Designed for NCBI E-utilities (transient network errors and rate limits are common).
    Raises on non-retriable errors.
    """
    # Conservative retry policy: transient network failures, 429, and 5xx get retried.
    max_attempts = 4
    backoff_base = 0.8
    backoff_max = 12.0

    last_err: Optional[Exception] = None
    for attempt in range(1, max_attempts + 1):
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "User-Agent": f"{__app_name__}/{__version__} ({__build__}; metadata enrichment)",
                    "Accept": "application/json",
                },
            )
            with urllib.request.urlopen(req, timeout=timeout_s) as resp:
                data = resp.read().decode("utf-8", errors="replace")
            return json.loads(data)
        except urllib.error.HTTPError as e:
            # Retriable status codes
            status = getattr(e, "code", None)
            retriable = status in (429,) or (isinstance(status, int) and 500 <= status <= 599)
            last_err = e
            if not retriable or attempt == max_attempts:
                raise
        except Exception as e:
            last_err = e
            if attempt == max_attempts:
                raise

        # Backoff with jitter (polite + reduces thundering herd)
        sleep_s = min(backoff_max, backoff_base * (2 ** (attempt - 1)))
        sleep_s *= (1.0 + random.random() * 0.25)
        time.sleep(sleep_s)

    # Should never get here
    if last_err:
        raise last_err
    raise RuntimeError("HTTP request failed unexpectedly")

def _infer_replicon_from_text(txt: str) -> Tuple[str, str]:
    t = (txt or "").lower()
    if "plasmid" in t:
        return "plasmid", "high"
    if "chromosome" in t or "complete genome" in t:
        return "chromosome", "high"
    if "bacteriophage" in t or "phage" in t or "virus" in t or "prophage" in t:
        return "phage", "high"
    if "contig" in t or "scaffold" in t:
        return "contig", "medium"
    return "unknown", "low"

def _split_genus_species(organism_name: str) -> Tuple[Optional[str], Optional[str]]:
    if not organism_name:
        return None, None
    parts = organism_name.split()
    if len(parts) == 0:
        return None, None
    genus = parts[0]
    if len(parts) < 2:
        return genus, None
    sp = parts[1]
    if sp.lower() in {"sp.", "sp"}:
        return genus, None
    return genus, sp

def load_meta_cache(path: Path) -> Dict[str, dict]:
    try:
        if path.exists():
            return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        pass
    return {}

def save_meta_cache(path: Path, cache: Dict[str, dict]) -> None:
    ensure_dir(path.parent)
    safe_write_text(path, json.dumps(cache, ensure_ascii=False, indent=2))

def fetch_ncbi_meta_for_accession(
    accession: str,
    email: str = "",
    api_key: str = "",
    polite_delay_s: float = 0.34,
) -> Optional[dict]:
    """Fetch minimal metadata for a nuccore accession via NCBI E-utilities.

    Returns dict with keys:
      organism_name, genus, species, taxid, replicon_type, replicon_confidence, replicon_source

    Notes:
    - This is best-effort; failures return None.
    - Replicon type is inferred from title/defline keywords when available.
    """
    acc = (accession or "").strip()
    if not acc:
        return None

    # ESearch: accession -> UID
    term = f"{acc}[Accession]"
    params = {"db": "nuccore", "term": term, "retmode": "json"}
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    url = f"{NCBI_EUTILS_BASE}/esearch.fcgi?{urllib.parse.urlencode(params)}"

    try:
        js = _http_get_json(url)
    except Exception:
        return None

    try:
        idlist = js.get("esearchresult", {}).get("idlist", [])
        uid = idlist[0] if idlist else None
    except Exception:
        uid = None

    time.sleep(polite_delay_s)

    if not uid:
        return None

    # ESummary: UID -> title/taxname/taxid (shape can vary)
    params2 = {"db": "nuccore", "id": uid, "retmode": "json"}
    if email:
        params2["email"] = email
    if api_key:
        params2["api_key"] = api_key
    url2 = f"{NCBI_EUTILS_BASE}/esummary.fcgi?{urllib.parse.urlencode(params2)}"

    try:
        js2 = _http_get_json(url2)
    except Exception:
        return None

    time.sleep(polite_delay_s)

    doc = {}
    try:
        doc = js2.get("result", {}).get(str(uid), {})
    except Exception:
        doc = {}

    title = doc.get("title", "") or ""
    organism = doc.get("organism", "") or doc.get("taxname", "") or ""
    taxid = doc.get("taxid", None)

    replicon_type, conf = _infer_replicon_from_text(title)
    genus, species = _split_genus_species(organism)

    return {
        "organism_name": organism or None,
        "genus": genus,
        "species": species,
        "taxid": taxid,
        "replicon_type": replicon_type,
        "replicon_confidence": conf,
        "replicon_source": "ncbi",
        "ncbi_title": title or None,
        "ncbi_uid": str(uid),
    }

class MetaResolver:
    def __init__(self, cache_path: Path, log, email: str = "", api_key: str = "", enabled: bool = False):
        self.cache_path = cache_path
        self.log = log
        self.email = (email or "").strip()
        self.api_key = (api_key or "").strip()
        self.enabled = bool(enabled)
        self.cache: Dict[str, dict] = load_meta_cache(cache_path)
        self.hits = 0
        self.misses = 0
        self.failures = 0

    def resolve(self, accession: str) -> dict:
        # Always return a dict with expected keys (even if disabled)
        empty = {
            "replicon_type": "unknown",
            "replicon_confidence": "low",
            "replicon_source": "",
            "organism_name": None,
            "genus": None,
            "species": None,
            "taxid": None,
        }
        acc = (accession or "").strip()
        if not acc:
            return empty

        # Cache hit (including temporary failures)
        if acc in self.cache:
            m = self.cache[acc] or {}
            # If we previously failed to fetch from NCBI, avoid permanent poison:
            # retry after a cooldown.
            if m.get("_status") == "failed":
                ts = int(m.get("_ts") or 0)
                if int(time.time()) - ts < NCBI_FAIL_COOLDOWN_S:
                    self.hits += 1
                    return empty
                # Cooldown elapsed: treat as a miss and retry
            else:
                self.hits += 1
                out = empty.copy()
                out.update({k: m.get(k, out.get(k)) for k in out.keys()})
                return out

        self.misses += 1
        if not self.enabled:
            return empty

        meta = fetch_ncbi_meta_for_accession(acc, email=self.email, api_key=self.api_key)
        if not meta:
            # Cache a *temporary* failure marker (do not poison the cache forever).
            self.failures += 1
            fail_rec = empty.copy()
            fail_rec.update({"_status": "failed", "_ts": int(time.time())})
            self.cache[acc] = fail_rec
            save_meta_cache(self.cache_path, self.cache)
            return empty

        # Normalize to our output keys
        out = empty.copy()
        out.update({
            "replicon_type": meta.get("replicon_type", "unknown"),
            "replicon_confidence": meta.get("replicon_confidence", "low"),
            "replicon_source": meta.get("replicon_source", "ncbi"),
            "organism_name": meta.get("organism_name", None),
            "genus": meta.get("genus", None),
            "species": meta.get("species", None),
            "taxid": meta.get("taxid", None),
        })
        # Store extra debug fields too, but they won't be written to TSV unless you add columns
        self.cache[acc] = {**out, "ncbi_uid": meta.get("ncbi_uid"), "ncbi_title": meta.get("ncbi_title")}
        save_meta_cache(self.cache_path, self.cache)
        return out

    def summarize(self) -> str:
        if not self.enabled:
            return "[META] NCBI enrichment disabled."
        return f"[META] Cache hits={self.hits}  misses={self.misses}  failures={self.failures}  cache_file={self.cache_path}"


# -----------------------------
# Offline metadata from FASTA headers (fast, less reliable)
# -----------------------------

def _extract_bracket_field(header: str, key: str) -> Optional[str]:
    """Extract values like [key=value] from an NCBI-style FASTA header."""
    if not header:
        return None
    m = re.search(r"\[" + re.escape(key) + r"=([^\]]+)\]", header)
    if not m:
        return None
    val = m.group(1).strip()
    return val or None

def offline_meta_from_header(header: str) -> dict:
    """Best-effort offline metadata: infer replicon + taxonomy from FASTA header text.

    Returns the same keys as MetaResolver.resolve().
    Confidence is lower than NCBI enrichment and depends on the header content.
    """
    empty = {
        "replicon_type": "unknown",
        "replicon_confidence": "low",
        "replicon_source": "",
        "organism_name": None,
        "genus": None,
        "species": None,
        "taxid": None,
    }

    h = (header or "").strip()
    if not h:
        return empty

    # Replicon inference from the full header text
    rep, conf = _infer_replicon_from_text(h)

    # Try NCBI-style bracket fields first
    organism = _extract_bracket_field(h, "organism") or _extract_bracket_field(h, "taxname")
    taxid_raw = _extract_bracket_field(h, "taxid") or _extract_bracket_field(h, "taxon")

    taxid: Optional[int] = None
    if taxid_raw:
        try:
            taxid = int(re.sub(r"\D+", "", taxid_raw))
        except Exception:
            taxid = None

    genus, species = _split_genus_species(organism or "")

    out = empty.copy()
    out.update({
        "replicon_type": rep,
        "replicon_confidence": conf if conf in {"high", "medium", "low"} else "low",
        "replicon_source": "header",
        "organism_name": organism or None,
        "genus": genus,
        "species": species,
        "taxid": taxid,
    })
    return out




def is_accession_like(s: str) -> bool:
    """Best-effort check for NCBI-like accessions (to avoid pointless network calls)."""
    if not s:
        return False
    s = s.strip()
    low = s.lower()
    # Common local assembly contig patterns
    if low.startswith(("node_", "contig", "scaffold", "unitig")):
        return False
    # Assembly accessions
    if s.startswith(("GCA_", "GCF_")):
        return True
    # GenBank/RefSeq nucleotide accessions (rough)
    return bool(re.match(r"^[A-Z]{1,4}\d{5,}(?:\.\d+)?$", s))

# -----------------------------
# FASTA helpers
# -----------------------------

def prodigal_version() -> str:
    try:
        out = run_cmd([which_or_raise("prodigal"), "-v"])
        for ln in out.splitlines():
            ln = ln.strip()
            if ln:
                return ln
        return "unknown"
    except Exception:
        return "unknown"

def parse_fasta_headers(fasta_path: Path) -> Dict[str, str]:
    headers: Dict[str, str] = {}
    with fasta_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                h = line[1:].strip()
                acc = h.split()[0]
                if acc not in headers:
                    headers[acc] = h
    return headers

def wrap_fasta_seq(seq: str, width: int = 70) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))

def write_fasta(path: Path, records: List[Tuple[str, str]]) -> None:
    """Write FASTA atomically (best-effort)."""
    ensure_dir(path.parent)
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    with tmp.open("w", encoding="utf-8", newline="\n") as out:
        for h, s in records:
            out.write(f">{h}\n")
            out.write(wrap_fasta_seq(s) + "\n")
    os.replace(tmp, path)

def write_tsv(path: Path, rows: List[dict], header: List[str]) -> None:
    """Write TSV atomically (best-effort)."""
    ensure_dir(path.parent)
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    with tmp.open("w", encoding="utf-8", newline="\n") as out:
        w = csv.DictWriter(out, fieldnames=header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    os.replace(tmp, path)

def _sort_rows(rows: List[dict], key_fields: List[str]) -> List[dict]:
    def _k(r: dict):
        return tuple(str(r.get(k, "")) for k in key_fields)
    return sorted(rows, key=_k)



# -----------------------------
# Prodigal parsing
# -----------------------------

@dataclass(frozen=True)
class CDSKey:
    contig: str
    start_1: int
    end_1: int
    strand: str

@dataclass
class CDSInfo:
    key: CDSKey
    partial: str
    nt: str
    aa: str

PRODIGAL_HEADER_RE = re.compile(
    r"^(?P<contig>\S+?)_\d+\s+#\s+(?P<s>\d+)\s+#\s+(?P<e>\d+)\s+#\s+(?P<strandint>-?\d+)\s+#\s+(?P<attrs>.*)$"
)

def parse_prodigal_fasta(path: Path) -> Dict[CDSKey, Tuple[str, str]]:
    d: Dict[CDSKey, Tuple[str, str]] = {}
    header: Optional[str] = None
    chunks: List[str] = []

    def flush():
        nonlocal header, chunks
        if header is None:
            return
        m = PRODIGAL_HEADER_RE.match(header)
        if not m:
            header = None
            chunks = []
            return
        contig = m.group("contig")
        s = int(m.group("s"))
        e = int(m.group("e"))
        strandint = int(m.group("strandint"))
        strand = "+" if strandint == 1 else "-"
        attrs = m.group("attrs")
        pm = re.search(r"\bpartial=([01]{2})\b", attrs)
        partial = pm.group(1) if pm else ""
        key = CDSKey(contig=contig, start_1=min(s, e), end_1=max(s, e), strand=strand)
        seq = "".join(chunks).strip()
        if seq:
            d[key] = (seq, partial)
        header = None
        chunks = []

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        flush()

    return d

def build_cds_index_from_prodigal(faa_path: Path, fnn_path: Path) -> Dict[str, CDSInfo]:
    aa_map = parse_prodigal_fasta(faa_path)
    nt_map = parse_prodigal_fasta(fnn_path)

    index: Dict[str, CDSInfo] = {}
    for key, (nt_seq, partial) in nt_map.items():
        aa_seq, partial2 = aa_map.get(key, ("", partial))
        partial_flag = partial or partial2 or ""
        cds_id = f"{key.contig}|{key.start_1}|{key.end_1}|{key.strand}"
        index[cds_id] = CDSInfo(key=key, partial=partial_flag, nt=nt_seq, aa=aa_seq)
    return index

def prodigal_complete(partial_flag: str) -> bool:
    return partial_flag.strip() == "00"


# -----------------------------
# BLAST TSV parsing
# -----------------------------

DEFAULT_COLS = [
    "qaccver", "saccver", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

def sniff_blast_tsv_columns(tsv_path: Path) -> List[str]:
    with tsv_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("# Fields:"):
                raw = line.split(":", 1)[1].strip()
                parts = [p.strip() for p in raw.split(",")]
                mapped: List[str] = []
                for p in parts:
                    pl = p.lower()
                    if "query acc" in pl:
                        mapped.append("qaccver")
                    elif "subject acc" in pl:
                        mapped.append("saccver")
                    elif "% identity" in pl:
                        mapped.append("pident")
                    elif "alignment length" in pl:
                        mapped.append("length")
                    elif "mismatches" in pl:
                        mapped.append("mismatch")
                    elif "gap opens" in pl:
                        mapped.append("gapopen")
                    elif "q. start" in pl or "q start" in pl:
                        mapped.append("qstart")
                    elif "q. end" in pl or "q end" in pl:
                        mapped.append("qend")
                    elif "s. start" in pl or "s start" in pl:
                        mapped.append("sstart")
                    elif "s. end" in pl or "s end" in pl:
                        mapped.append("send")
                    elif "evalue" in pl:
                        mapped.append("evalue")
                    elif "bit score" in pl:
                        mapped.append("bitscore")
                    # NOTE: If your BLAST TSV contains qcovs/qlen, we detect it but we do NOT use it.
                    elif "qcovs" in pl:
                        mapped.append("qcovs")
                    elif "query length" in pl or "qlen" == pl:
                        mapped.append("qlen")
                    else:
                        mapped.append(re.sub(r"\s+", "_", pl))
                return mapped
            if not line.startswith("#") and line.strip():
                break
    return DEFAULT_COLS

@dataclass
class BlastHit:
    qacc: str
    sacc: str
    pident: float
    aln_len: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: Optional[int] = None
    qcovs: Optional[float] = None

def parse_blast_tsv_best_by_subject(tsv_path: Path, cols: List[str]) -> Dict[str, BlastHit]:
    col_idx = {c: i for i, c in enumerate(cols)}
    required = ["qaccver", "saccver", "pident", "length", "sstart", "send", "evalue", "bitscore"]
    for r in required:
        if r not in col_idx:
            raise RuntimeError(f"BLAST TSV missing required column '{r}'. Detected columns: {cols}")

    best: Dict[str, BlastHit] = {}
    with tsv_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(cols):
                parts = line.rstrip("\n").split()
            if len(parts) < len(cols):
                continue

            qacc = parts[col_idx["qaccver"]]
            sacc = parts[col_idx["saccver"]]
            try:
                pident = float(parts[col_idx["pident"]])
                aln_len = int(float(parts[col_idx["length"]]))
                sstart = int(float(parts[col_idx["sstart"]]))
                send = int(float(parts[col_idx["send"]]))
                evalue = float(parts[col_idx["evalue"]].replace("e", "E"))
                bitscore = float(parts[col_idx["bitscore"]])

                qlen = None
                if "qlen" in col_idx:
                    qlen = int(float(parts[col_idx["qlen"]]))

                qcovs = None
                if "qcovs" in col_idx:
                    qcovs = float(parts[col_idx["qcovs"]])
            except Exception:
                continue

            hit = BlastHit(
                qacc=qacc,
                sacc=sacc,
                pident=pident,
                aln_len=aln_len,
                sstart=sstart,
                send=send,
                evalue=evalue,
                bitscore=bitscore,
                qlen=qlen,
                qcovs=qcovs,
            )

            cur = best.get(sacc)
            if (cur is None) or (hit.bitscore > cur.bitscore) or (hit.bitscore == cur.bitscore and hit.aln_len > cur.aln_len):
                best[sacc] = hit

    return best




# -----------------------------
# BLAST subject -> contig mapping (robust to composite IDs)
# -----------------------------

_BLAST_ID_PREFIXES = {"lcl", "gnl", "ref", "gb", "emb", "dbj", "tpg", "tpe", "tpd", "sp", "tr"}

def _unique_preserve(xs: List[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for x in xs:
        if x and x not in seen:
            seen.add(x)
            out.append(x)
    return out

def possible_contig_ids(subject_id: str) -> List[str]:
    """
    Generate plausible contig identifiers from a BLAST subject id / header.

    Handles common patterns such as:
      - lcl|contig0001
      - ref|NC_000000.1| some description
      - contig0001|123|456
      - contig0001_12345:67890
      - IDs followed by whitespace description

    Returns candidates in best-first order; caller should pick the first present in cds_by_contig.
    """
    raw = (subject_id or "").strip()
    if not raw:
        return []
    tok = raw.split()[0]  # drop description
    cands = [raw, tok]

    # Split pipe-delimited tokens
    if "|" in tok:
        parts = [p for p in tok.split("|") if p]
        cands.extend(parts)

        # Common prefix|ID|... patterns: pick the element after prefix, and also last element
        if len(parts) >= 2 and parts[0].lower() in _BLAST_ID_PREFIXES:
            cands.append(parts[1])
        if parts:
            cands.append(parts[-1])

    # Strip common coordinate suffixes (very conservative)
    for base in list(cands):
        if not base:
            continue
        # contig:123-456 or contig_123:456
        base2 = re.split(r"[:]", base, maxsplit=1)[0]
        if base2 != base:
            cands.append(base2)
        base3 = re.split(r"[;]", base, maxsplit=1)[0]
        if base3 != base:
            cands.append(base3)

    return _unique_preserve(cands)

def map_subject_to_contig(subject_id: str, cds_by_contig: Dict[str, List[Tuple[str, 'CDSInfo']]]) -> Optional[str]:
    """
    Map a BLAST subject identifier to a contig key present in cds_by_contig.
    Returns the matched contig, or None if no mapping is possible.
    """
    if not subject_id:
        return None

    # Fast path: exact match
    if subject_id in cds_by_contig:
        return subject_id

    # Try candidate normalizations
    for cand in possible_contig_ids(subject_id):
        if cand in cds_by_contig:
            return cand

    # Last resort: try case-insensitive match (some tools normalize case)
    subj_lower = subject_id.lower()
    for contig in cds_by_contig.keys():
        if contig.lower() == subj_lower:
            return contig
    for cand in possible_contig_ids(subject_id):
        cand_lower = cand.lower()
        for contig in cds_by_contig.keys():
            if contig.lower() == cand_lower:
                return contig

    return None

# -----------------------------
# Matching BLAST hit -> CDS
# -----------------------------

def overlap_bp(a1: int, a2: int, b1: int, b2: int) -> int:
    lo = max(min(a1, a2), min(b1, b2))
    hi = min(max(a1, a2), max(b1, b2))
    return max(0, hi - lo + 1)

def build_cds_by_contig(index: Dict[str, CDSInfo]) -> Dict[str, List[Tuple[str, CDSInfo]]]:
    d: Dict[str, List[Tuple[str, CDSInfo]]] = {}
    for cds_id, info in index.items():
        d.setdefault(info.key.contig, []).append((cds_id, info))
    for contig in d:
        d[contig].sort(key=lambda x: (x[1].key.start_1, x[1].key.end_1, x[1].key.strand))
    return d

def best_cds_for_hit(cds_list: List[Tuple[str, CDSInfo]], sstart: int, send: int, strand: str) -> Tuple[Optional[Tuple[str, CDSInfo]], int]:
    bstart, bend = (sstart, send) if sstart <= send else (send, sstart)
    best_item: Optional[Tuple[str, CDSInfo]] = None
    best_ov = 0
    best_score = -1

    for cds_id, info in cds_list:
        ov = overlap_bp(info.key.start_1, info.key.end_1, bstart, bend)
        if ov <= 0:
            continue
        same = 1 if info.key.strand == strand else 0
        score = same * 1_000_000 + ov
        if score > best_score:
            best_score = score
            best_item = (cds_id, info)
            best_ov = ov

    return best_item, best_ov


# -----------------------------
# Cache layer + label
# -----------------------------

def cache_dir_for(fasta_hash: str, mode: str, prod_ver: str, cache_root: Path) -> Path:
    safe_ver = re.sub(r"[^A-Za-z0-9._-]+", "_", prod_ver)[:80]
    return cache_root / fasta_hash / f"prodigal_{mode}_{safe_ver}"

def build_or_load_cds_index(
    fasta_path: Path,
    mode: str,
    cache_root: Path,
    log,
    cache_label: str = "",
    cache_dir_override: Optional[Path] = None,
) -> Tuple[Path, Dict[str, CDSInfo]]:
    """Load an existing CDS index from cache, or build it with Prodigal.

    Behavior:
    - If `cache_dir_override` is provided and contains metadata.json + index.json, load it directly
      (NO FASTA hashing, NO Prodigal calls).
    - Otherwise, hash the FASTA to locate the cache namespace.
    - If any cache exists for this hash+mode, load it WITHOUT calling `prodigal -v`.
    - Only if no cache exists do we call Prodigal to build a new one.
    """

    def _load_index(cdir: Path) -> Dict[str, CDSInfo]:
        idx_path = cdir / "index.json"
        data = json.loads(idx_path.read_text(encoding="utf-8"))
        index: Dict[str, CDSInfo] = {}
        for cds_id, rec in data.items():
            key = CDSKey(contig=rec["contig"], start_1=rec["start_1"], end_1=rec["end_1"], strand=rec["strand"])
            index[cds_id] = CDSInfo(
                key=key,
                partial=rec.get("partial", ""),
                nt=rec.get("nt", ""),
                aa=rec.get("aa", ""),
            )
        return index

    def _maybe_update_label(cdir: Path, label: str) -> None:
        label = (label or "").strip()
        if not label:
            return
        meta_path = cdir / "metadata.json"
        if not meta_path.exists():
            return
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except Exception:
            meta = {}
        cur = (meta.get("cache_label") or "").strip()
        if cur != label:
            meta["cache_label"] = label
            safe_write_text(meta_path, json.dumps(meta, ensure_ascii=False, indent=2))

    # 1) User-selected cache override (fast path)
    if cache_dir_override is not None:
        cdir = Path(cache_dir_override)
        meta_path = cdir / "metadata.json"
        idx_path = cdir / "index.json"
        if meta_path.exists() and idx_path.exists():
            log(f"[CACHE] Using user-selected CDS index (no FASTA hashing): {cdir}")
            _maybe_update_label(cdir, cache_label)
            return cdir, _load_index(cdir)
        else:
            log(f"[CACHE] Selected cache invalid (missing metadata.json/index.json). Falling back to auto-cache lookup: {cdir}")

    # 2) Auto-cache lookup (requires hashing the FASTA)
    fasta_hash = sha256_file(fasta_path)
    hash_root = cache_root / fasta_hash

    if hash_root.exists() and hash_root.is_dir():
        # Load any existing cache for this hash+mode, WITHOUT calling Prodigal.
        candidates = sorted([p for p in hash_root.glob(f"prodigal_{mode}_*") if p.is_dir()])
        for cdir in candidates:
            meta_path = cdir / "metadata.json"
            idx_path = cdir / "index.json"
            if meta_path.exists() and idx_path.exists():
                log(f"[CACHE] Using existing CDS index: {cdir}")
                _maybe_update_label(cdir, cache_label)
                return cdir, _load_index(cdir)

    # 3) No cache found -> build new (ONLY place Prodigal is invoked)
    prod_ver = prodigal_version()
    cdir = cache_dir_for(fasta_hash, mode, prod_ver, cache_root)
    meta_path = cdir / "metadata.json"
    idx_path = cdir / "index.json"

    ensure_dir(cdir)
    log(f"[CACHE] No CDS index found. Building new index at: {cdir}")

    prod = which_or_raise("prodigal")
    gff = cdir / "prodigal.gff"
    faa = cdir / "prodigal.faa"
    fnn = cdir / "prodigal.fnn"

    cmd = [prod, "-i", str(fasta_path), "-f", "gff", "-o", str(gff), "-a", str(faa), "-d", str(fnn), "-p", mode, "-q"]
    out = run_cmd(cmd)
    if out.strip():
        for line in out.splitlines():
            log(line)

    log("[CACHE] Parsing Prodigal outputs into index...")
    raw_index = build_cds_index_from_prodigal(faa, fnn)

    serial: Dict[str, dict] = {}
    for cds_id, info in raw_index.items():
        serial[cds_id] = {
            "contig": info.key.contig,
            "start_1": info.key.start_1,
            "end_1": info.key.end_1,
            "strand": info.key.strand,
            "partial": info.partial,
            "nt": info.nt,
            "aa": info.aa,
        }

    safe_write_text(idx_path, json.dumps(serial, ensure_ascii=False))
    safe_write_text(meta_path, json.dumps({
        "fasta_path": str(fasta_path),
        "fasta_sha256": fasta_hash,
        "prodigal_version": prod_ver,
        "mode": mode,
        "cache_label": (cache_label or "").strip(),
    }, ensure_ascii=False, indent=2))

    log(f"[CACHE] Built CDS index: {len(raw_index)} CDS records.")
    return cdir, raw_index


# -----------------------------
# Analysis
# -----------------------------

def run_analysis(
    fasta_path: Path,
    blast_tsv_path: Path,
    output_dir: Path,
    prefix: str,
    prod_mode: str,
    cache_root: Path,
    min_cds_len_nt: int,     # <-- HARD MIN LENGTH FILTER
    cache_label: str,
    log,
    meta_mode: str,          # "none" | "ncbi" | "offline"
    ncbi_email: str,
    ncbi_api_key: str,
    max_cds_len_nt: int,     # <-- HARD MAX LENGTH FILTER (0 disables)
    min_qcov_pct: float,     # <-- MIN QUERY COVERAGE % (0 disables; requires qcovs or qlen)
    split_by_replicon: bool,  # <-- OPTIONAL: split PASS outputs by replicon type
    cache_dir_override: Optional[Path] = None,
    progress_cb=None,        # optional progress callback(current:int,total:int)
) -> Path:
    out_base = output_dir / f"{prefix}_OUTPUT"
    ensure_dir(out_base)

    filtered_base = out_base / "FILTERED"
    ensure_dir(filtered_base)

    fasta_headers = parse_fasta_headers(fasta_path)

    meta_mode = (meta_mode or "none").strip().lower()
    use_ncbi = meta_mode == "ncbi"
    use_offline = meta_mode == "offline"

    meta_resolver = MetaResolver(
        cache_path=default_meta_cache_path(),
        log=log,
        email=ncbi_email,
        api_key=ncbi_api_key,
        enabled=bool(use_ncbi),
    )

    if use_ncbi:
        log("[INFO] NCBI enrichment: contig IDs must be valid NCBI accessions for reliable lookups. "
            "Local assembly IDs (e.g., NODE_/scaffold_/contig_) typically will not resolve.")

    cache_dir, cds_index = build_or_load_cds_index(
        fasta_path=fasta_path,
        mode=prod_mode,
        cache_root=cache_root,
        log=log,
        cache_label=cache_label,
        cache_dir_override=cache_dir_override,
    )

        # Release guardrail: warn and stop if Prodigal produced no parseable CDS
    if not cds_index:
        raise RuntimeError(
            "No CDS features were parsed from Prodigal output. "
            "This usually indicates an incompatible/empty FASTA, a Prodigal failure, "
            "or unexpected FASTA headers/output format. "
            "Check the log above and try Prodigal mode 'meta' vs 'single'."
        )

    cds_by_contig = build_cds_by_contig(cds_index)

    cols = sniff_blast_tsv_columns(blast_tsv_path)
    hits = parse_blast_tsv_best_by_subject(blast_tsv_path, cols)
    log(f"[INFO] Parsed best HSPs for {len(hits)} subjects from BLAST TSV.")

    matched_nt: Dict[str, List[Tuple[str, str]]] = {}
    matched_aa: Dict[str, List[Tuple[str, str]]] = {}

    def norm_replicon(rt: str) -> str:
        rt = (rt or "").strip().lower()
        if rt in {"chromosome", "plasmid", "phage"}:
            return rt
        return "unknown"

    # Optional: outputs split by replicon type (PASS-only)
    matched_nt_rep: Dict[str, Dict[str, List[Tuple[str, str]]]] = {}
    matched_aa_rep: Dict[str, Dict[str, List[Tuple[str, str]]]] = {}
    summary_complete_rep: Dict[str, List[dict]] = {}
    summary_incomplete_rep: Dict[str, List[dict]] = {}

    summary_complete: List[dict] = []
    summary_incomplete: List[dict] = []
    filtered_complete: List[dict] = []
    filtered_incomplete: List[dict] = []

    log_rows: List[dict] = []

    ok = 0
    no_contig = 0
    no_cds = 0
    filtered_n = 0
    total = len(hits)
    logged_ncbi_skip_note = False

    if progress_cb:
        try:
            progress_cb(0, total)
        except Exception:
            pass

    for i, sacc in enumerate(sorted(hits.keys()), start=1):
        hit = hits[sacc]

        # Compute query-coverage once per hit (used in all branches)
        cov_pct = None
        if getattr(hit, "qcovs", None) is not None:
            try:
                cov_pct = float(hit.qcovs)
            except Exception:
                cov_pct = None
        elif getattr(hit, "qlen", None):
            try:
                qlen_val = int(hit.qlen)
                if qlen_val > 0:
                    cov_pct = (100.0 * float(hit.aln_len) / float(qlen_val))
            except Exception:
                cov_pct = None

        if min_qcov_pct > 0 and cov_pct is None:
            cov_filter = "NO_COV_INFO"
        elif min_qcov_pct <= 0:
            cov_filter = "NA"
        else:
            cov_filter = ("PASS" if cov_pct >= float(min_qcov_pct) else f"FILTERED_QCOV_LT_{min_qcov_pct}")

        passes_cov = (min_qcov_pct <= 0) or (cov_pct is not None and cov_pct >= float(min_qcov_pct))

        contig = map_subject_to_contig(sacc, cds_by_contig)
        if not contig:
            no_contig += 1

            # Coverage already computed for this hit (cov_pct/cov_filter)

            sstrand = "+" if hit.sstart <= hit.send else "-"
            log_rows.append({
                "subject": sacc,
                "contig": "",
                "qacc": hit.qacc,
                "pident": hit.pident,
                "aln_len": hit.aln_len,
                "cov_pct": ("" if cov_pct is None else round(float(cov_pct), 3)),
                "cov_filter": cov_filter,
                "blast_start": min(hit.sstart, hit.send),
                "blast_end": max(hit.sstart, hit.send),
                "blast_strand": sstrand,
                "cds_start": "",
                "cds_end": "",
                "cds_strand": "",
                "cds_len_nt": "",
                "orf_aa_len": "",
                "prodigal_partial": "",
                "status": "NO_CONTIG_MATCH",
                "overlap_bp": "",
                "bitscore": hit.bitscore,
                "evalue": hit.evalue,
                "fasta_header": sacc,
                "len_filter": "NA",
                "replicon_type": "unknown",
                "replicon_confidence": "low",
                "genus": "",
                "species": "",
                "taxid": "",
            })
            if progress_cb:
                try:
                    progress_cb(i, total)
                except Exception:
                    pass
            continue

        sstrand = "+" if hit.sstart <= hit.send else "-"
        cds_list = cds_by_contig.get(contig, [])

        best_item, ov = best_cds_for_hit(cds_list, hit.sstart, hit.send, sstrand)
        if best_item is None:
            no_cds += 1
            log_rows.append({
                "subject": sacc,
                "contig": contig,
                "qacc": hit.qacc,
                "pident": hit.pident,
                "aln_len": hit.aln_len,
                "cov_pct": ("" if cov_pct is None else round(float(cov_pct), 3)),
                "cov_filter": cov_filter,
                "blast_start": min(hit.sstart, hit.send),
                "blast_end": max(hit.sstart, hit.send),
                "blast_strand": sstrand,
                "cds_start": "",
                "cds_end": "",
                "cds_strand": "",
                "cds_len_nt": "",
                "orf_aa_len": "",
                "prodigal_partial": "",
                "status": "NO_CDS_MATCH",
                "overlap_bp": "",
                "bitscore": hit.bitscore,
                "evalue": hit.evalue,
                "fasta_header": fasta_headers.get(contig, sacc),
                "len_filter": "NA",
                "replicon_type": "unknown",
                "replicon_confidence": "low",
                "genus": "",
                "species": "",
                "taxid": "",
            })
            if progress_cb:
                try:
                    progress_cb(i, total)
                except Exception:
                    pass
            continue

        _, info = best_item
        complete = prodigal_complete(info.partial)
        status = "COMPLETE" if complete else "INCOMPLETE"
        ok += 1

        cds_len_nt_val = len(info.nt) if info.nt else (info.key.end_1 - info.key.start_1 + 1)
        orf_aa_len_val = cds_len_nt_val // 3

        passes_min = (min_cds_len_nt <= 0) or (cds_len_nt_val >= min_cds_len_nt)
        passes_max = (max_cds_len_nt <= 0) or (cds_len_nt_val <= max_cds_len_nt)
        passes_len = passes_min and passes_max

        orig_hdr = fasta_headers.get(contig, sacc)
        if use_offline:
            meta = offline_meta_from_header(fasta_headers.get(contig, ""))
        elif use_ncbi:
            if not is_accession_like(contig):
                if not logged_ncbi_skip_note:
                    log("[META] NCBI enrichment skipped for non-accession contig IDs (e.g., NODE_/scaffold). Use accession-like headers to enable.")
                    logged_ncbi_skip_note = True
                meta = meta_resolver.resolve("")  # returns empty/unknown structure
            else:
                meta = meta_resolver.resolve(contig)
        else:
            meta = meta_resolver.resolve("")  # disabled -> empty/unknown

        hdr = f"{orig_hdr} | CDS:{info.key.start_1}-{info.key.end_1}({info.key.strand}) partial={info.partial or 'NA'}"

        log_rows.append({
            "subject": sacc,
            "contig": contig,
            "qacc": hit.qacc,
            "pident": hit.pident,
            "aln_len": hit.aln_len,
            "cov_pct": ("" if cov_pct is None else round(float(cov_pct), 3)),
            "cov_filter": cov_filter,
            "blast_start": min(hit.sstart, hit.send),
            "blast_end": max(hit.sstart, hit.send),
            "blast_strand": sstrand,
            "cds_start": info.key.start_1,
            "cds_end": info.key.end_1,
            "cds_strand": info.key.strand,
            "cds_len_nt": cds_len_nt_val,
            "orf_aa_len": orf_aa_len_val,
            "prodigal_partial": info.partial,
            "status": status,
            "overlap_bp": ov,
            "bitscore": hit.bitscore,
            "evalue": hit.evalue,
            "fasta_header": fasta_headers.get(contig, sacc),
            "len_filter": (
                "PASS" if passes_len else
                (f"FILTERED_LT_{min_cds_len_nt}" if not passes_min else f"FILTERED_GT_{max_cds_len_nt}")
            ),
            "replicon_type": meta.get("replicon_type", "unknown"),
            "replicon_confidence": meta.get("replicon_confidence", "low"),
            "genus": meta.get("genus") or "",
            "species": meta.get("species") or "",
            "taxid": "" if meta.get("taxid") is None else meta.get("taxid"),
        })

        base_summary = {
            "accession": contig,
            "qacc": hit.qacc,
            "pident": hit.pident,
            "bitscore": hit.bitscore,
            "evalue": hit.evalue,
            "cov_pct": ("" if cov_pct is None else round(float(cov_pct), 3)),
            "cov_filter": cov_filter,
            "cds_start": info.key.start_1,
            "cds_end": info.key.end_1,
            "strand": info.key.strand,
            "cds_len_nt": cds_len_nt_val,
            "orf_aa_len": orf_aa_len_val,
            "prodigal_partial": info.partial,
            "fasta_header": fasta_headers.get(contig, sacc),
            "replicon_type": meta.get("replicon_type", "unknown"),
            "replicon_confidence": meta.get("replicon_confidence", "low"),
            "genus": meta.get("genus") or "",
            "species": meta.get("species") or "",
            "taxid": "" if meta.get("taxid") is None else meta.get("taxid"),
        }

        rep_type = norm_replicon(meta.get("replicon_type", "unknown"))

        if (not passes_len) or (not passes_cov):
            filtered_n += 1
            if status == "COMPLETE":
                filtered_complete.append(base_summary)
            else:
                filtered_incomplete.append(base_summary)
            if progress_cb:
                try:
                    progress_cb(i, total)
                except Exception:
                    pass
            continue

        if info.nt:
            matched_nt.setdefault(status, []).append((hdr, info.nt))
        if info.aa:
            matched_aa.setdefault(status, []).append((hdr, info.aa))

        if status == "COMPLETE":
            summary_complete.append(base_summary)
        else:
            summary_incomplete.append(base_summary)

        if split_by_replicon:
            if info.nt:
                matched_nt_rep.setdefault(rep_type, {}).setdefault(status, []).append((hdr, info.nt))
            if info.aa:
                matched_aa_rep.setdefault(rep_type, {}).setdefault(status, []).append((hdr, info.aa))

            if status == "COMPLETE":
                summary_complete_rep.setdefault(rep_type, []).append(base_summary)
            else:
                summary_incomplete_rep.setdefault(rep_type, []).append(base_summary)

        if progress_cb:
            try:
                progress_cb(i, total)
            except Exception:
                pass

    log_tsv = out_base / "CDS_MATCH_LOG.tsv"
    log_header = [
        "subject", "contig", "qacc", "pident", "aln_len", "cov_pct", "cov_filter",
        "blast_start", "blast_end", "blast_strand",
        "cds_start", "cds_end", "cds_strand",
        "cds_len_nt", "orf_aa_len",
        "prodigal_partial", "status",
        "overlap_bp", "bitscore", "evalue", "fasta_header",
        "len_filter",
        "replicon_type", "replicon_confidence", "genus", "species", "taxid",
    ]
    log_rows = _sort_rows(log_rows, ['subject','qacc','status'])
    write_tsv(log_tsv, log_rows, log_header)
    log(f"[INFO] Wrote: {log_tsv}")
    log(f"[INFO] OK: {ok}  NO_CONTIG_MATCH: {no_contig}  NO_CDS_MATCH: {no_cds}  FILTERED_LEN_OR_COV: {filtered_n}")

    if use_offline:
        log("[META] Offline header enrichment enabled (no network calls).")
    else:
        log(meta_resolver.summarize())

    summary_header = [
        "accession", "qacc", "pident", "bitscore", "evalue", "cov_pct", "cov_filter",
        "cds_start", "cds_end", "strand", "cds_len_nt", "orf_aa_len",
        "prodigal_partial", "fasta_header",
        "replicon_type", "replicon_confidence", "genus", "species", "taxid",
    ]
    summary_complete = _sort_rows(summary_complete, ['accession','qacc'])
    summary_incomplete = _sort_rows(summary_incomplete, ['accession','qacc'])
    write_tsv(out_base / "SUMMARY_COMPLETE.tsv", summary_complete, summary_header)
    write_tsv(out_base / "SUMMARY_INCOMPLETE.tsv", summary_incomplete, summary_header)
    log(f"[INFO] Wrote: {out_base/'SUMMARY_COMPLETE.tsv'} ({len(summary_complete)} rows)")
    log(f"[INFO] Wrote: {out_base/'SUMMARY_INCOMPLETE.tsv'} ({len(summary_incomplete)} rows)")

    filtered_complete = _sort_rows(filtered_complete, ['accession','qacc'])
    filtered_incomplete = _sort_rows(filtered_incomplete, ['accession','qacc'])
    write_tsv(filtered_base / "SUMMARY_COMPLETE.tsv", filtered_complete, summary_header)
    write_tsv(filtered_base / "SUMMARY_INCOMPLETE.tsv", filtered_incomplete, summary_header)
    log(f"[INFO] Wrote: {filtered_base/'SUMMARY_COMPLETE.tsv'} ({len(filtered_complete)} rows)")
    log(f"[INFO] Wrote: {filtered_base/'SUMMARY_INCOMPLETE.tsv'} ({len(filtered_incomplete)} rows)")

    for status, recs in matched_nt.items():
        out_dir2 = out_base / status
        write_fasta(out_dir2 / f"{prefix}_{status}_CDS_nt.fasta", recs)

    for status, recs in matched_aa.items():
        out_dir2 = out_base / status
        write_fasta(out_dir2 / f"{prefix}_{status}_CDS_aa.fasta", recs)

    # Optional: Replicon-split outputs (PASS-only)
    if split_by_replicon:
        split_root = out_base / "REPLICON_SPLIT"
        ensure_dir(split_root)

        for rep in ["chromosome", "plasmid", "phage", "unknown"]:
            rep_dir = split_root / rep
            ensure_dir(rep_dir)

            # TSV summaries per replicon
            write_tsv(rep_dir / "SUMMARY_COMPLETE.tsv", summary_complete_rep.get(rep, []), summary_header)
            write_tsv(rep_dir / "SUMMARY_INCOMPLETE.tsv", summary_incomplete_rep.get(rep, []), summary_header)

            # FASTA outputs per replicon + status
            for status, recs in matched_nt_rep.get(rep, {}).items():
                out_dir2 = rep_dir / status
                write_fasta(out_dir2 / f"{prefix}_{status}_CDS_nt.fasta", recs)

            for status, recs in matched_aa_rep.get(rep, {}).items():
                out_dir2 = rep_dir / status
                write_fasta(out_dir2 / f"{prefix}_{status}_CDS_aa.fasta", recs)

        log(f"[INFO] Replicon-split outputs written to: {split_root}")

    log(f"[INFO] Done. Outputs in: {out_base}")

    # Reproducibility: write run metadata
    run_info = {
        "app": __app_name__,
        "version": __version__,
        "build": __build__,
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "fasta_path": str(fasta_path),
        "blast_tsv_path": str(blast_tsv_path),
        "output_dir": str(output_dir),
        "prefix": prefix,
        "prodigal_mode": prod_mode,
        "cache_root": str(cache_root),
        "cache_dir": str(cache_dir),
        "min_cds_len_nt": min_cds_len_nt,
        "max_cds_len_nt": max_cds_len_nt,
        "min_qcov_pct": min_qcov_pct,
        "meta_mode": meta_mode,
        "split_by_replicon": bool(split_by_replicon),
    }
    try:
        safe_write_text(out_base / "RUN_INFO.json", json.dumps(run_info, indent=2))
        log(f"[INFO] Wrote: {out_base / 'RUN_INFO.json'}")
    except Exception as _e:
        log(f"[WARN] Could not write RUN_INFO.json: {_e}")


    # Minimal README to make the output self-describing (release polish)
    try:
        readme_txt = f"""{__app_name__} v{__version__}

What this run does:
  • Uses Prodigal to call CDS on the input genome FASTA (cached by FASTA SHA256).
  • Uses a BLAST outfmt6 TSV to locate your gene-of-interest per subject/contig and extracts the best-overlap CDS.
  • Applies hard CDS length filters (min/max, nt) and optional minimum query coverage (qcov).

Required inputs:
  1) Genome FASTA (complete genomes recommended).
  2) BLAST TSV (outfmt 6). If using qcov filtering, include qcovs or qlen in the outfmt fields.

Key outputs in this folder:
  • COMPLETE/ and INCOMPLETE/: PASS CDS FASTA outputs
  • FILTERED/: filtered-out sequences and summary TSVs
  • CDS_MATCH_LOG.tsv: per-subject matching log (contig mapping, overlap, filters)
  • RUN_INFO.json: run parameters + cache provenance

Notes:
  • NCBI enrichment only works reliably if contig IDs are valid NCBI accessions.
"""
        safe_write_text(out_base / "README.txt", readme_txt)
    except Exception as _e:
        log(f"[WARN] Could not write README.txt: {_e}")

    return out_base




# -----------------------------
# Background worker (prevents GUI freeze)
# -----------------------------

class AnalysisWorker(QtCore.QObject):
    log = QtCore.pyqtSignal(str)
    progress = QtCore.pyqtSignal(int, int)   # cur, total
    finished = QtCore.pyqtSignal(str)        # output dir as string
    error = QtCore.pyqtSignal(str)

    def __init__(
        self,
        *,
        fasta_path: Path,
        blast_tsv_path: Path,
        output_dir: Path,
        prefix: str,
        prod_mode: str,
        cache_root: Path,
        min_cds_len_nt: int,
        cache_label: str,
        meta_mode: str,
        ncbi_email: str,
        ncbi_api_key: str,
        max_cds_len_nt: int,
        min_qcov_pct: float,
        split_by_replicon: bool,
        cache_dir_override: Optional[Path],
    ):
        super().__init__()
        self.fasta_path = fasta_path
        self.blast_tsv_path = blast_tsv_path
        self.output_dir = output_dir
        self.prefix = prefix
        self.prod_mode = prod_mode
        self.cache_root = cache_root
        self.min_cds_len_nt = min_cds_len_nt
        self.cache_label = cache_label
        self.meta_mode = meta_mode
        self.ncbi_email = ncbi_email
        self.ncbi_api_key = ncbi_api_key
        self.max_cds_len_nt = max_cds_len_nt
        self.min_qcov_pct = float(min_qcov_pct)
        self.split_by_replicon = split_by_replicon
        self.cache_dir_override = cache_dir_override

    @QtCore.pyqtSlot()
    def run(self) -> None:
        try:
            def _log(msg: str) -> None:
                self.log.emit(str(msg))

            def _progress(cur: int, total: int) -> None:
                self.progress.emit(int(cur), int(total))

            out_base = run_analysis(
                fasta_path=self.fasta_path,
                blast_tsv_path=self.blast_tsv_path,
                output_dir=self.output_dir,
                prefix=self.prefix,
                prod_mode=self.prod_mode,
                cache_root=self.cache_root,
                min_cds_len_nt=self.min_cds_len_nt,
                cache_label=self.cache_label,
                log=_log,
                meta_mode=self.meta_mode,
                ncbi_email=self.ncbi_email,
                ncbi_api_key=self.ncbi_api_key,
                max_cds_len_nt=self.max_cds_len_nt,
                min_qcov_pct=self.min_qcov_pct,
                split_by_replicon=self.split_by_replicon,
                cache_dir_override=self.cache_dir_override,
                progress_cb=_progress,
            )
            self.finished.emit(str(out_base))
        except Exception as e:
            self.error.emit(str(e))


# -----------------------------
# GUI
# -----------------------------

class App(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"{__app_name__} v{__version__} — Blast Assisted Knowledge Organizer")
        self.resize(1080, 780)
        self.cache_root = default_cache_root()
        self._build_ui()
        self._refresh_cache_list()

    def log(self, msg: str) -> None:
        self.log_box.appendPlainText(msg)
        self.log_box.verticalScrollBar().setValue(self.log_box.verticalScrollBar().maximum())


    def _show_about(self) -> None:
        details = [
            f"{__app_name__} v{__version__}",
            "Blast-Assisted Knowledge Organizer",
            "",
            "What it does:",
            "• Builds a cached CDS index from FASTA inputs using Prodigal",
            "• Reconciles candidate CDSs using BLAST tabular (outfmt 6) hits",
            "• Applies length and coverage filters to retain full-length homologs",
            "• Produces PASS-only FASTA outputs plus FILTERED and TSV summary logs",
            "",
            "Dependencies:",
            "• Python 3.x, PyQt6",
            "• Prodigal (on PATH)",
            "",
            f"Build: {__build__}",
        ]
        if __repo__:
            details.append(f"Repo: {__repo__}")
        QtWidgets.QMessageBox.information(self, "About BAKO", "\n".join(details))




    def _set_running(self, running: bool) -> None:
        # Disable/enable UI controls during background work
        try:
            self.btn_run.setEnabled(not running)
        except Exception:
            pass
        try:
            self.btn_build.setEnabled(not running)
        except Exception:
            pass
        try:
            # Browsing is fine but can be confusing mid-run; disable for cleanliness
            self.fasta_edit.setEnabled(not running)
            self.blast_tsv_edit.setEnabled(not running)
            self.out_dir_edit.setEnabled(not running)
        except Exception:
            pass

    def _on_worker_log(self, msg: str) -> None:
        self.log(str(msg))

    def _on_worker_progress(self, cur: int, total: int) -> None:
        if not hasattr(self, "progress"):
            return
        try:
            self.progress.setMaximum(max(1, int(total)))
            self.progress.setValue(int(cur))
        except Exception:
            pass

    def _on_worker_finished(self, out_dir: str) -> None:
        self._set_running(False)
        self.log(f"[OK] Finished. Outputs in: {out_dir}")
        QtWidgets.QMessageBox.information(self, "BAKO", f"Finished. Outputs in:\n{out_dir}")

    def _on_worker_error(self, err: str) -> None:
        self._set_running(False)
        self.log(f"[ERROR] {err}")
        QtWidgets.QMessageBox.critical(self, "BAKO error", str(err))


    def _build_ui(self):
        layout = QtWidgets.QVBoxLayout(self)

        # Header
        header = QtWidgets.QHBoxLayout()
        lbl = QtWidgets.QLabel(f"{__app_name__} v{__version__}")
        fnt = lbl.font()
        fnt.setPointSize(max(10, fnt.pointSize()+2))
        fnt.setBold(True)
        lbl.setFont(fnt)
        header.addWidget(lbl)
        header.addStretch(1)
        btn_about = QtWidgets.QPushButton("About")
        btn_about.clicked.connect(self._show_about)
        header.addWidget(btn_about)
        layout.addLayout(header)

        ds = QtWidgets.QGroupBox("Dataset and CDS Index (cached Prodigal)")
        ds_layout = QtWidgets.QGridLayout()
        ds.setLayout(ds_layout)
        ds_layout.setColumnStretch(1, 1)

        self.fasta_edit = QtWidgets.QLineEdit()
        btn_fasta = QtWidgets.QPushButton("Browse…")
        btn_fasta.clicked.connect(self._browse_fasta)

        self.mode_combo = QtWidgets.QComboBox()
        self.mode_combo.addItems(["single", "meta"])

        self.cache_label_edit = QtWidgets.QLineEdit()
        self.cache_label_edit.setPlaceholderText("Optional cache name (e.g., 'msrA Dec15 meta')")

        self.cache_root_edit = QtWidgets.QLineEdit(str(self.cache_root))
        btn_cache = QtWidgets.QPushButton("Change…")
        btn_cache.clicked.connect(self._browse_cache_root)

        self.cache_combo = QtWidgets.QComboBox()
        btn_refresh = QtWidgets.QPushButton("Refresh list")
        btn_refresh.clicked.connect(self._refresh_cache_list)
        btn_delete_cache = QtWidgets.QPushButton("Delete selected cache")
        btn_delete_cache.clicked.connect(self._delete_selected_cache)

        self.btn_build = QtWidgets.QPushButton("Build CDS Index (Prodigal)")
        self.btn_build.clicked.connect(self._build_index_clicked)

        ds_layout.addWidget(QtWidgets.QLabel("FASTA dataset:"), 0, 0)
        ds_layout.addWidget(self.fasta_edit, 0, 1)
        ds_layout.addWidget(btn_fasta, 0, 2)

        ds_layout.addWidget(QtWidgets.QLabel("Prodigal mode:"), 1, 0)
        ds_layout.addWidget(self.mode_combo, 1, 1)

        ds_layout.addWidget(QtWidgets.QLabel("Cache name (label):"), 2, 0)
        ds_layout.addWidget(self.cache_label_edit, 2, 1)

        ds_layout.addWidget(QtWidgets.QLabel("Cache root:"), 3, 0)
        ds_layout.addWidget(self.cache_root_edit, 3, 1)
        ds_layout.addWidget(btn_cache, 3, 2)

        ds_layout.addWidget(QtWidgets.QLabel("Existing CDS indexes:"), 4, 0)
        ds_layout.addWidget(self.cache_combo, 4, 1)
        ds_layout.addWidget(btn_refresh, 4, 2)
        ds_layout.addWidget(btn_delete_cache, 4, 3)
        ds_layout.setColumnStretch(3, 0)

        ds_layout.addWidget(self.btn_build, 5, 0, 1, 3)

        layout.addWidget(ds)

        an = QtWidgets.QGroupBox("Analysis (PASS-only FASTA outputs + FILTERED folder)")
        an_layout = QtWidgets.QGridLayout()
        an.setLayout(an_layout)
        an_layout.setColumnStretch(1, 1)

        self.blast_tsv_edit = QtWidgets.QLineEdit()
        btn_tsv = QtWidgets.QPushButton("Browse…")
        btn_tsv.clicked.connect(self._browse_blast_tsv)

        self.out_dir_edit = QtWidgets.QLineEdit()
        btn_out = QtWidgets.QPushButton("Browse…")
        btn_out.clicked.connect(self._browse_out_dir)

        self.prefix_edit = QtWidgets.QLineEdit("msrA")

        self.min_len_nt = QtWidgets.QSpinBox()
        self.min_len_nt.setRange(0, 1000000)
        self.min_len_nt.setValue(0)
        self.min_len_nt.setToolTip("Minimum CDS length (nt). 0 disables filtering. Below-min goes to OUTPUT/FILTERED.")

        self.max_len_nt = QtWidgets.QSpinBox()
        self.max_len_nt.setRange(0, 1000000)
        self.max_len_nt.setValue(0)
        self.max_len_nt.setToolTip("Maximum CDS length (nt). 0 disables filtering. Above-max goes to OUTPUT/FILTERED.")


        self.min_qcov = QtWidgets.QDoubleSpinBox()
        self.min_qcov.setRange(0.0, 100.0)
        self.min_qcov.setDecimals(2)
        self.min_qcov.setSingleStep(1.0)
        self.min_qcov.setValue(0.0)
        self.min_qcov.setToolTip("Minimum query coverage (%) for BLAST hits. 0 disables. Requires TSV to include qcovs or qlen.")

        self.progress = QtWidgets.QProgressBar()
        self.progress.setMinimum(0)
        self.progress.setValue(0)
        self.progress.setFormat("%p%")

        # Metadata enrichment mode (mutually exclusive)
        self.rb_meta_none = QtWidgets.QRadioButton("No metadata enrichment")
        self.rb_meta_none.setChecked(True)
        self.rb_meta_ncbi = QtWidgets.QRadioButton("Enrich metadata (NCBI): replicon + genus/species")
        self.rb_meta_offline = QtWidgets.QRadioButton("Offline metadata from headers (faster, less reliable)")

        self.meta_group = QtWidgets.QButtonGroup(self)
        self.meta_group.setExclusive(True)
        self.meta_group.addButton(self.rb_meta_none)
        self.meta_group.addButton(self.rb_meta_ncbi)
        self.meta_group.addButton(self.rb_meta_offline)

        self.ncbi_email_edit = QtWidgets.QLineEdit()
        self.ncbi_email_edit.setPlaceholderText("NCBI email (recommended, optional)")

        self.ncbi_key_edit = QtWidgets.QLineEdit()
        self.ncbi_key_edit.setPlaceholderText("NCBI API key (optional)")

        self.btn_run = QtWidgets.QPushButton("Run Phase 1")
        self.btn_run.clicked.connect(self._run_clicked)

        self.chk_split_rep = QtWidgets.QCheckBox("Split PASS outputs by replicon type (chromosome/plasmid/phage/unknown)")
        self.chk_split_rep.setChecked(False)

        an_layout.addWidget(QtWidgets.QLabel("BLAST TSV (outfmt 6):"), 0, 0)
        an_layout.addWidget(self.blast_tsv_edit, 0, 1)
        an_layout.addWidget(btn_tsv, 0, 2)

        an_layout.addWidget(QtWidgets.QLabel("Output directory:"), 1, 0)
        an_layout.addWidget(self.out_dir_edit, 1, 1)
        an_layout.addWidget(btn_out, 1, 2)

        an_layout.addWidget(QtWidgets.QLabel("Prefix:"), 2, 0)
        an_layout.addWidget(self.prefix_edit, 2, 1)

        an_layout.addWidget(QtWidgets.QLabel("Min CDS len (nt):"), 3, 0)
        an_layout.addWidget(self.min_len_nt, 3, 1)

        an_layout.addWidget(QtWidgets.QLabel("Max CDS len (nt):"), 4, 0)
        an_layout.addWidget(self.max_len_nt, 4, 1)

        an_layout.addWidget(QtWidgets.QLabel("Progress:"), 4, 2)
        an_layout.addWidget(self.progress, 4, 3)

        an_layout.addWidget(QtWidgets.QLabel("Min query cov (%):"), 5, 0)
        an_layout.addWidget(self.min_qcov, 5, 1)

        an_layout.addWidget(self.rb_meta_none, 6, 0, 1, 4)
        an_layout.addWidget(self.rb_meta_ncbi, 7, 0, 1, 4)
        an_layout.addWidget(self.rb_meta_offline, 8, 0, 1, 4)

        an_layout.addWidget(QtWidgets.QLabel("NCBI email:"), 9, 0)
        an_layout.addWidget(self.ncbi_email_edit, 9, 1, 1, 3)
        an_layout.addWidget(QtWidgets.QLabel("NCBI API key:"), 10, 0)
        an_layout.addWidget(self.ncbi_key_edit, 10, 1, 1, 3)

        an_layout.addWidget(self.btn_run, 11, 0, 1, 4)
        an_layout.addWidget(self.chk_split_rep, 12, 0, 1, 4)

        layout.addWidget(an)

        self.log_box = QtWidgets.QPlainTextEdit()
        self.log_box.setReadOnly(True)
        layout.addWidget(self.log_box, stretch=1)

    def _delete_selected_cache(self):
        try:
            pdir = self.cache_combo.currentData()
        except Exception:
            pdir = None
        if not pdir:
            self.log("[CACHE] No cache selected.")
            return

        cache_path = Path(str(pdir))
        if not cache_path.exists() or not cache_path.is_dir():
            self.log(f"[CACHE] Selected cache does not exist: {cache_path}")
            self._refresh_cache_list()
            return

        reply = QtWidgets.QMessageBox.question(
            self,
            "Delete cache",
            f"Delete this cached CDS index?\n\n{cache_path}\n\nThis cannot be undone.",
            QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No,
            QtWidgets.QMessageBox.StandardButton.No,
        )
        if reply != QtWidgets.QMessageBox.StandardButton.Yes:
            self.log("[CACHE] Deletion cancelled.")
            return

        try:
            shutil.rmtree(cache_path)
            self.log(f"[CACHE] Deleted: {cache_path}")
        except Exception as e:
            self.log(f"[CACHE][ERROR] Failed to delete cache: {e}")

        self._refresh_cache_list()

    def _browse_fasta(self):
        p, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select FASTA", "", "FASTA (*.fa *.fasta *.fna *.ffn *.txt);;All files (*)")
        if p:
            self.fasta_edit.setText(p)

    def _browse_blast_tsv(self):
        p, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select BLAST TSV", "", "TSV (*.tsv *.txt);;All files (*)")
        if p:
            self.blast_tsv_edit.setText(p)

    def _browse_out_dir(self):
        p = QtWidgets.QFileDialog.getExistingDirectory(self, "Select output directory", "")
        if p:
            self.out_dir_edit.setText(p)

    def _browse_cache_root(self):
        p = QtWidgets.QFileDialog.getExistingDirectory(self, "Select cache root directory", "")
        if p:
            self.cache_root_edit.setText(p)
            self.cache_root = Path(p)
            self._refresh_cache_list()

    def _refresh_cache_list(self):
        root = Path(self.cache_root_edit.text().strip() or str(default_cache_root()))
        self.cache_root = root
        ensure_dir(root)

        items: List[Tuple[str, str]] = []
        for hdir in sorted(root.glob("*")):
            if not hdir.is_dir():
                continue
            for pdir in sorted(hdir.glob("prodigal_*")):
                meta = pdir / "metadata.json"
                if not meta.exists():
                    continue
                try:
                    m = json.loads(meta.read_text(encoding="utf-8"))
                    user_label = (m.get("cache_label") or "").strip()
                    if user_label:
                        label = f"{user_label}  |  {hdir.name[:12]}...  mode={m.get('mode','?')}  ver={str(m.get('prodigal_version','?'))[:18]}"
                    else:
                        label = f"{hdir.name[:12]}...  mode={m.get('mode','?')}  ver={str(m.get('prodigal_version','?'))[:18]}"
                    items.append((str(pdir), label))
                except Exception:
                    items.append((str(pdir), f"{hdir.name[:12]}...  {pdir.name}"))

        self.cache_combo.clear()
        for pdir, label in items:
            self.cache_combo.addItem(label, userData=pdir)

        self.log(f"[CACHE] Found {len(items)} cached CDS indexes under {root}")

    def _build_index_clicked(self):
        self.log_box.clear()
        fasta = Path(self.fasta_edit.text().strip())
        if not fasta.exists():
            self.log("[ERROR] Select a valid FASTA dataset.")
            return

        mode = self.mode_combo.currentText().strip()
        label = (self.cache_label_edit.text() or "").strip()
        self.cache_root = Path(self.cache_root_edit.text().strip() or str(default_cache_root()))

        try:
            _ = which_or_raise("prodigal")
            out_dir, idx = build_or_load_cds_index(fasta, mode, self.cache_root, self.log, cache_label=label, cache_dir_override=None)
            self.log(f"[OK] CDS index ready: {out_dir} ({len(idx)} CDS)")
            self._refresh_cache_list()
        except Exception as e:
            self.log(f"[ERROR] {e}")

    def _run_clicked(self):
        self.log_box.clear()

        fasta = Path(self.fasta_edit.text().strip())
        tsv = Path(self.blast_tsv_edit.text().strip())
        out_dir = Path(self.out_dir_edit.text().strip() or str(Path.cwd()))
        prefix = self.prefix_edit.text().strip() or "project"
        mode = self.mode_combo.currentText().strip()
        label = (self.cache_label_edit.text() or "").strip()

        if hasattr(self, "progress"):
            self.progress.setValue(0)
            self.progress.setMaximum(1)

        self.cache_root = Path(self.cache_root_edit.text().strip() or str(default_cache_root()))

        # If the user selected an existing cache, use it directly (fast; skips FASTA hashing).
        selected_cache_dir = None
        try:
            pdir = self.cache_combo.currentData()
        except Exception:
            pdir = None
        if pdir:
            selected_cache_dir = Path(str(pdir))

        # Determine metadata enrichment mode (mutually exclusive)
        meta_mode = "none"
        if hasattr(self, "rb_meta_ncbi") and self.rb_meta_ncbi.isChecked():
            meta_mode = "ncbi"
        elif hasattr(self, "rb_meta_offline") and self.rb_meta_offline.isChecked():
            meta_mode = "offline"

        ncbi_email = (self.ncbi_email_edit.text().strip() if hasattr(self, "ncbi_email_edit") else "")
        ncbi_api_key = (self.ncbi_key_edit.text().strip() if hasattr(self, "ncbi_key_edit") else "")
        max_cds_len_nt = int(self.max_len_nt.value()) if hasattr(self, "max_len_nt") else 0
        min_qcov_pct = float(self.min_qcov.value()) if hasattr(self, "min_qcov") else 0.0
        split_by_replicon = self.chk_split_rep.isChecked() if hasattr(self, "chk_split_rep") else False

        if not fasta.exists():
            self.log("[ERROR] Select a valid FASTA dataset.")
            return
        if not tsv.exists():
            self.log("[ERROR] Select a valid BLAST TSV.")
            return

        # Start background worker (prevents GUI freeze)
        self._set_running(True)

        self._thread = QtCore.QThread(self)
        self._worker = AnalysisWorker(
            fasta_path=fasta,
            blast_tsv_path=tsv,
            output_dir=out_dir,
            prefix=prefix,
            prod_mode=mode,
            cache_root=self.cache_root,
            min_cds_len_nt=int(self.min_len_nt.value()),
            cache_label=label,
            meta_mode=meta_mode,
            ncbi_email=ncbi_email,
            ncbi_api_key=ncbi_api_key,
            max_cds_len_nt=max_cds_len_nt,
            min_qcov_pct=min_qcov_pct,
            split_by_replicon=split_by_replicon,
            cache_dir_override=selected_cache_dir,
        )
        self._worker.moveToThread(self._thread)

        self._thread.started.connect(self._worker.run)
        self._worker.log.connect(self._on_worker_log)
        self._worker.progress.connect(self._on_worker_progress)
        self._worker.finished.connect(self._on_worker_finished)
        self._worker.error.connect(self._on_worker_error)

        # Always stop thread at end
        self._worker.finished.connect(self._thread.quit)
        self._worker.error.connect(self._thread.quit)
        self._thread.finished.connect(self._worker.deleteLater)
        self._thread.finished.connect(self._thread.deleteLater)

        self._thread.start()



def main():
    import sys
    app = QtWidgets.QApplication(sys.argv)
    w = App()
    w.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
