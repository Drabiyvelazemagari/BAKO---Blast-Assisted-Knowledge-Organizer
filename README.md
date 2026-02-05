# BAKO v2.0.0 — Blast-Assisted Knowledge Organizer

**Release file:** `BAKO_V2.py`  
**Interface:** PyQt6 desktop GUI  
**Scope:** Gene-of-interest extraction from genome FASTA using external BLAST TSV results

---

## What BAKO is

BAKO is a **standalone desktop application** that organizes and extracts the **best matching coding sequence (CDS)** for a gene of interest from bacterial genomes.

It is designed to:
- take a **genome FASTA**
- take a **BLAST tabular output (outfmt 6) produced elsewhere**
- predict CDSs using **Prodigal**
- map BLAST hits onto predicted CDSs
- apply strict, explicit filtering
- output reproducible FASTA and TSV summaries

BAKO is **not a BLAST pipeline** and **does not run BLAST**.  
It consumes BLAST results you already generated.

---

## What BAKO does (exactly)

For a single run, BAKO performs the following steps:

1. **Parse genome FASTA**
   - Reads contig headers and sequences
   - Optionally parses metadata directly from FASTA headers

2. **Predict CDS features using Prodigal**
   - User-selectable mode:
     - `single` (complete genomes)
     - `meta` (fragmented assemblies)
   - Extracts:
     - CDS coordinates
     - strand
     - nucleotide sequence
     - amino-acid sequence
     - Prodigal completeness flag (`partial=`)

3. **Cache CDS predictions**
   - CDS results are cached using a SHA256 hash of the FASTA
   - Prevents unnecessary re-running of Prodigal
   - Cache is persistent across runs

4. **Parse BLAST TSV (outfmt 6)**
   - Robust tab-delimited parsing
   - Handles Windows CRLF safely
   - Retains the **best hit per subject** (by bitscore)

5. **Map BLAST hits to CDSs**
   - BLAST subject IDs are normalized to match FASTA/Prodigal contigs
   - For each hit:
     - CDS overlap is computed
     - strand consistency is prioritized
     - the CDS with the **maximum overlap** is selected

6. **Apply filters**
   - CDS length (nt):
     - minimum
     - maximum
   - Query coverage (%), optional
   - Filters are explicit and logged

7. **Classify CDSs**
   - `COMPLETE` → Prodigal `partial=00`
   - `INCOMPLETE` → any other `partial` flag
   - Filtered hits are separated but retained in logs

8. **Write structured outputs**
   - FASTA files (nt + aa)
   - TSV summaries
   - Full decision log
   - Machine-readable run metadata

---

## What BAKO does NOT do

- Does **not** run BLAST
- Does **not** assemble genomes
- Does **not** annotate entire genomes
- Does **not** upload data anywhere
- Does **not** require BLAST executables to be installed

---

## Required inputs

### 1. Genome FASTA
Accepted extensions:
- `.fa`, `.fasta`, `.fna`

Requirements:
- Valid FASTA format
- Contig IDs must match those used when generating the BLAST TSV

### 2. BLAST TSV (outfmt 6)
Produced externally by the user.

Required fields:
- `qaccver`
- `saccver`
- `pident`
- `length`
- `sstart`
- `send`
- `evalue`
- `bitscore`

Optional fields:
- `qcovs`
- `qlen`

If `qcovs` or `qlen` are absent, query coverage filtering cannot be applied.

---

## Output structure

```
<PREFIX>_OUTPUT/
├── COMPLETE/
├── INCOMPLETE/
├── FILTERED/
├── SUMMARY_COMPLETE.tsv
├── SUMMARY_INCOMPLETE.tsv
├── CDS_MATCH_LOG.tsv
├── RUN_INFO.json
└── README.txt
```

---

## Installation requirements

### Python
- Python 3.9+

### Python packages
- PyQt6

```bash
pip install PyQt6
```

### External software
- **Prodigal** must be installed and available on PATH

---

## Running BAKO

```bash
python BAKO_V2.py
```

---

## Version information

- Application: BAKO
- Version: 2.0.0
- Build: PyQt6 cross-platform
