BAKO: Blast-Assisted Knowledge Organizer for Targeted CDS Extraction

BAKO (Blast-Assisted Knowledge Organizer) is a standalone desktop application designed for the reproducible extraction of coding sequences (CDSs) corresponding to a gene of interest from bacterial genome assemblies. The tool integrates external sequence similarity results with de novo gene prediction to enable precise localization and standardized retrieval of homologous loci across genomes.

Overview and Design Principles

BAKO operates under a modular design in which sequence similarity search is decoupled from downstream genomic interpretation. Specifically, the tool does not perform sequence alignment internally but instead consumes precomputed BLAST tabular outputs (outfmt 6). This separation ensures flexibility in upstream analysis while maintaining strict control over downstream processing and filtering.

The core functionality of BAKO consists of three stages: (i) CDS prediction, (ii) mapping of BLAST hits to genomic features, and (iii) application of explicit filtering and classification criteria. All steps are deterministic and logged to ensure reproducibility.

CDS Prediction and Caching

Coding sequences are predicted from input genome FASTA files using Prodigal (Hyatt et al., 2010), with support for both single-genome and metagenomic modes. For each predicted CDS, genomic coordinates, strand orientation, nucleotide sequence, amino acid sequence, and Prodigal completeness flags are retained.

To avoid redundant computation, CDS predictions are cached using a SHA256 hash of the input FASTA file. Cached results are reused across runs, ensuring both computational efficiency and reproducibility. The cache additionally records Prodigal version and execution parameters to maintain provenance.

BLAST Parsing and Hit Selection

BAKO parses BLAST tabular outputs (outfmt 6) using a robust column-detection strategy that supports both standard and annotated BLAST outputs. For each subject sequence (i.e., genome contig), only the highest-scoring alignment is retained, prioritizing bit score and alignment length.

Subject identifiers are normalized to ensure compatibility between BLAST outputs and genome FASTA headers. This includes handling of common identifier formats (e.g., lcl|, ref|, and composite identifiers) and removal of coordinate suffixes when present.

Mapping of BLAST Hits to CDS Features

BLAST hits are mapped onto predicted CDSs based on genomic coordinates. For each alignment, all overlapping CDSs are evaluated, and the CDS with the maximum overlap is selected. Strand concordance between the BLAST hit and CDS is strongly prioritized to ensure biologically consistent assignments.

If no overlapping CDS is identified, the hit is classified accordingly and retained in diagnostic logs.

Filtering and Classification

BAKO applies explicit and user-defined filtering criteria to ensure high-confidence sequence extraction. These include:

CDS length constraints (minimum and optional maximum, in nucleotides)
Query coverage thresholds, derived either from BLAST-reported coverage (qcovs) or computed from alignment length and query length (qlen)

Each CDS is further classified based on Prodigal completeness flags:

COMPLETE: full-length CDS (partial=00)
INCOMPLETE: truncated or partial CDS

Filtered sequences are excluded from primary outputs but retained in separate summaries to ensure transparency.

Metadata Enrichment (Optional)

BAKO supports optional enrichment of sequence metadata through two complementary approaches:

NCBI-based enrichment, using Entrez E-utilities to retrieve organismal and replicon information (e.g., chromosome, plasmid, phage), with built-in caching and rate-limit handling.
Header-based inference, which extracts taxonomic and replicon information directly from FASTA headers.

These annotations enable downstream stratification of results and facilitate comparative analyses across genomic contexts.

Output Structure and Reproducibility

BAKO produces a structured set of outputs organized into distinct categories:

FASTA files of extracted CDSs (nucleotide and amino acid sequences)
Tabular summaries for complete and incomplete CDSs
Comprehensive per-hit logs detailing mapping decisions, overlap metrics, and filtering outcomes
Filtered result summaries
A machine-readable metadata file (RUN_INFO.json) containing all run parameters and cache provenance

Optionally, outputs can be further partitioned by replicon type (e.g., chromosome, plasmid, phage), enabling stratified analyses.

All outputs are deterministically ordered and generated using atomic write operations to ensure data integrity.

GUI

<img width="1078" height="815" alt="image" src="https://github.com/user-attachments/assets/e4096545-ee45-4331-9534-1c96b749d753" />


Implementation

BAKO is implemented in Python (version ≥3.9) with a PyQt6-based graphical user interface for cross-platform usability. External dependencies are limited to Prodigal for gene prediction. The application is designed to operate locally without reliance on cloud services or external pipelines.

Scope and Limitations

BAKO is intended for targeted gene-of-interest extraction and does not perform genome assembly, annotation, or sequence alignment. Its performance depends on the quality of external BLAST results and the consistency between BLAST subject identifiers and genome FASTA headers.
