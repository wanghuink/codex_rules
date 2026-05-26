# SpliceAI Sequence Scanner

This directory contains a command-line SpliceAI scanner, the bundled model
files, and a Docker image definition for reproducible runs.

## Contents

- `spliceai_sequence_scanner.py`: scans a single FASTA or GenBank record and
  writes per-position donor and acceptor scores
- `splice_models/`: the five official SpliceAI model files
- `Dockerfile`: container image definition for running the scanner
- `USAGE.md`: command-line examples and option reference

## What the scanner writes

For a given `--output_prefix`, the script writes:

- `{output_prefix}.spliceai_scores.tsv`
- `{output_prefix}.spliceai_scores.filtered.tsv`
- `{output_prefix}.spliceai_scores.png`
- `{output_prefix}.log`

Positions in the TSV files are 1-based relative to the selected input record.
The filtered table keeps rows where `max(donor_score, acceptor_score)` is at
least `--min_score`.

## Build

Build the Docker image from the repository root:

```bash
docker build -t spliceai-scanner -f SpliceAI/Dockerfile .
```

## Example

Run the scanner against one of the bundled test GenBank files:

```bash
docker run --rm \
  -v "$PWD":/work \
  -w /work \
  spliceai-scanner \
  --input /work/test_data/pARK.gb \
  --output_prefix /work/test_outputs/pARK \
  --format genbank
```

The script defaults to `SpliceAI/splice_models` when `--model_dir` is omitted.
Relative `--output_prefix` values are resolved against the script directory.
