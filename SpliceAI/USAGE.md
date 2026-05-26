# Usage

## Command

```bash
python SpliceAI/spliceai_sequence_scanner.py \
  --input INPUT \
  --output_prefix OUTPUT_PREFIX \
  --format {genbank,fasta,auto}
```

Or run it in Docker with the bundled image:

```bash
docker run --rm \
  -v "$PWD":/work \
  -w /work \
  spliceai-scanner \
  --input /work/test_data/pARK.gb \
  --output_prefix /work/test_outputs/pARK \
  --format genbank
```

## Required arguments

- `--input`: input sequence file in FASTA or GenBank format
- `--output_prefix`: prefix used for all output files
- `--format`: one of `genbank`, `fasta`, or `auto`

## Optional arguments

- `--model_dir`: directory containing `spliceai1.h5` through `spliceai5.h5`
- `--seq_id`: record ID to use when the input contains multiple records
- `--min_score`: minimum score threshold for the filtered TSV
- `--plot_width`: plot width in inches
- `--plot_height`: plot height in inches
- `--title`: custom plot title
- `--dpi`: output plot resolution
- `--force`: overwrite existing outputs

## Output behavior

- The script creates the output directory if needed.
- Existing output files are not overwritten unless `--force` is supplied.
- Logging is written to `{output_prefix}.log`.

## Notes

- Coordinates in the TSV files are 1-based.
- The filtered TSV keeps rows where either donor or acceptor score meets
  `--min_score`.
- If `--model_dir` is omitted, the script uses `SpliceAI/splice_models`.
- Relative `--output_prefix` values are resolved against the script directory.

## Examples

FASTA input:

```bash
python SpliceAI/spliceai_sequence_scanner.py \
  --input example.fa \
  --output_prefix results/example \
  --format fasta
```

GenBank input with a specific record:

```bash
python SpliceAI/spliceai_sequence_scanner.py \
  --input example.gbk \
  --output_prefix results/example \
  --format genbank \
  --seq_id record_name \
  --min_score 0.2 \
  --force
```
