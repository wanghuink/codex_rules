#!/usr/bin/env python3
"""Scan a single FASTA/GenBank sequence with SpliceAI.

This script reads one sequence record from a FASTA or GenBank file, pads it
with neutral context, runs the five official SpliceAI models from a user
supplied directory, and writes:

* a full per-position TSV
* a filtered TSV using a configurable minimum score threshold
* a donor/acceptor line plot
* a log file in the output directory

Coordinates in the TSV are 1-based positions relative to the selected input
sequence record.
"""

from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Tuple


@dataclass(frozen=True)
class Record:
    record_id: str
    description: str
    sequence: str


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run SpliceAI on a single sequence record from a FASTA or GenBank "
            "file and report donor/acceptor scores for every position."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", required=True, help="Input sequence file.")
    parser.add_argument(
        "--output_prefix",
        required=True,
        help="Prefix for all output files. The parent directory is created if needed.",
    )
    parser.add_argument(
        "--model_dir",
        help=(
            "Directory containing the five SpliceAI model files. If omitted, "
            "the script uses splice_models next to this script."
        ),
    )
    parser.add_argument(
        "--format",
        required=True,
        choices=("genbank", "fasta", "auto"),
        help="Input format or auto-detection mode.",
    )
    parser.add_argument(
        "--seq_id",
        help=(
            "Sequence record ID to use if the input contains multiple records. "
            "If omitted and multiple records exist, the first record is used."
        ),
    )
    parser.add_argument(
        "--min_score",
        type=float,
        default=0.0,
        help=(
            "Minimum score threshold for the filtered TSV. A row is kept when "
            "either donor or acceptor score is at least this value."
        ),
    )
    parser.add_argument(
        "--plot_width",
        type=float,
        default=12.0,
        help="Plot width in inches.",
    )
    parser.add_argument(
        "--plot_height",
        type=float,
        default=4.0,
        help="Plot height in inches.",
    )
    parser.add_argument(
        "--title",
        help="Optional plot title. Defaults to the selected record ID.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Plot resolution in dots per inch.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files if they already exist.",
    )
    return parser


def configure_logging(log_path: Path) -> logging.Logger:
    logger = logging.getLogger("spliceai_sequence_scanner")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    formatter = logging.Formatter(
        fmt="%(asctime)s\t%(levelname)s\t%(message)s",
        datefmt="%Y-%m-%dT%H:%M:%SZ",
    )

    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    logger.propagate = False
    return logger


def detect_format(input_path: Path) -> str:
    suffix = input_path.suffix.lower()
    if suffix in {".gb", ".gbk"}:
        return "genbank"
    if suffix in {".fa", ".fasta"}:
        return "fasta"

    with input_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                return "fasta"
            if stripped.startswith("LOCUS"):
                return "genbank"
            break

    raise ValueError(
        f"Could not infer input format from {input_path}; provide --format explicitly."
    )


def read_records(input_path: Path, file_format: str) -> List[Record]:
    try:
        from Bio import SeqIO  # type: ignore
    except ImportError as exc:  # pragma: no cover - depends on runtime env
        raise RuntimeError(
            "Biopython is required to read FASTA/GenBank inputs."
        ) from exc

    records: List[Record] = []
    with input_path.open("r", encoding="utf-8", errors="replace") as handle:
        for seq_record in SeqIO.parse(handle, file_format):
            records.append(
                Record(
                    record_id=str(seq_record.id),
                    description=str(seq_record.description),
                    sequence=str(seq_record.seq).upper(),
                )
            )

    if not records:
        raise ValueError(f"No sequence records found in {input_path}.")
    return records


def choose_record(records: Sequence[Record], seq_id: Optional[str], logger: logging.Logger) -> Record:
    if seq_id is None:
        if len(records) > 1:
            logger.warning(
                "Input contains %d records; --seq_id was not provided, so the first record (%s) will be used.",
                len(records),
                records[0].record_id,
            )
        return records[0]

    for record in records:
        if record.record_id == seq_id:
            return record
    available = ", ".join(record.record_id for record in records)
    raise ValueError(
        f"Record ID {seq_id!r} was not found. Available record IDs: {available}"
    )


def resolve_model_paths(model_dir: Path) -> List[Path]:
    expected = [model_dir / f"spliceai{i}.h5" for i in range(1, 6)]
    if all(path.exists() for path in expected):
        return expected

    alternatives = []
    for i in range(1, 6):
        candidates = [
            model_dir / f"spliceai{i}.keras",
            model_dir / f"spliceai{i}.hdf5",
            model_dir / f"spliceai{i}",
        ]
        found = [path for path in candidates if path.exists()]
        if len(found) == 1:
            alternatives.append(found[0])
            continue
        if len(found) > 1:
            raise ValueError(
                f"Multiple model files matched index {i} in {model_dir}: "
                + ", ".join(str(path) for path in found)
            )
        raise FileNotFoundError(
            f"Could not find a model file for index {i} in {model_dir}. "
            "Expected spliceai1.h5 ... spliceai5.h5."
        )
    return alternatives


def import_model_loader():
    try:
        from tensorflow.keras.models import load_model  # type: ignore
        return load_model
    except ImportError:
        try:
            from keras.models import load_model  # type: ignore
            return load_model
        except ImportError as exc:  # pragma: no cover - depends on runtime env
            raise RuntimeError(
                "TensorFlow/Keras is required to load SpliceAI models."
            ) from exc


def one_hot_encode(sequence: str) -> np.ndarray:
    import numpy as np

    mapping = {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3,
        "U": 3,
    }
    encoded = np.zeros((len(sequence), 4), dtype=np.float32)
    for idx, base in enumerate(sequence.upper()):
        col = mapping.get(base)
        if col is not None:
            encoded[idx, col] = 1.0
    return encoded


def load_spliceai_models(model_paths: Sequence[Path], logger: logging.Logger):
    load_model = import_model_loader()
    models = []
    for path in model_paths:
        logger.info("Loading model: %s", path)
        models.append(load_model(str(path), compile=False))
    return models


def average_model_predictions(models, encoded_sequence: np.ndarray) -> np.ndarray:
    import numpy as np

    x = encoded_sequence[None, :, :]
    predictions = [np.asarray(model.predict(x, verbose=0)) for model in models]
    return np.mean(predictions, axis=0)


def extract_position_scores(
    prediction: np.ndarray, sequence_length: int
) -> Tuple[np.ndarray, np.ndarray]:
    import numpy as np

    if prediction.ndim != 3 or prediction.shape[0] != 1 or prediction.shape[2] < 3:
        raise ValueError(f"Unexpected prediction shape: {prediction.shape}")

    core = prediction[0]
    if core.shape[0] < sequence_length:
        raise ValueError(
            f"Prediction length {core.shape[0]} is shorter than sequence length {sequence_length}."
        )

    if core.shape[0] > sequence_length:
        start = (core.shape[0] - sequence_length) // 2
        core = core[start : start + sequence_length]

    acceptor = np.asarray(core[:, 1], dtype=np.float64)
    donor = np.asarray(core[:, 2], dtype=np.float64)
    return acceptor, donor


def write_tsv(
    output_path: Path,
    record: Record,
    acceptor: np.ndarray,
    donor: np.ndarray,
    logger: logging.Logger,
) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(
            "\t".join(
                [
                    "record_id",
                    "position_1based",
                    "base",
                    "acceptor_score",
                    "donor_score",
                    "max_score",
                ]
            )
            + "\n"
        )
        for idx, base in enumerate(record.sequence, start=1):
            a = float(acceptor[idx - 1])
            d = float(donor[idx - 1])
            handle.write(
                f"{record.record_id}\t{idx}\t{base}\t{a:.6f}\t{d:.6f}\t{max(a, d):.6f}\n"
            )
    logger.info("Wrote %s", output_path)


def write_filtered_tsv(
    output_path: Path,
    record: Record,
    acceptor: np.ndarray,
    donor: np.ndarray,
    min_score: float,
    logger: logging.Logger,
) -> int:
    kept = 0
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(
            "\t".join(
                [
                    "record_id",
                    "position_1based",
                    "base",
                    "acceptor_score",
                    "donor_score",
                    "max_score",
                ]
            )
            + "\n"
        )
        for idx, base in enumerate(record.sequence, start=1):
            a = float(acceptor[idx - 1])
            d = float(donor[idx - 1])
            m = max(a, d)
            if m < min_score:
                continue
            kept += 1
            handle.write(
                f"{record.record_id}\t{idx}\t{base}\t{a:.6f}\t{d:.6f}\t{m:.6f}\n"
            )
    logger.info(
        "Filtered rows passing max(score) >= %.6f: %d of %d",
        min_score,
        kept,
        len(record.sequence),
    )
    logger.info("Wrote %s", output_path)
    return kept


def write_plot(
    output_path: Path,
    record: Record,
    acceptor: np.ndarray,
    donor: np.ndarray,
    width: float,
    height: float,
    dpi: int,
    title: Optional[str],
    logger: logging.Logger,
) -> None:
    import numpy as np
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    positions = np.arange(1, len(record.sequence) + 1, dtype=np.int64)

    plt.figure(figsize=(width, height), dpi=dpi)
    plt.plot(positions, acceptor, label="Acceptor", linewidth=1.0, color="#1f77b4")
    plt.plot(positions, donor, label="Donor", linewidth=1.0, color="#d62728")
    plt.xlabel("Position (1-based)")
    plt.ylabel("SpliceAI score")
    plt.ylim(0.0, 1.0)
    plt.xlim(1, max(1, len(record.sequence)))
    plt.legend(frameon=False, loc="upper right")
    plt.title(title if title else record.record_id)
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()
    logger.info("Wrote %s", output_path)


def ensure_writable_outputs(paths: Sequence[Path], force: bool) -> None:
    existing = [path for path in paths if path.exists()]
    if existing and not force:
        raise FileExistsError(
            "Refusing to overwrite existing outputs without --force: "
            + ", ".join(str(path) for path in existing)
        )


def collect_versions() -> str:
    import numpy as np

    parts = [f"Python={sys.version.split()[0]}"]
    module_versions = {
        "numpy": np.__version__,
    }
    try:
        import matplotlib as mpl  # type: ignore

        module_versions["matplotlib"] = mpl.__version__
    except Exception:
        pass
    try:
        import Bio  # type: ignore

        module_versions["biopython"] = Bio.__version__
    except Exception:
        pass
    try:
        import tensorflow as tf  # type: ignore

        module_versions["tensorflow"] = tf.__version__
    except Exception:
        pass
    try:
        import keras  # type: ignore

        module_versions["keras"] = getattr(keras, "__version__", "unknown")
    except Exception:
        pass

    for name, version in sorted(module_versions.items()):
        parts.append(f"{name}={version}")
    return ", ".join(parts)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    script_dir = Path(__file__).resolve().parent
    input_path = Path(args.input)
    output_prefix = Path(args.output_prefix)
    if not output_prefix.is_absolute():
        output_prefix = script_dir / output_prefix
    model_dir = Path(args.model_dir) if args.model_dir else script_dir / "splice_models"
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    log_path = Path(f"{output_prefix}.log")
    scores_path = Path(f"{output_prefix}.spliceai_scores.tsv")
    filtered_path = Path(f"{output_prefix}.spliceai_scores.filtered.tsv")
    plot_path = Path(f"{output_prefix}.spliceai_scores.png")

    ensure_writable_outputs(
        [log_path, scores_path, filtered_path, plot_path], force=args.force
    )

    logger = configure_logging(log_path)
    logger.info("Starting SpliceAI sequence scan")
    logger.info("Arguments: %s", vars(args))
    logger.info("Resolved model directory: %s", model_dir)
    logger.info("Software versions: %s", collect_versions())

    file_format = args.format
    if file_format == "auto":
        file_format = detect_format(input_path)
    logger.info("Using input format: %s", file_format)

    records = read_records(input_path, file_format)
    record = choose_record(records, args.seq_id, logger)
    logger.info(
        "Selected record: id=%s length=%d description=%s",
        record.record_id,
        len(record.sequence),
        record.description,
    )

    if len(record.sequence) == 0:
        raise ValueError("Selected record has zero length.")

    model_paths = resolve_model_paths(model_dir)
    logger.info("Resolved model files: %s", ", ".join(str(path) for path in model_paths))
    models = load_spliceai_models(model_paths, logger)

    context = 10000
    padded_sequence = "N" * (context // 2) + record.sequence + "N" * (context // 2)
    encoded = one_hot_encode(padded_sequence)
    logger.info("Encoded padded sequence with length %d", encoded.shape[0])

    prediction = average_model_predictions(models, encoded)
    acceptor, donor = extract_position_scores(prediction, len(record.sequence))
    logger.info("Prediction shape: %s", prediction.shape)

    write_tsv(scores_path, record, acceptor, donor, logger)
    write_filtered_tsv(
        filtered_path,
        record,
        acceptor,
        donor,
        float(args.min_score),
        logger,
    )
    write_plot(
        plot_path,
        record,
        acceptor,
        donor,
        float(args.plot_width),
        float(args.plot_height),
        int(args.dpi),
        args.title,
        logger,
    )

    logger.info("Completed successfully")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
