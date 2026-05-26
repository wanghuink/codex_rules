FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

RUN pip install --no-cache-dir \
    biopython \
    matplotlib \
    numpy \
    tensorflow-cpu

COPY SpliceAI/spliceai_sequence_scanner.py /app/spliceai_sequence_scanner.py
COPY SpliceAI/splice_models /app/splice_models

ENTRYPOINT ["python", "/app/spliceai_sequence_scanner.py"]
