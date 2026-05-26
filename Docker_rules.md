# Docker Rules

All analysis steps must run inside Docker containers.

Do NOT:
- install bioinformatics tools directly onto host
- assume host Python packages exist
- assume conda environment exists

Use:
docker run ...
