# Core Rules

- Never hardcode paths.
- All host paths must come from argparse/sys.argv.
- All analysis steps must run in Docker containers.
- Never silently overwrite outputs.
- Ask for clarification if requirements are ambiguous.
- Scientific correctness is more important than speed.
- Never fabricate successful execution.

# Strict follow thest rules
- Data_reules.md
- Docker_rules.md
- Output_convention.md
- Reproducibility.md
- Safety_rules.md
- Scientific.md

# Clarification Policy

Do NOT make assumptions when requirements are ambiguous.

If any of the following are unclear:
- biological intent
- file format
- coordinate system
- strand orientation
- threshold values
- expected outputs
- overwrite behavior
- filtering logic

STOP and ask for clarification.

Never invent:
- reference coordinates
- sample names
- expected biology
- statistical thresholds
- pipeline behavior

Prefer asking questions over making assumptions.

Do not claim code works unless it was actually tested successfully.

When uncertain, explain uncertainty explicitly.

Accuracy is more important than speed.
Scientific correctness takes priority over convenience.
Do not optimize for task completion at the expense of correctness.

