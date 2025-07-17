# run_phase1.py

import os
import subprocess

print("=== Phase 1: Degron Data Pipeline ===")

script_dir = os.path.join("phase1", "scripts")

steps = [
    ("Parsing source files...", "parse_sources.py"),
    ("Extracting sequence features...", "extract_features.py"),
    ("Labeling proteins...", "label_dataset.py"),
    ("Merging mutation-related features...", "merge_features.py"),
    ("Generating sequence logos...", "generate_sequence_logos.py")
]

for i, (message, script) in enumerate(steps, 1):
    print(f"{message}")
    result = subprocess.run(["python", os.path.join(script_dir, script)])
    if result.returncode != 0:
        print(f"{message} failed: {script}")
        break
else:
    print("Phase 1 complete.")
