#!/usr/bin/env python3
"""Run Chapter 26's actual WeChat code from scratch and compare its outputs.

Requires the chapter's declared R packages, Rscript, lxml, PyYAML and Pillow.
The working directory must not exist. No project data or R workspace is copied
there: all three inputs are downloaded by the code visible in the article.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
import tempfile
from pathlib import Path
from xml.etree import ElementTree

from PIL import Image

from reader_reproducibility import public_blocks, script_text, validate_reader_contract


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table_values(path: Path) -> dict[tuple[str, str], float]:
    with path.open() as stream:
        reader = csv.reader(stream, delimiter="\t")
        samples = next(reader)[1:]
        return {(row[0], sample): float(value)
                for row in reader for sample, value in zip(samples, row[1:])}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--draft-json", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    project = args.project_root.resolve()
    qmd = project / "chapters/26-community-composition.qmd"
    source = qmd.read_text(encoding="utf-8")
    content = json.loads(args.draft_json.read_text(encoding="utf-8"))["content"]
    contract = validate_reader_contract(qmd, content, project)
    blocks = public_blocks(source, content)
    # No R profiles, saved workspace, or pre-existing input files are used.
    args.work_dir.mkdir(parents=True, exist_ok=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    assert not list(args.work_dir.iterdir())
    with tempfile.TemporaryDirectory(prefix="ch26-public-code-") as directory:
        script = Path(directory) / "public.R"
        script.write_text(script_text(source, blocks), encoding="utf-8")
        result = subprocess.run(["Rscript", "--vanilla", str(script)],
                                cwd=args.work_dir, capture_output=True, text=True, timeout=900)
    args.output.with_suffix(".R.log").write_text(result.stdout + result.stderr, encoding="utf-8")
    checks = [{"check": "public_code_executes", "passed": result.returncode == 0}]
    report = {"chapter": 26, "contract": contract, "fresh_working_directory": True,
              "R_vanilla": True, "preinstalled_declared_packages": True, "checks": checks}
    if result.returncode == 0:
        for filename in ("otutab.tsv", "taxonomy.tsv", "metadata.tsv"):
            downloaded = args.work_dir / "data/small" / filename
            checks.append({"check": "download_" + filename,
                           "passed": sha256(downloaded) == sha256(project / "data/small" / filename)})
        for rank in ("phylum", "genus"):
            relative = Path("results/26-community-composition")
            values = table_values(args.work_dir / relative / f"reader-{rank}-relative.tsv")
            with (project / relative / f"{rank}-relative-abundance.tsv").open() as stream:
                reference = {(r["Taxon"], r["SampleID"]): float(r["RelativeAbundance"])
                             for r in csv.DictReader(stream, delimiter="\t")}
            same_keys = values.keys() == reference.keys()
            difference = max(abs(values[key] - reference[key]) for key in values) if same_keys else math.inf
            checks.append({"check": rank + "_matches_independent_reference",
                           "passed": same_keys and difference < 1e-12,
                           "values_compared": len(values), "maximum_absolute_difference": difference})
        sizes = {"26-phylum-stacked": (183, 112), "26-genus-bubble": (140, 118),
                 "26-genus-heatmap": (183, 122), "26-group-phylum-alluvial": (183, 112)}
        for stem, (width, height) in sizes.items():
            for extension in ("pdf", "svg", "png", "tiff"):
                path = args.work_dir / "figures" / f"{stem}.{extension}"
                checks.append({"check": path.name + "_exported", "passed": path.is_file() and path.stat().st_size > 1000})
            for extension in ("png", "tiff"):
                path = args.work_dir / "figures" / f"{stem}.{extension}"
                with Image.open(path) as im:
                    size_ok = abs(im.width - width / 25.4 * 600) <= 1 and abs(im.height - height / 25.4 * 600) <= 1
                    dpi_ok = all(abs(float(d) - 600) < 1 for d in im.info.get("dpi", (0, 0)))
                checks.append({"check": path.name + "_physical_dimensions_and_600ppi", "passed": size_ok and dpi_ok})
            svg = ElementTree.parse(args.work_dir / "figures" / f"{stem}.svg")
            labels = " ".join("".join(node.itertext()) for node in svg.iter() if node.tag.endswith("}text"))
            checks.append({"check": stem + "_english_labels", "passed": bool(labels) and not re.search(r"[\u3400-\u9fff]", labels)})
            clean_png = args.work_dir / "figures" / f"{stem}.png"
            rendered_png = project / "figures" / clean_png.name
            checks.append({"check": stem + "_matches_rendered_export", "passed": sha256(clean_png) == sha256(rendered_png)})
    report["status"] = "passed" if all(check["passed"] for check in checks) else "failed"
    report["checks_total"] = len(checks)
    report["checks_passed"] = sum(check["passed"] for check in checks)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    raise SystemExit(0 if report["status"] == "passed" else 1)


if __name__ == "__main__":
    main()
