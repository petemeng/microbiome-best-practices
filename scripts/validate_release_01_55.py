#!/usr/bin/env python3
"""Combine established QA with release gates for Articles 01–55."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument(
        "--output", type=Path,
        default=Path("results/articles-01-55-release-validation.json"),
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    manifest_path = root / "tutorial.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    chapters = manifest["series"]["chapters"][:55]
    components = {
        "articles_01_29_full_qa": root / "qa_report.json",
        "series_contract": root / "results" / "scaffold_validation.json",
        "articles_30_35": root / "results" / "articles-30-35-validation.json",
        "articles_36_40": root / "results" / "articles-36-40-validation.json",
        "articles_41_45": root / "results" / "articles-41-45-validation.json",
        "articles_46_50": root / "results" / "articles-46-50-validation.json",
        "articles_51_55": root / "results" / "articles-51-55-validation.json",
    }
    checks: list[dict[str, object]] = []

    def check(identifier: str, passed: bool, detail: object) -> None:
        checks.append({
            "id": identifier,
            "status": "PASS" if passed else "FAIL",
            "detail": detail,
        })

    component_records: dict[str, dict[str, object]] = {}
    for name, path in components.items():
        exists = path.is_file()
        check(f"component-{name}-exists", exists, str(path))
        if not exists:
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        component_records[name] = {
            "path": str(path),
            "sha256": sha256(path),
            "status": payload.get("status"),
            "checks_passed": payload.get("checks_passed"),
            "checks_failed": payload.get("checks_failed"),
        }
        check(
            f"component-{name}-passed",
            payload.get("status") == "passed",
            payload.get("status"),
        )

    scaffold_path = components["series_contract"]
    scaffold = (
        json.loads(scaffold_path.read_text(encoding="utf-8"))
        if scaffold_path.is_file() else {}
    )
    check(
        "formal-series-coverage",
        scaffold.get("pilot_articles") == list(range(1, 56)),
        scaffold.get("pilot_articles"),
    )
    check("manifest-formal-count", len(chapters) == 55, len(chapters))
    check(
        "manifest-total-count",
        manifest["series"]["total_articles"] == 55,
        manifest["series"]["total_articles"],
    )

    html_count = 0
    for chapter in chapters:
        source = Path(chapter["file"])
        html = root / "_site" / (
            "index.html" if source.name == "index.qmd" else source.with_suffix(".html")
        )
        exists = html.is_file() and html.stat().st_size > 30_000
        check(f"article-{int(chapter['number']):02d}-html", exists, str(html))
        html_count += int(exists)
    check("formal-html-count", html_count == 55, html_count)

    failed = [row for row in checks if row["status"] == "FAIL"]
    run_material = "|".join([
        sha256(manifest_path),
        *(str(record["sha256"]) for record in component_records.values()),
    ])
    payload = {
        "status": "passed" if not failed else "failed",
        "scope": "Articles 01–55 local release assets",
        "manifest": str(manifest_path),
        "manifest_hash": sha256(manifest_path),
        "run_key": hashlib.sha256(run_material.encode("utf-8")).hexdigest()[:16],
        "formal_count": 55,
        "html_count": html_count,
        "components": component_records,
        "checks_total": len(checks),
        "checks_passed": len(checks) - len(failed),
        "checks_failed": len(failed),
        "failed_checks": failed,
        "checks": checks,
    }
    output = args.output if args.output.is_absolute() else root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(
        {key: value for key, value in payload.items() if key != "checks"},
        ensure_ascii=False,
        indent=2,
    ))
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
