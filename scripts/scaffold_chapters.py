#!/usr/bin/env python3
"""Create missing planned chapter files from tutorial.yaml without overwriting content."""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml


STUB = """---
title: "{display_title}"
draft: true
execute:
  eval: false
---

<!--
Planned chapter. It must not be published until its public data, literature,
commands, outputs, and assertions are added to tutorial.yaml and pass QA.
-->
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    manifest = yaml.safe_load((root / "tutorial.yaml").read_text(encoding="utf-8"))
    created: list[Path] = []

    for chapter in manifest["series"]["chapters"]:
        path = root / chapter["file"]
        if chapter["number"] in {1, 22, 24} or path.exists():
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        display_title = f"第 {chapter['number']:02d} 篇 · {chapter['title']}"
        path.write_text(STUB.format(display_title=display_title), encoding="utf-8")
        created.append(path)

    print(f"created={len(created)}")
    for path in created:
        print(path.relative_to(root))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
