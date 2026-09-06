#!/usr/bin/env python3
"""Keep a chapter's executable R code identical across QMD, script and WeChat.

This is a packaging check, not a substitute for running the extracted public
code in a fresh R session and an empty working directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import yaml
from lxml import html


def frontmatter(source: str) -> dict:
    match = re.match(r"\A---\n(.*?)\n---(?:\n|$)", source, re.S)
    return (yaml.safe_load(match.group(1)) or {}) if match else {}


def normalize_code(code: str) -> str:
    return "\n".join(line.rstrip() for line in code.replace("\r\n", "\n").splitlines()).strip()


def reader_blocks(source: str) -> list[tuple[str, str]]:
    blocks = []
    required = frontmatter(source).get("reader-reproduction", {}).get("required", False)
    for match in re.finditer(r"^```\{r(?:[ ,][^}]*)?\}\n(.*?)^```[ \t]*$", source, re.M | re.S):
        raw = match.group(1)
        option_lines = re.findall(r"^#\| ?(.*)$", raw, re.M)
        options = yaml.safe_load("\n".join(option_lines)) or {}
        if options.get("eval") is False:
            continue
        label = str(options.get("label", ""))
        if not label.startswith(("reader-", "fig-")):
            if required:
                raise ValueError("Every executed R chunk needs a reader- or fig- label")
            continue
        if options.get("echo") is False or options.get("include") is False:
            raise ValueError(f"Required reader code is hidden: {label}")
        code = normalize_code(re.sub(r"^#\|.*\n?", "", raw, flags=re.M))
        if not code:
            raise ValueError(f"Empty reader code: {label}")
        blocks.append((label, code))
    labels = [label for label, _ in blocks]
    if len(labels) != len(set(labels)):
        raise ValueError("Duplicate reader code labels")
    if required and not blocks:
        raise ValueError("Required reader code is missing")
    return blocks


def script_text(source: str, blocks: list[tuple[str, str]]) -> str:
    title = frontmatter(source).get("title", "R analysis")
    header = f"# {title}\n# Run sequentially in a new working directory.\n# Required packages: ggplot2, ggalluvial, ragg, svglite.\n"
    return header + "\n" + "\n\n".join(code for _, code in blocks) + "\n"


def public_blocks(source: str, content: str) -> list[tuple[str, str]]:
    expected = reader_blocks(source)
    displayed = [normalize_code(pre.text_content()) for pre in html.fromstring(content).xpath(".//pre")]
    selected = []
    cursor = 0
    for label, code in expected:
        try:
            position = displayed.index(code, cursor)
        except ValueError as exc:
            raise ValueError(f"Required code missing, changed, or out of order in WeChat: {label}") from exc
        selected.append((label, displayed[position]))
        cursor = position + 1
    return selected


def validate_reader_contract(source_qmd: Path, content: str, project: Path | None = None) -> dict:
    source = source_qmd.read_text(encoding="utf-8")
    config = frontmatter(source).get("reader-reproduction", {})
    if not config.get("required", False):
        return {"required": False}
    blocks = public_blocks(source, content)
    expected_script = script_text(source, blocks)
    report = {
        "required": True,
        "status": "passed",
        "code_blocks": len(blocks),
        "code_sha256": hashlib.sha256(expected_script.encode()).hexdigest(),
    }
    if project is not None:
        script_path = (project / config["script"]).resolve()
        if not script_path.is_relative_to(project.resolve()):
            raise ValueError("Reader script must be inside the project")
        if script_path.read_text(encoding="utf-8") != expected_script:
            raise ValueError("Downloadable R script differs from the public code")
        report["script"] = config["script"]
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qmd", type=Path, required=True)
    parser.add_argument("--article", type=Path)
    parser.add_argument("--project-root", type=Path)
    parser.add_argument("--emit-script", action="store_true")
    args = parser.parse_args()
    source = args.qmd.read_text(encoding="utf-8")
    if args.article:
        content = args.article.read_text(encoding="utf-8")
        report = validate_reader_contract(args.qmd, content, args.project_root)
        blocks = public_blocks(source, content)
    else:
        blocks = reader_blocks(source)
        report = {"code_blocks": len(blocks)}
    print(script_text(source, blocks), end="") if args.emit_script else print(json.dumps(report))


if __name__ == "__main__":
    main()
