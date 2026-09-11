#!/usr/bin/env python3
"""Keep the full R workflow faithful to QMD and WeChat excerpts faithful to it.

This is a packaging check, not a substitute for running the extracted public
code in a fresh R session and an empty working directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from functools import lru_cache
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
    config = frontmatter(source).get("reader-reproduction", {})
    required = config.get("required", False)
    all_executed = config.get("label-policy") == "all-executed"
    for match in re.finditer(r"^```\{r(?:[ ,][^}]*)?\}\n(.*?)^```[ \t]*$", source, re.M | re.S):
        raw = match.group(1)
        option_lines = re.findall(r"^#\| ?(.*)$", raw, re.M)
        options = yaml.safe_load("\n".join(option_lines)) or {}
        if options.get("eval") is False:
            continue
        label = str(options.get("label", ""))
        if not label or (not all_executed and not label.startswith(("reader-", "fig-"))):
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
    metadata = frontmatter(source)
    title = metadata.get("title", "R analysis")
    packages = metadata.get("reader-reproduction", {}).get("packages", ["ggplot2", "ggalluvial", "ragg", "svglite"])
    header = f"# {title}\n# Run sequentially in a new working directory.\n# Required packages: {', '.join(packages)}.\n"
    return header + "\n" + "\n\n".join(code for _, code in blocks) + "\n"


GENERIC_HELPERS = ('font_pub', 'pal_pub', 'scale_color_pub', 'scale_fill_pub', 'theme_pub', 'save_pub')


@lru_cache(maxsize=256)
def decisive_code(code: str) -> str:
    """Omit only named, top-level plotting bootstrap assignments from WeChat.

    R parses expression boundaries; nested functions, strings and scientific
    input/analysis functions cannot accidentally be removed by a text regex.
    The complete executable code remains unchanged in QMD and its R download.
    """
    if not re.search(r'^\s*(?:' + '|'.join(GENERIC_HELPERS) + r')\s*(?:<-|=)', code, re.M):
        return code
    parser = r'''
code <- paste(readLines(file("stdin"), warn = FALSE), collapse = "\n")
parsed <- parse(text = code, keep.source = TRUE)
refs <- attr(parsed, "srcref")
names <- c("font_pub", "pal_pub", "scale_color_pub", "scale_fill_pub", "theme_pub", "save_pub")
for (i in seq_along(parsed)) {
  expr <- parsed[[i]]
  if (is.call(expr) && is.symbol(expr[[1]]) && as.character(expr[[1]]) %in% c("<-", "=") &&
      is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% names) {
    ref <- refs[[i]]
    cat(ref[1], ref[2], ref[3], ref[4], sep = "\t"); cat("\n")
  }
}
'''
    run = subprocess.run(['Rscript', '--vanilla', '-e', parser], input=code,
                         text=True, capture_output=True, check=True)
    lines = code.splitlines()
    omitted = set()
    for row in run.stdout.splitlines():
        first, start, last, end = map(int, row.split('\t'))
        prefix = lines[first - 1].encode()[:start - 1].decode()
        suffix = lines[last - 1].encode()[end:].decode().strip()
        if prefix.strip() or (suffix and not suffix.startswith('#')):
            raise ValueError('A generic helper shares a line with another expression')
        omitted.update(range(first - 1, last))
    return normalize_code(re.sub(r'\n{3,}', '\n\n', '\n'.join(
        line for index, line in enumerate(lines) if index not in omitted)))


def omit_generic_reader_helpers(source: str, main) -> int:
    """Transform verified rendered source blocks, never the runnable R download."""
    if not frontmatter(source).get('reader-reproduction', {}).get('required'):
        return 0
    replacements = {code: decisive_code(code) for _, code in reader_blocks(source)}
    changed = 0
    for pre in list(main.xpath('.//pre')):
        original = normalize_code(pre.text_content())
        if original not in replacements or replacements[original] == original:
            continue
        replacement = replacements[original]
        if replacement:
            for child in list(pre):pre.remove(child)
            pre.text = replacement
        else:
            pre.getparent().remove(pre)
        changed += 1
    return changed


def public_blocks(source: str, content: str) -> list[tuple[str, str]]:
    expected = [(label, decisive_code(code)) for label, code in reader_blocks(source)]
    displayed = [normalize_code(pre.text_content()) for pre in html.fromstring(content).xpath(".//pre")]
    selected = []
    cursor = 0
    for label, code in expected:
        if not code:
            continue
        try:
            position = displayed.index(code, cursor)
        except ValueError as exc:
            raise ValueError(f"Required code missing, changed, or out of order in WeChat: {label}") from exc
        selected.append((label, displayed[position]))
        cursor = position + 1
    return selected


def coalesce_rendered_code(source: str, main) -> int:
    """Join code split by knitr plot events, only after verifying every segment.

    Never insert code into a stale cell: each existing source segment must match
    the QMD exactly, allowing whitespace only at boundaries between segments.
    Computed outputs remain in place and are not used to reconstruct source.
    """
    if not frontmatter(source).get('reader-reproduction', {}).get('required'):
        return 0
    cells = main.xpath('.//*[contains(concat(" ", normalize-space(@class), " "), " cell ")]')
    cursor = 0
    merged = 0
    for label, expected in reader_blocks(source):
        for index in range(cursor, len(cells)):
            pres = cells[index].xpath('.//pre[contains(concat(" ", normalize-space(@class), " "), " sourceCode ")]')
            segments = [normalize_code(pre.text_content()) for pre in pres]
            position = 0
            matched = bool(segments)
            for segment in segments:
                start = expected.find(segment, position)
                if start < 0 or expected[position:start].strip():
                    matched = False
                    break
                position = start + len(segment)
            if not matched or expected[position:].strip():
                continue
            if len(pres) > 1:
                first = pres[0]
                for child in list(first):
                    first.remove(child)
                first.text = expected
                for pre in pres[1:]:
                    pre.getparent().remove(pre)
                merged += 1
            cursor = index + 1
            break
    return merged


def validate_reader_contract(source_qmd: Path, content: str, project: Path | None = None) -> dict:
    source = source_qmd.read_text(encoding="utf-8")
    config = frontmatter(source).get("reader-reproduction", {})
    if not config.get("required", False):
        return {"required": False}
    blocks = public_blocks(source, content)
    full_blocks = reader_blocks(source)
    expected_script = script_text(source, full_blocks)
    report = {
        "required": True,
        "status": "passed",
        "code_blocks": len(blocks),
        "public_code_policy": "analysis code with complete downloadable R workflow; generic plotting bootstrap omitted",
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
        blocks = reader_blocks(source)
    else:
        blocks = reader_blocks(source)
        report = {"code_blocks": len(blocks)}
    print(script_text(source, blocks), end="") if args.emit_script else print(json.dumps(report))


if __name__ == "__main__":
    main()
