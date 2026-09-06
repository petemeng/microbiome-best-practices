#!/usr/bin/env python3
"""Build WeChat-ready review articles from the validated Quarto HTML site.

The script is intentionally offline: it sanitizes and inlines styles, optimizes
article images, creates deterministic covers from each article's representative
figure, and writes resumable JSON payloads.
Uploading images and creating official-account drafts remain separate actions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import urljoin

import yaml
from lxml import etree, html
from PIL import Image, ImageOps

from reader_reproducibility import validate_reader_contract


MAX_TITLE_CHARS = 64
MAX_DIGEST_CHARS = 120
MAX_THUMB_BYTES = 64 * 1024
MAX_ARTICLE_IMAGE_BYTES = 950 * 1024
MAX_ARTICLE_IMAGE_WIDTH = 1600

ROOT_STYLE = (
    "font-family:-apple-system,BlinkMacSystemFont,'PingFang SC',"
    "'Hiragino Sans GB','Noto Sans CJK SC',sans-serif;"
    "color:#2d2d2d;font-size:16px;line-height:1.85;"
    "word-break:break-word;"
)
STYLES = {
    "h2": (
        "margin:38px 0 18px;padding:10px 16px;border-left:4px solid #7c9970;"
        "background:#eef3ea;color:#203124;font-size:22px;line-height:1.45;"
    ),
    "h3": "margin:30px 0 14px;color:#314735;font-size:19px;line-height:1.5;",
    "h4": "margin:24px 0 12px;color:#4b5f4e;font-size:17px;line-height:1.5;",
    "p": "margin:0 0 18px;color:#2d2d2d;font-size:16px;line-height:1.85;",
    "ul": "margin:0 0 18px;padding-left:1.35em;color:#2d2d2d;line-height:1.85;",
    "ol": "margin:0 0 18px;padding-left:1.35em;color:#2d2d2d;line-height:1.85;",
    "li": "margin:0 0 10px;",
    "pre": (
        "margin:20px 0;padding:16px;background:#1f2a20;border-radius:12px;"
        "white-space:pre-wrap;word-break:break-all;color:#eef5ea;"
        "font-family:Menlo,Consolas,monospace;font-size:12px;line-height:1.7;"
    ),
    "code": (
        "padding:2px 5px;background:#f1ede4;border-radius:5px;"
        "font-family:Menlo,Consolas,monospace;font-size:0.9em;color:#8a5a1f;"
    ),
    "a": "color:#8a6428;text-decoration:none;border-bottom:1px solid #cfb17b;",
    "blockquote": (
        "margin:22px 0;padding:16px 18px 2px;background:#f6f1e7;"
        "border-left:4px solid #d7a35b;border-radius:10px;"
    ),
    "table": (
        "width:100%;margin:18px 0;border-collapse:collapse;table-layout:auto;"
        "font-size:13px;line-height:1.55;"
    ),
    "th": "padding:8px 6px;border:1px solid #d8d8d8;background:#eef3ea;text-align:left;",
    "td": "padding:8px 6px;border:1px solid #dedede;vertical-align:top;",
    "figure": "margin:22px 0;text-align:center;",
    "figcaption": "margin:8px 8px 20px;color:#666;font-size:13px;line-height:1.65;text-align:left;",
    "img": "display:block;width:100%;max-width:100%;height:auto;margin:18px auto;border-radius:10px;",
    "hr": "margin:28px auto;width:72px;height:1px;border:0;background:#d8ccb7;",
}
CALLOUT_STYLES = {
    "caution": "background:#fff4ed;border-left:4px solid #d9734e;",
    "warning": "background:#fff8e8;border-left:4px solid #d7a35b;",
    "important": "background:#f4eef8;border-left:4px solid #8c6ca8;",
    "tip": "background:#edf6f1;border-left:4px solid #5f9a7d;",
    "note": "background:#eef4f7;border-left:4px solid #5f879a;",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--manifest", default="tutorial.yaml")
    parser.add_argument("--qa-report", default="qa_report.json")
    parser.add_argument("--site-dir", default="_site")
    parser.add_argument("--output-dir", default="rendered/wechat_review_01_55")
    parser.add_argument(
        "--fallback-bundle",
        help=(
            "Optional previously verified review bundle used only when a "
            "rendered site article is missing. The report records every "
            "fallback chapter."
        ),
    )
    parser.add_argument("--formal-count", type=int, default=55)
    parser.add_argument(
        "--chapters", type=int, nargs="+",
        help="Build only these chapter numbers, in manifest order (for scoped revisions).",
    )
    parser.add_argument("--author", default="Peter")
    parser.add_argument(
        "--review-url",
        default="https://github.com/petemeng/microbiome-best-practices",
    )
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def truncate(text: str, limit: int) -> str:
    normalized = re.sub(r"\s+", " ", text).strip()
    if len(normalized) <= limit:
        return normalized
    return normalized[: max(1, limit - 1)].rstrip() + "…"


def select_chapters(
    chapters: list[dict[str, Any]],
    formal_count: int,
    selected: list[int] | None,
) -> list[dict[str, Any]]:
    if not 1 <= formal_count <= len(chapters):
        raise ValueError("formal-count must be within the manifest chapter count")
    formal = chapters[:formal_count]
    if selected is None:
        return formal
    if not selected or len(set(selected)) != len(selected):
        raise ValueError("Selected chapter numbers must be nonempty and unique")
    available = {int(chapter["number"]) for chapter in formal}
    if set(selected) - available:
        raise ValueError("Selected chapters must be within the formal manifest scope")
    return [chapter for chapter in formal if int(chapter["number"]) in selected]


def source_description(source_qmd: Path) -> str | None:
    source = source_qmd.read_text(encoding="utf-8")
    match = re.match(r"\A---\s*\n(.*?)\n---(?:\s*\n|$)", source, re.S)
    metadata = yaml.safe_load(match.group(1)) if match else None
    description = metadata.get("description") if isinstance(metadata, dict) else None
    return description.strip() if isinstance(description, str) and description.strip() else None


def article_readability(content: str) -> dict[str, int | float | None]:
    """Report reading load without treating a word-count threshold as quality."""
    document = html.fromstring(content)
    visible_chars = len(normalized_text(document))
    code_chars = sum(len(normalized_text(node)) for node in document.xpath(".//pre"))
    first_image = re.search(r"<img\b", content, re.I)
    before_image = (
        len(normalized_text(html.fromstring(content[:first_image.start()])))
        if first_image else None
    )
    return {
        "visible_chars": visible_chars,
        "preformatted_block_count": len(document.xpath(".//pre")),
        "preformatted_chars": code_chars,
        "preformatted_share": round(code_chars / max(1, visible_chars), 4),
        "chars_before_first_image": before_image,
    }


def class_tokens(element: etree._Element) -> set[str]:
    return set((element.get("class") or "").split())


def site_html_path(site_dir: Path, qmd_path: str) -> Path:
    qmd = Path(qmd_path)
    if qmd.name == "index.qmd":
        return site_dir / "index.html"
    return site_dir / qmd.with_suffix(".html")


def missing_local_images(source_html: Path) -> list[str]:
    if not source_html.exists():
        return []
    document = html.parse(str(source_html)).getroot()
    missing: list[str] = []
    for image_element in document.xpath("//img"):
        src = image_element.get("src") or ""
        if not src or src.startswith(("http://", "https://", "data:")):
            continue
        if not (source_html.parent / src).resolve().exists():
            missing.append(src)
    return missing


def create_cover(representative_image: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    image = Image.new("RGB", (900, 383), "white")
    with Image.open(representative_image) as opened:
        figure = ImageOps.exif_transpose(opened)
        if figure.mode in {"RGBA", "LA"} or "transparency" in figure.info:
            rgba = figure.convert("RGBA")
            background = Image.new("RGBA", rgba.size, "white")
            background.alpha_composite(rgba)
            figure = background.convert("RGB")
        else:
            figure = figure.convert("RGB")
        figure = ImageOps.contain(figure, image.size, Image.Resampling.LANCZOS)
        image.paste(
            figure,
            ((image.width - figure.width) // 2, (image.height - figure.height) // 2),
        )
    for quality in (82, 74, 68, 62, 56, 50, 44, 38, 34):
        image.save(output, "JPEG", quality=quality, optimize=True, progressive=True, subsampling=2)
        if output.stat().st_size <= MAX_THUMB_BYTES:
            return
    if output.stat().st_size > MAX_THUMB_BYTES:
        raise RuntimeError(f"Cover remains above {MAX_THUMB_BYTES} bytes: {output}")


def optimize_article_image(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(source) as opened:
        image = ImageOps.exif_transpose(opened)
        if image.mode in {"RGBA", "LA"} or "transparency" in image.info:
            rgba = image.convert("RGBA")
            background = Image.new("RGBA", rgba.size, "white")
            background.alpha_composite(rgba)
            image = background.convert("RGB")
        else:
            image = image.convert("RGB")
        if image.width > MAX_ARTICLE_IMAGE_WIDTH:
            height = round(image.height * MAX_ARTICLE_IMAGE_WIDTH / image.width)
            image = image.resize((MAX_ARTICLE_IMAGE_WIDTH, height), Image.Resampling.LANCZOS)
        quality = 90
        while True:
            image.save(destination, "JPEG", quality=quality, optimize=True, progressive=True, subsampling=2)
            if destination.stat().st_size <= MAX_ARTICLE_IMAGE_BYTES:
                return
            if quality > 64:
                quality -= 8
                continue
            new_width = max(900, round(image.width * 0.82))
            if new_width >= image.width:
                raise RuntimeError(f"Could not optimize article image: {source}")
            new_height = round(image.height * new_width / image.width)
            image = image.resize((new_width, new_height), Image.Resampling.LANCZOS)
            quality = 82


def remove_unwanted(main: etree._Element) -> None:
    # These labels are hidden on the website, but inline-only WeChat HTML
    # would expose them next to the actual callout title. Preserve the tail.
    for element in main.xpath(
        './/*[contains(concat(" ", normalize-space(@class), " "), " screen-reader-only ")]'
    ):
        element.drop_tree()
    selectors = [
        ".//script",
        ".//style",
        ".//button",
        ".//nav",
        './/*[@id="title-block-header"]',
        './/*[contains(concat(" ", normalize-space(@class), " "), " code-copy-button ")]',
        './/*[contains(concat(" ", normalize-space(@class), " "), " code-annotation-gutter ")]',
        './/*[contains(concat(" ", normalize-space(@class), " "), " anchorjs-link ")]',
        './/*[contains(concat(" ", normalize-space(@class), " "), " header-section-number ")]',
        './/*[contains(concat(" ", normalize-space(@class), " "), " quarto-title-meta ")]',
    ]
    seen: set[etree._Element] = set()
    for selector in selectors:
        for element in main.xpath(selector):
            if element in seen:
                continue
            seen.add(element)
            parent = element.getparent()
            if parent is not None:
                parent.remove(element)


def normalized_text(element: etree._Element) -> str:
    return re.sub(r"\s+", " ", "".join(element.itertext())).strip()


def localize_figure_labels(main: etree._Element) -> None:
    """Use article-local figure numbers; keep chapter prefixes on the website."""
    pattern = re.compile(r"^(\s*)(图|Figure)\s+(\d+(?:\.\d+)+)(?=\s*[:：]|\s*$)")
    labels: dict[str, str] = {}
    for caption in main.xpath(".//figcaption"):
        match = pattern.match(caption.text or "")
        if match:
            local = str(len(labels) + 1)
            labels[match.group(3)] = local
            caption.text = pattern.sub(
                lambda m: f"{m.group(1)}{m.group(2)} {local}",
                caption.text, count=1,
            )
    for link in main.xpath('.//a[starts-with(@href, "#fig-")]'):
        match = pattern.match(link.text or "")
        if match and match.group(3) in labels:
            local = labels[match.group(3)]
            link.text = pattern.sub(
                lambda m: f"{m.group(1)}{m.group(2)} {local}",
                link.text, count=1,
            )


SOURCE_H2 = re.compile(
    r"(?m)^##\s+(.+?)\s+\{#(sec-[A-Za-z0-9_-]+)\}\s*$"
)


def source_section_headings(source_qmd: Path) -> dict[str, str]:
    """Read stable section anchors and their reader-facing titles from QMD."""
    source = source_qmd.read_text(encoding="utf-8")
    return {
        anchor: re.sub(r"[`*_]", "", heading).strip()
        for heading, anchor in SOURCE_H2.findall(source)
    }


def replace_element_text(element: etree._Element, value: str) -> None:
    for child in list(element):
        element.remove(child)
    element.text = value


def sync_topic_heading(main: etree._Element, source_qmd: Path) -> int:
    """Keep a verified fallback body aligned with the current topic heading."""
    expected = source_section_headings(source_qmd).get("sec-theory")
    if not expected:
        return 0
    for heading in main.xpath(".//h2"):
        anchor = heading.get("data-anchor-id") or heading.get("id")
        current = normalized_text(heading)
        if (
            anchor == "sec-theory"
            or current.startswith("理论：")
            or current == "为什么这么做"
            or current == expected
        ):
            if current != expected:
                replace_element_text(heading, expected)
                return 1
            return 0
    raise RuntimeError(f"Could not locate sec-theory heading in {source_qmd}")


def remove_wechat_bootstrap(main: etree._Element) -> int:
    """Remove website-only environment bootstrap blocks from WeChat prose."""
    removed = 0
    for section in list(main.xpath(".//section")):
        headings = section.xpath("./h2[1]")
        if not headings or normalized_text(headings[0]) != "准备工作":
            continue
        parent = section.getparent()
        if parent is not None:
            parent.remove(section)
            removed += 1

    # Keep the filter robust if a future chapter changes the outer heading but
    # retains the generic dependency/data/theme bootstrap disclosure.
    for details in list(main.xpath(".//details")):
        summaries = details.xpath("./summary[1]")
        if not summaries:
            continue
        summary = normalized_text(summaries[0])
        if not (
            summary.startswith("展开：")
            and (
                "安装依赖" in summary
                or "定义作图函数" in summary
                or "出版级函数" in summary
            )
        ):
            continue
        parent = details.getparent()
        if parent is not None:
            parent.remove(details)
            removed += 1
    return removed


def remove_explicit_wechat_omissions(main: etree._Element) -> int:
    removed = 0
    selector = (
        './/*[contains(concat(" ", normalize-space(@class), " "), '
        '" wechat-omit ")]'
    )
    elements = list(main.xpath(selector))
    marked = set(elements)
    for element in elements:
        if any(ancestor in marked for ancestor in element.iterancestors()):
            continue
        parent = element.getparent()
        if parent is not None:
            parent.remove(element)
            removed += 1
    return removed


INSTALL_CALL = re.compile(
    r"^\s*(?:install\.packages|BiocManager::install|"
    r"remotes::install_github|pak::pkg_install)\s*\("
)


def strip_leading_install_calls(code_text: str) -> tuple[str, int]:
    lines = code_text.splitlines()
    prefix: list[str] = []
    index = 0
    while index < len(lines) and (
        not lines[index].strip() or lines[index].lstrip().startswith("#|")
    ):
        prefix.append(lines[index])
        index += 1

    removed = 0
    while index < len(lines) and INSTALL_CALL.match(lines[index]):
        depth = lines[index].count("(") - lines[index].count(")")
        index += 1
        while index < len(lines) and depth > 0:
            depth += lines[index].count("(") - lines[index].count(")")
            index += 1
        removed += 1
        while index < len(lines) and not lines[index].strip():
            index += 1

    if not removed:
        return code_text, 0
    return "\n".join(prefix + lines[index:]).rstrip(), removed


def flatten_code(main: etree._Element) -> int:
    stripped_install_calls = 0
    for pre in main.xpath(".//pre"):
        code_text = "".join(pre.itertext()).rstrip()
        code_text, removed = strip_leading_install_calls(code_text)
        stripped_install_calls += removed
        for child in list(pre):
            pre.remove(child)
        pre.text = code_text
        pre.set("style", STYLES["pre"])
    return stripped_install_calls


def transform_special_blocks(main: etree._Element) -> None:
    for details in main.xpath(".//details"):
        details.tag = "section"
        details.set(
            "style",
            "margin:22px 0;padding:16px 18px;background:#faf8f2;"
            "border:1px solid #e6dfd1;border-radius:10px;",
        )
    for summary in main.xpath(".//summary"):
        summary.tag = "p"
        summary.set(
            "style",
            "margin:0 0 14px;color:#314735;font-size:17px;line-height:1.6;font-weight:bold;",
        )
    for element in main.xpath(
        './/*[contains(concat(" ", normalize-space(@class), " "), " callout ")]'
    ):
        tokens = class_tokens(element)
        kind = next((name for name in CALLOUT_STYLES if f"callout-{name}" in tokens), "note")
        element.tag = "section"
        element.set(
            "style",
            "margin:22px 0;padding:16px 18px 4px;border-radius:10px;"
            + CALLOUT_STYLES[kind],
        )


def apply_inline_styles(main: etree._Element) -> None:
    for tag, style in STYLES.items():
        for element in main.xpath(f".//{tag}"):
            if tag == "code" and element.getparent() is not None and element.getparent().tag == "pre":
                continue
            element.set("style", style)
    for element in main.xpath(".//strong"):
        element.set("style", "color:#203124;font-weight:700;")
    for element in main.xpath(".//em"):
        element.set("style", "color:#555;font-style:italic;")
    for element in main.xpath(".//div"):
        if not element.get("style"):
            element.tag = "section"


def resolve_and_optimize_images(
    main: etree._Element,
    source_html: Path,
    article_dir: Path,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    cache: dict[Path, tuple[str, Path]] = {}
    for index, image_element in enumerate(main.xpath(".//img"), start=1):
        src = image_element.get("src") or ""
        if src.startswith(("http://", "https://", "data:")):
            continue
        source = (source_html.parent / src).resolve()
        if not source.exists():
            raise FileNotFoundError(f"Missing rendered image {src} referenced by {source_html}")
        if source not in cache:
            relative = f"images/image-{len(cache) + 1:02d}.jpg"
            destination = article_dir / relative
            optimize_article_image(source, destination)
            cache[source] = (relative, destination)
        relative, destination = cache[source]
        image_element.set("src", relative)
        image_element.set("style", STYLES["img"])
        with Image.open(destination) as optimized:
            width, height = optimized.size
        record = {
            "source_path": str(source),
            "local_path": str(destination.resolve()),
            "relative_src": relative,
            "sha256": sha256(destination),
            "size_bytes": destination.stat().st_size,
            "width": width,
            "height": height,
        }
        if not any(item["local_path"] == record["local_path"] for item in records):
            records.append(record)
    return records


def strip_unsupported_attributes(main: etree._Element) -> None:
    allowed = {"style", "href", "src", "alt", "title", "colspan", "rowspan"}
    for element in main.iter():
        for key in list(element.attrib):
            if key not in allowed:
                del element.attrib[key]
        if element.tag == "a":
            href = element.get("href") or ""
            # The official-account draft API rejects document-fragment links
            # (for example Quarto's ``href="#fig-..."``) with error 45166.
            # External citation links remain valid; only local anchors and
            # executable pseudo-URLs are stripped.
            if href.startswith(("#", "javascript:")):
                element.attrib.pop("href", None)


def sanitize_article(
    source_html: Path,
    article_dir: Path,
    source_qmd: Path,
) -> tuple[str, list[dict[str, Any]], int, int, int, int]:
    document = html.parse(str(source_html)).getroot()
    mains = document.xpath(
        '//main[contains(concat(" ",normalize-space(@class)," ")," content ")]'
        ' | //main[@id="quarto-document-content"] | //main'
    )
    if not mains:
        # A prior local review bundle wraps its sanitized article in the first
        # section under body. Accept that surface only when explicitly chosen
        # as a fallback source.
        mains = document.xpath("//body/section[1]")
    if not mains:
        raise RuntimeError(f"No article content found in {source_html}")
    main = deepcopy(mains[0])
    main.tag = "section"
    remove_unwanted(main)
    synced_topic_heading_count = sync_topic_heading(main, source_qmd)
    removed_wechat_omit_blocks = remove_explicit_wechat_omissions(main)
    removed_bootstrap_blocks = remove_wechat_bootstrap(main)
    localize_figure_labels(main)
    stripped_install_calls = flatten_code(main)
    transform_special_blocks(main)
    apply_inline_styles(main)
    images = resolve_and_optimize_images(main, source_html, article_dir)
    strip_unsupported_attributes(main)
    main.set("style", ROOT_STYLE)
    content = etree.tostring(main, encoding="unicode", method="html")
    validate_reader_contract(source_qmd, content)
    return (
        content,
        images,
        removed_bootstrap_blocks,
        removed_wechat_omit_blocks,
        stripped_install_calls,
        synced_topic_heading_count,
    )


def local_preview_html(title: str, content: str) -> str:
    return (
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        f"<title>{title}</title></head><body style='margin:0 auto;padding:24px;max-width:760px;'>"
        f"{content}</body></html>"
    )


def build(args: argparse.Namespace) -> dict[str, Any]:
    project = Path(args.project_root).resolve()
    manifest_path = (project / args.manifest).resolve()
    qa_path = (project / args.qa_report).resolve()
    site_dir = (project / args.site_dir).resolve()
    output_dir = (project / args.output_dir).resolve()
    fallback_bundle = (
        (project / args.fallback_bundle).resolve()
        if args.fallback_bundle
        else None
    )
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    qa_report = json.loads(qa_path.read_text(encoding="utf-8"))
    if qa_report.get("status") != "passed":
        raise RuntimeError("qa_report.json.status must be passed before draft generation")
    chapters = manifest.get("series", {}).get("chapters", [])
    if len(chapters) != 55:
        raise RuntimeError(f"Expected 55 manifest chapters, found {len(chapters)}")
    formal = select_chapters(chapters, args.formal_count, getattr(args, "chapters", None))
    output_dir.mkdir(parents=True, exist_ok=True)
    items: list[dict[str, Any]] = []
    for chapter in formal:
        number = int(chapter["number"])
        raw_title = str(chapter["title"])
        qmd_path = str(chapter["file"])
        source_html = site_html_path(site_dir, qmd_path)
        source_surface = "rendered_site"
        source_surface_reason = "current_render"
        site_problem: str | None = None
        if not source_html.exists():
            site_problem = "missing_site_html"
        else:
            missing_images = missing_local_images(source_html)
            if missing_images:
                site_problem = f"missing_site_images:{len(missing_images)}"
        if site_problem is not None:
            if fallback_bundle is None:
                raise FileNotFoundError(
                    f"{site_problem} for rendered article {source_html}"
                )
            fallback_html = fallback_bundle / f"{number:02d}" / "article.html"
            if not fallback_html.exists():
                raise FileNotFoundError(
                    f"Missing site HTML {source_html} and fallback {fallback_html}"
                )
            fallback_missing_images = missing_local_images(fallback_html)
            if fallback_missing_images:
                raise FileNotFoundError(
                    f"Fallback {fallback_html} has missing images: "
                    f"{fallback_missing_images[:3]}"
                )
            source_html = fallback_html
            source_surface = "verified_fallback_bundle"
            source_surface_reason = site_problem
        article_dir = output_dir / f"{number:02d}"
        article_dir.mkdir(parents=True, exist_ok=True)
        title = truncate(f"16S最佳实践｜{number}. {raw_title}", MAX_TITLE_CHARS)
        digest_lead = raw_title if raw_title.endswith(("。", "！", "？", "!", "?")) else f"{raw_title}。"
        digest_text = (
            f"{digest_lead}真实数据、关键分析步骤、结果解释与发表级重绘图。"
        )
        if number == 1:
            digest_text = (
                f"{digest_lead}用七张代表性结果图理解多样性、组成、差异、"
                "预测与因果证据的边界。"
            )
        if number == 55:
            digest_text = (
                f"{digest_lead}用七项一手研究比较不同设计的证据边界、"
                "残余偏倚与可辩护措辞。"
            )
        digest_text = source_description(project / qmd_path) or digest_text
        digest = truncate(digest_text, MAX_DIGEST_CHARS)
        (
            content,
            images,
            removed_bootstrap_blocks,
            removed_wechat_omit_blocks,
            stripped_install_calls,
            synced_topic_heading_count,
        ) = sanitize_article(
            source_html=source_html,
            article_dir=article_dir,
            source_qmd=(project / qmd_path).resolve(),
        )
        reader_reproduction = validate_reader_contract(project / qmd_path, content, project)
        if not images:
            raise RuntimeError(f"Article {number:02d} has no representative figure for its cover")
        cover = article_dir / "cover.jpg"
        cover_source = Path(images[0]["local_path"])
        create_cover(cover_source, cover)
        payload = {
            "title": title,
            "author": args.author,
            "digest": digest,
            "content": content,
            "content_source_url": args.review_url,
            "thumb_media_id": None,
            "show_cover_pic": 1,
            "need_open_comment": 0,
            "only_fans_can_comment": 0,
        }
        article_html = article_dir / "article.html"
        draft_json = article_dir / "draft.json"
        article_html.write_text(local_preview_html(title, content), encoding="utf-8")
        draft_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
        items.append(
            {
                "chapter_id": f"{number:02d}",
                "title": title,
                "title_order_prefix": f"{number}.",
                "source_qmd": str((project / qmd_path).resolve()),
                "source_html": str(source_html),
                "source_surface": source_surface,
                "source_surface_reason": source_surface_reason,
                "article_html": str(article_html),
                "draft_json": str(draft_json),
                "cover_image": str(cover),
                "cover_layout": "representative_figure_only",
                "cover_source_image": str(cover_source),
                "cover_size_bytes": cover.stat().st_size,
                "cover_sha256": sha256(cover),
                "html_chars": len(content),
                "readability": article_readability(content),
                "reader_reproduction": reader_reproduction,
                "removed_bootstrap_block_count": removed_bootstrap_blocks,
                "removed_wechat_omit_block_count": removed_wechat_omit_blocks,
                "stripped_install_call_count": stripped_install_calls,
                "synced_topic_heading_count": synced_topic_heading_count,
                "embedded_image_count": len(images),
                "embedded_images": images,
            }
        )
    errors: list[str] = []
    if len(items) != len(formal):
        errors.append(f"Expected {len(formal)} selected items, found {len(items)}")
    if len({item["title"] for item in items}) != len(items):
        errors.append("Draft titles are not unique")
    for item in items:
        if item["cover_size_bytes"] > MAX_THUMB_BYTES:
            errors.append(f"{item['chapter_id']}: cover exceeds {MAX_THUMB_BYTES} bytes")
        if item["html_chars"] < 3000:
            errors.append(f"{item['chapter_id']}: article content is unexpectedly short")
        draft = json.loads(Path(item["draft_json"]).read_text(encoding="utf-8"))
        expected_prefix = f"16S最佳实践｜{int(item['chapter_id'])}. "
        if not draft["title"].startswith(expected_prefix):
            errors.append(f"{item['chapter_id']}: title order prefix is missing")
        if len(draft["title"]) > MAX_TITLE_CHARS:
            errors.append(f"{item['chapter_id']}: title is too long")
        if len(draft["digest"]) > MAX_DIGEST_CHARS:
            errors.append(f"{item['chapter_id']}: digest is too long")
        if re.search(r"<(script|style|button|nav)\b", draft["content"], flags=re.I):
            errors.append(f"{item['chapter_id']}: unsupported HTML remains")
        if re.search(r"\{#sec-[^}]+\}", draft["content"]):
            errors.append(f"{item['chapter_id']}: unrendered section anchor remains")
        if re.search(
            r"审阅草稿|开放审阅|GitHub Draft PR|草稿箱继续查看|"
            r"header-section-number|data-local-image|"
            r"16S最佳实践[（(]\d{1,2}/\d{1,2}[）)]|"
            r"第\s*\d{1,2}\s*/\s*\d{1,2}\s*篇",
            draft["content"],
            flags=re.I,
        ):
            errors.append(f"{item['chapter_id']}: internal review or numbering metadata remains")
        if re.search(
            r"<h2[^>]*>\s*准备工作\s*</h2>|"
            r"展开：[^<]{0,80}(?:安装依赖|定义作图函数|出版级函数)|"
            r"整仓库(?:运行时|使用者|用户|的一次性|的验收器)|"
            r"只复制(?:本页|本文)|单篇复现|独立运行以上",
            draft["content"],
            flags=re.I,
        ):
            errors.append(f"{item['chapter_id']}: website-only bootstrap or maintainer prose remains")
        if re.search(
            r"<h2[^>]*>\s*(?:理论：|为什么这么做)|"
            r"隐藏决定|"
            r"(?:这里|本页|本文|本篇|我们)[^。<]{0,40}"
            r"不复制[^。<]{0,24}(?:原图|成图)",
            draft["content"],
            flags=re.I,
        ):
            errors.append(
                f"{item['chapter_id']}: template heading or figure-production narration remains"
            )
        for code_block in re.findall(
            r"<pre\b[^>]*>.*?</pre>",
            draft["content"],
            flags=re.I | re.S,
        ):
            if re.search(
                r"(?:install\.packages|BiocManager::install|"
                r"remotes::install_github|pak::pkg_install)\s*\(",
                code_block,
            ):
                errors.append(
                    f"{item['chapter_id']}: package-install command remains in WeChat code"
                )
                break
        if "/pull/" in str(draft.get("content_source_url", "")):
            errors.append(f"{item['chapter_id']}: content_source_url points to a pull request")
        for image_record in item["embedded_images"]:
            if image_record["size_bytes"] > MAX_ARTICLE_IMAGE_BYTES:
                errors.append(f"{item['chapter_id']}: article image exceeds upload budget")
    report = {
        "generated_at": utc_now(),
        "status": "failed" if errors else "passed",
        "project_root": str(project),
        "manifest": str(manifest_path),
        "manifest_sha256": sha256(manifest_path),
        "qa_report": str(qa_path),
        "qa_status": qa_report.get("status"),
        "qa_run_key": qa_report.get("run_key"),
        "qa_manifest_hash": qa_report.get("manifest_hash"),
        "formal_count": args.formal_count,
        "selected_chapters": [int(chapter["number"]) for chapter in formal],
        "author": args.author,
        "title_style": "series_then_ordinal_dot",
        "review_url": args.review_url,
        "item_count": len(items),
        "fallback_item_count": sum(
            item["source_surface"] == "verified_fallback_bundle"
            for item in items
        ),
        "fallback_chapters": [
            item["chapter_id"]
            for item in items
            if item["source_surface"] == "verified_fallback_bundle"
        ],
        "embedded_image_count": sum(item["embedded_image_count"] for item in items),
        "errors": errors,
        "items": items,
    }
    report_path = output_dir / "report.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    if errors:
        raise RuntimeError("; ".join(errors))
    return report


def main() -> int:
    report = build(parse_args())
    print(
        json.dumps(
            {
                "status": report["status"],
                "item_count": report["item_count"],
                "embedded_image_count": report["embedded_image_count"],
                "review_url": report["review_url"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
