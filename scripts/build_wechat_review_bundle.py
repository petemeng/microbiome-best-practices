#!/usr/bin/env python3
"""Build WeChat-ready review articles from the validated Quarto HTML site.

The script is intentionally offline: it sanitizes and inlines styles, optimizes
article images, creates deterministic covers, and writes resumable JSON payloads.
Uploading images and creating official-account drafts remain separate actions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import urljoin

import yaml
from lxml import etree, html
from PIL import Image, ImageDraw, ImageFont, ImageOps


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
    parser.add_argument("--formal-count", type=int, default=55)
    parser.add_argument("--author", default="Songlab")
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


def class_tokens(element: etree._Element) -> set[str]:
    return set((element.get("class") or "").split())


def site_html_path(site_dir: Path, qmd_path: str) -> Path:
    qmd = Path(qmd_path)
    if qmd.name == "index.qmd":
        return site_dir / "index.html"
    return site_dir / qmd.with_suffix(".html")


def _font_path(family: str) -> str:
    result = subprocess.run(
        ["fc-match", "-f", "%{file}\n", family],
        check=True,
        capture_output=True,
        text=True,
    )
    path = result.stdout.strip().splitlines()[0]
    if not path:
        raise RuntimeError(f"No font found for {family}")
    return path


def _wrap_cjk(
    text: str,
    draw: ImageDraw.ImageDraw,
    font: ImageFont.FreeTypeFont,
    max_width: int = 610,
    max_lines: int = 3,
) -> list[str]:
    """Wrap mixed Chinese/Latin titles without splitting technical terms."""

    def wrap_segment(segment: str) -> list[str]:
        normalized = re.sub(r"\s+", " ", segment).strip()
        tokens = re.findall(
            r"[A-Za-z0-9]+(?:[./+–—-][A-Za-z0-9]+)*|\s+|.",
            normalized,
        )
        output: list[str] = []
        current = ""
        for token in tokens:
            candidate = (current + token).lstrip()
            if current and draw.textlength(candidate, font=font) > max_width:
                output.append(current.rstrip())
                current = token.lstrip()
            else:
                current = candidate
        if current.strip():
            output.append(current.rstrip())
        return output

    normalized = re.sub(r"\s+", " ", text).strip()
    if "：" in normalized:
        head, tail = normalized.split("：", 1)
        lines = wrap_segment(head) + wrap_segment(tail)
    else:
        lines = wrap_segment(normalized)
    if len(lines) > max_lines:
        lines = lines[:max_lines]
        while lines[-1] and draw.textlength(lines[-1] + "…", font=font) > max_width:
            lines[-1] = lines[-1][:-1].rstrip()
        lines[-1] = lines[-1] + "…"
    return lines or ["16S 微生物组最佳实践"]


def create_cover(raw_title: str, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    image = Image.new("RGB", (900, 383), "#f7f1e6")
    draw = ImageDraw.Draw(image)
    for y in range(image.height):
        ratio = y / max(1, image.height - 1)
        start = (247, 241, 230)
        end = (238, 244, 234)
        colour = tuple(round(a + (b - a) * ratio) for a, b in zip(start, end))
        draw.line((0, y, image.width, y), fill=colour)
    draw.rectangle((0, 0, 900, 14), fill="#7c9970")
    draw.ellipse((735, 52, 833, 150), fill="#e1b36a")
    serif = ImageFont.truetype(_font_path("Noto Serif CJK SC"), 44)
    sans = ImageFont.truetype(_font_path("Noto Sans CJK SC"), 22)
    small = ImageFont.truetype(_font_path("Noto Sans CJK SC"), 19)
    lines = _wrap_cjk(raw_title, draw, serif)
    draw.multiline_text((54, 76), "\n".join(lines), font=serif, fill="#203124", spacing=10)
    draw.rounded_rectangle((54, 300, 330, 344), radius=14, fill="#f9f4ea")
    draw.text((72, 307), "MICROBIOME TUTORIAL", font=sans, fill="#6f8561")
    draw.text((690, 326), "16S · BEST PRACTICES", font=small, fill="#6f8561")
    for quality in (82, 74, 68, 62, 56, 50, 44):
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


def flatten_code(main: etree._Element) -> None:
    for pre in main.xpath(".//pre"):
        code_text = "".join(pre.itertext()).rstrip()
        for child in list(pre):
            pre.remove(child)
        pre.text = code_text
        pre.set("style", STYLES["pre"])


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
        record = {
            "source_path": str(source),
            "local_path": str(destination.resolve()),
            "relative_src": relative,
            "sha256": sha256(destination),
            "size_bytes": destination.stat().st_size,
            "width": Image.open(destination).width,
            "height": Image.open(destination).height,
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
) -> tuple[str, list[dict[str, Any]]]:
    document = html.parse(str(source_html)).getroot()
    mains = document.xpath(
        '//main[contains(concat(" ",normalize-space(@class)," ")," content ")]'
        ' | //main[@id="quarto-document-content"] | //main'
    )
    if not mains:
        raise RuntimeError(f"No main content found in {source_html}")
    main = deepcopy(mains[0])
    main.tag = "section"
    remove_unwanted(main)
    flatten_code(main)
    transform_special_blocks(main)
    apply_inline_styles(main)
    images = resolve_and_optimize_images(main, source_html, article_dir)
    strip_unsupported_attributes(main)
    main.set("style", ROOT_STYLE)
    return etree.tostring(main, encoding="unicode", method="html"), images


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
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    qa_report = json.loads(qa_path.read_text(encoding="utf-8"))
    if qa_report.get("status") != "passed":
        raise RuntimeError("qa_report.json.status must be passed before draft generation")
    chapters = manifest.get("series", {}).get("chapters", [])
    if len(chapters) != 55:
        raise RuntimeError(f"Expected 55 manifest chapters, found {len(chapters)}")
    formal = chapters[: args.formal_count]
    output_dir.mkdir(parents=True, exist_ok=True)
    items: list[dict[str, Any]] = []
    for chapter in formal:
        number = int(chapter["number"])
        raw_title = str(chapter["title"])
        qmd_path = str(chapter["file"])
        source_html = site_html_path(site_dir, qmd_path)
        if not source_html.exists():
            raise FileNotFoundError(source_html)
        article_dir = output_dir / f"{number:02d}"
        article_dir.mkdir(parents=True, exist_ok=True)
        title = truncate(f"16S最佳实践｜{raw_title}", MAX_TITLE_CHARS)
        digest_lead = raw_title if raw_title.endswith(("。", "！", "？", "!", "?")) else f"{raw_title}。"
        digest_text = (
            f"{digest_lead}真实数据、完整复现代码、结果解释与发表级重绘图。"
        )
        if number == 55:
            digest_text = (
                f"{digest_lead}用七项一手研究比较不同设计的证据边界、"
                "残余偏倚与可辩护措辞。"
            )
        digest = truncate(digest_text, MAX_DIGEST_CHARS)
        cover = article_dir / "cover.jpg"
        create_cover(raw_title, cover)
        content, images = sanitize_article(
            source_html=source_html,
            article_dir=article_dir,
        )
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
                "source_qmd": str((project / qmd_path).resolve()),
                "source_html": str(source_html),
                "article_html": str(article_html),
                "draft_json": str(draft_json),
                "cover_image": str(cover),
                "cover_size_bytes": cover.stat().st_size,
                "cover_sha256": sha256(cover),
                "html_chars": len(content),
                "embedded_image_count": len(images),
                "embedded_images": images,
            }
        )
    errors: list[str] = []
    if len(items) != args.formal_count:
        errors.append(f"Expected {args.formal_count} items, found {len(items)}")
    if len({item["title"] for item in items}) != len(items):
        errors.append("Draft titles are not unique")
    for item in items:
        if item["cover_size_bytes"] > MAX_THUMB_BYTES:
            errors.append(f"{item['chapter_id']}: cover exceeds {MAX_THUMB_BYTES} bytes")
        if item["html_chars"] < 3000:
            errors.append(f"{item['chapter_id']}: article content is unexpectedly short")
        draft = json.loads(Path(item["draft_json"]).read_text(encoding="utf-8"))
        if len(draft["title"]) > MAX_TITLE_CHARS:
            errors.append(f"{item['chapter_id']}: title is too long")
        if len(draft["digest"]) > MAX_DIGEST_CHARS:
            errors.append(f"{item['chapter_id']}: digest is too long")
        if re.search(r"<(script|style|button|nav)\b", draft["content"], flags=re.I):
            errors.append(f"{item['chapter_id']}: unsupported HTML remains")
        if re.search(
            r"审阅草稿|开放审阅|GitHub Draft PR|草稿箱继续查看|"
            r"header-section-number|data-local-image|"
            r"16S最佳实践[（(]\d{1,2}/\d{1,2}[）)]|"
            r"第\s*\d{1,2}\s*/\s*\d{1,2}\s*篇",
            draft["content"],
            flags=re.I,
        ):
            errors.append(f"{item['chapter_id']}: internal review or numbering metadata remains")
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
        "review_url": args.review_url,
        "item_count": len(items),
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
