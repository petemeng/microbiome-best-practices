#!/usr/bin/env python3
"""Offline regressions for scoped, reader-facing WeChat bundles."""

import tempfile
import unittest
from html import escape
from pathlib import Path

from lxml import html
from PIL import Image

from reader_reproducibility import (
    public_blocks, reader_blocks, script_text, validate_reader_contract,
)

from build_wechat_review_bundle import (
    article_readability, localize_figure_labels, remove_unwanted,
    sanitize_article, select_chapters, source_description,
)


class ChapterSelectionTests(unittest.TestCase):
    def setUp(self):
        self.chapters = [{"number": number} for number in range(1, 56)]

    def test_default_preserves_full_series(self):
        self.assertEqual(select_chapters(self.chapters, 55, None), self.chapters)

    def test_selected_chapters_follow_manifest_order(self):
        self.assertEqual(select_chapters(self.chapters, 55, [55, 26]),
                         [{"number": 26}, {"number": 55}])

    def test_rejects_empty_duplicate_missing_and_out_of_scope(self):
        for selection in ([], [26, 26], [0], [56]):
            with self.subTest(selection=selection), self.assertRaises(ValueError):
                select_chapters(self.chapters, 55, selection)
        with self.assertRaises(ValueError):
            select_chapters(self.chapters, 25, [26])

    def test_rejects_invalid_formal_count(self):
        for count in (0, 56):
            with self.subTest(count=count), self.assertRaises(ValueError):
                select_chapters(self.chapters, count, None)


class EditorialSurfaceTests(unittest.TestCase):
    def test_hidden_callout_label_is_removed_without_losing_title(self):
        document = html.fromstring('<section><div><span class="screen-reader-only">提示</span>按研究问题选图</div></section>')
        remove_unwanted(document)
        self.assertEqual(document.text_content(), "按研究问题选图")

    def test_figure_numbers_and_local_references_are_article_local(self):
        document = html.fromstring('''<section>
            <a href="#fig-first">图 26.1</a><a href="#fig-last">图 26.4</a>
            <a href="other.html#fig-other">图 27.1</a>
            <figure><figcaption>图 26.1: First result</figcaption></figure>
            <figure><figcaption>图 26.4: Last result <em>with emphasis</em></figcaption></figure>
            </section>''')
        localize_figure_labels(document)
        self.assertEqual(document.xpath(".//a")[0].text, "图 1")
        self.assertEqual(document.xpath(".//a")[1].text, "图 2")
        self.assertEqual(document.xpath(".//a")[2].text, "图 27.1")
        self.assertEqual(document.xpath(".//figcaption")[1].text, "图 2: Last result ")
        self.assertEqual(len(document.xpath(".//figcaption/em")), 1)

    def test_explicit_omission_keeps_adjacent_figure_and_reasoning(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            qmd = root / "chapter.qmd"
            qmd.write_text('---\ntitle: Composition\n---\n', encoding="utf-8")
            Image.new("RGB", (100, 80), "white").save(root / "result.png")
            source = root / "chapter.html"
            source.write_text('''<html><body><main>
                <p>Use each sample's total count as the denominator.</p>
                <details class="wechat-omit"><summary>Full plotting code</summary>
                <pre>save_pub(plot)</pre><pre>Internal export status</pre></details>
                <pre>relative = counts / sample_total</pre>
                <figure><img src="result.png"><figcaption>Observed composition</figcaption></figure>
                <p>The difference is descriptive, not a significance test.</p>
                </main></body></html>''', encoding="utf-8")
            article_dir = root / "article"
            article_dir.mkdir()
            content, images, _, omitted, _, _ = sanitize_article(source, article_dir, qmd)
            document = html.fromstring(content)
            self.assertEqual(omitted, 1)
            self.assertEqual(len(images), 1)
            self.assertEqual(len(document.xpath(".//img")), 1)
            self.assertEqual(len(document.xpath(".//pre")), 1)
            self.assertIn("not a significance test", document.text_content())
            self.assertNotIn("save_pub", content)
            self.assertNotIn("Internal export status", content)

    def test_description_is_optional_and_has_no_installation_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            qmd = Path(directory) / "chapter.qmd"
            qmd.write_text('---\ntitle: Composition\ndescription: Compare samples first.\n---\nBody', encoding="utf-8")
            self.assertEqual(source_description(qmd), "Compare samples first.")
            qmd.write_text('---\ntitle: Composition\n---\nBody', encoding="utf-8")
            self.assertIsNone(source_description(qmd))

    def test_reading_metrics_measure_code_and_first_figure(self):
        metrics = article_readability('<section><p>hello</p><pre>abc</pre><img src="x"><p>world</p></section>')
        self.assertEqual(metrics["preformatted_block_count"], 1)
        self.assertEqual(metrics["preformatted_chars"], 3)
        self.assertEqual(metrics["chars_before_first_image"], 8)
        self.assertIsNone(article_readability("<p>text</p>")["chars_before_first_image"])


class ReaderReproductionTests(unittest.TestCase):
    source = '''---
title: Composition
reader-reproduction:
  required: true
  script: examples/analysis.R
---
```{r}
#| label: reader-download
x <- read.delim("counts.tsv")
```
```{r}
#| label: fig-result
#| fig-cap: Result
plot(x)
```
'''

    def surface(self, codes):
        return "<section>" + "".join("<pre>" + escape(code) + "</pre>" for code in codes) + "</section>"

    def test_source_public_code_and_script_agree(self):
        blocks = reader_blocks(self.source)
        content = self.surface([blocks[0][1], "console output", blocks[1][1]])
        self.assertEqual(public_blocks(self.source, content), blocks)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            qmd = root / "chapter.qmd"
            qmd.write_text(self.source)
            (root / "examples").mkdir()
            script = root / "examples/analysis.R"
            script.write_text(script_text(self.source, blocks))
            self.assertEqual(validate_reader_contract(qmd, content, root)["code_blocks"], 2)
            script.write_text("plot(x)\n")
            with self.assertRaisesRegex(ValueError, "script differs"):
                validate_reader_contract(qmd, content, root)

    def test_rejects_missing_changed_and_out_of_order_code(self):
        codes = [code for _, code in reader_blocks(self.source)]
        for content in (self.surface(codes[1:]), self.surface(codes[::-1]),
                        self.surface([codes[0], "plot(y)"])):
            with self.subTest(content=content), self.assertRaises(ValueError):
                public_blocks(self.source, content)

    def test_rejects_hidden_or_unlabelled_required_code(self):
        for source in (self.source.replace("#| label: reader-download", "#| echo: false\n#| label: reader-download"),
                       self.source.replace("#| label: reader-download\n", "")):
            with self.assertRaises(ValueError):
                reader_blocks(source)

    def test_existing_labels_and_chapter_specific_packages(self):
        source = self.source.replace('  required: true',
            '  required: true\n  label-policy: all-executed\n  packages: [vegan, ggplot2]')
        source = source.replace('reader-download', 'read-study-inputs')
        blocks = reader_blocks(source)
        self.assertEqual(blocks[0][0], 'read-study-inputs')
        self.assertIn('Required packages: vegan, ggplot2.', script_text(source, blocks))
        with self.assertRaises(ValueError):
            reader_blocks(source.replace('#| label: read-study-inputs\n', ''))

    def test_sanitizer_cannot_remove_required_download(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            qmd = root / "chapter.qmd"
            qmd.write_text(self.source)
            blocks = reader_blocks(self.source)
            source_html = root / "chapter.html"
            source_html.write_text('<main><details class="wechat-omit">' +
                self.surface([blocks[0][1]]) + '</details>' + self.surface([blocks[1][1]]) + '</main>')
            with self.assertRaisesRegex(ValueError, "reader-download"):
                sanitize_article(source_html, root, qmd)


if __name__ == "__main__":
    unittest.main()
