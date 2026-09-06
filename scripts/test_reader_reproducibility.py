import unittest
from lxml import html
from reader_reproducibility import coalesce_rendered_code, public_blocks

SOURCE = '''---
reader-reproduction:
  required: true
---
```{r}
#| label: reader-demo
x <- 1
x
save_plot(x)
```
'''


class RenderedCodeTests(unittest.TestCase):
    def test_join_verified_plot_split(self):
        main = html.fromstring('<main><div class="cell"><pre class="sourceCode r">x &lt;- 1\nx</pre><div class="cell-output">result</div><pre class="sourceCode r">save_plot(x)</pre></div></main>')
        self.assertEqual(coalesce_rendered_code(SOURCE, main), 1)
        self.assertEqual(len(public_blocks(SOURCE, html.tostring(main, encoding='unicode'))), 1)
        self.assertIn('result', main.text_content())

    def test_never_repair_missing_code(self):
        main = html.fromstring('<main><div class="cell"><pre class="sourceCode r">x &lt;- 1</pre><pre class="sourceCode r">save_plot(x)</pre></div></main>')
        self.assertEqual(coalesce_rendered_code(SOURCE, main), 0)
        with self.assertRaises(ValueError):
            public_blocks(SOURCE, html.tostring(main, encoding='unicode'))

    def test_reject_stale_computation(self):
        main = html.fromstring('<main><div class="cell"><pre class="sourceCode r">x &lt;- 2\nx</pre><pre class="sourceCode r">save_plot(x)</pre></div></main>')
        self.assertEqual(coalesce_rendered_code(SOURCE, main), 0)
        with self.assertRaises(ValueError):
            public_blocks(SOURCE, html.tostring(main, encoding='unicode'))


if __name__ == '__main__':
    unittest.main()
