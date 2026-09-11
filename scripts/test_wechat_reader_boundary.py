import re,unittest
from lxml import html
from reader_reproducibility import decisive_code,omit_generic_reader_helpers,public_blocks,reader_blocks,script_text
from sync_wechat_revision import SafeError,verify_article,validate_plan
from validate_series import PROHIBITED_PUBLIC_REGEXES

class ReaderBoundaryTests(unittest.TestCase):
    def test_figure_production_narration_variants_are_rejected(self):
        for text in ('这里不嵌入原图', '以下图不复制原论文图', '这里不拼贴论文原图'):
            with self.subTest(text=text):
                self.assertTrue(any(re.search(pattern,text) for pattern in PROHIBITED_PUBLIC_REGEXES))
    def test_only_generic_top_level_helpers_are_omitted(self):
        code='library(ggplot2)\n\nfont_pub <- "sans"\ntheme_pub <- function() {\n  ggplot2::theme_bw(base_family = font_pub)\n}\n\notutab <- read.delim("otutab.tsv")\nclean_taxon <- function(x) trimws(x)\nggplot2::theme_set(theme_pub())'
        visible=decisive_code(code)
        self.assertNotIn('theme_pub <-',visible)
        self.assertNotIn('font_pub <-',visible)
        self.assertIn('otutab <-',visible)
        self.assertIn('clean_taxon <-',visible)
        self.assertIn('ggplot2::theme_set(theme_pub())',visible)

    def test_full_script_is_not_changed_when_public_bootstrap_is_omitted(self):
        code='theme_pub <- function() 1\n\nx <- 2\nprint(x)'
        source='---\nreader-reproduction:\n  required: true\n---\n```{r}\n#| label: reader-data\n'+code+'\n```\n'
        before=script_text(source,reader_blocks(source))
        tree=html.fromstring('<section><pre></pre></section>')
        tree.xpath('.//pre')[0].text=code
        self.assertEqual(omit_generic_reader_helpers(source,tree),1)
        public=html.tostring(tree,encoding='unicode')
        self.assertEqual(public_blocks(source,public),[('reader-data','x <- 2\nprint(x)')])
        self.assertEqual(before,script_text(source,reader_blocks(source)))

    def test_shared_statement_line_is_not_silently_removed(self):
        with self.assertRaises(ValueError):decisive_code('theme_pub <- function() 1; x <- 2')

class RemoteReadbackTests(unittest.TestCase):
    def article(self):
        return dict(title='16S最佳实践｜2. 范围',author='Peter',content_source_url='https://example.org/source',
                    thumb_media_id='test-cover',content='<section><p>真实结果</p><pre>x &lt;- 1</pre><img src="https://mmbiz.qpic.cn/a"/></section>')
    def test_equivalent_html_text_and_images_are_accepted(self):
        expected=self.article();actual=dict(expected,content=expected['content'].replace('<p>','<p class="wechat">'))
        verify_article(expected,actual)
    def test_changed_code_is_rejected(self):
        expected=self.article();actual=dict(expected,content=expected['content'].replace('x &lt;- 1','x &lt;- 2'))
        with self.assertRaises(SafeError):verify_article(expected,actual)
    def test_wechat_display_size_rewrite_is_accepted_but_another_image_is_not(self):
        expected=self.article();expected['content']=expected['content'].replace('/a"','/a/0?from=appmsg"')
        actual=dict(expected,content=expected['content'].replace('/a/0?from=appmsg','/a/640?from=appmsg'))
        verify_article(expected,actual)
        with self.assertRaises(SafeError):verify_article(expected,dict(actual,content=actual['content'].replace('/a/640','/b/640')))
    def test_changed_cover_is_rejected(self):
        expected=self.article()
        with self.assertRaises(SafeError):verify_article(expected,dict(expected,thumb_media_id='different'))
    def test_changed_digest_is_rejected(self):
        expected=self.article();expected['digest']='真实数据与实际结果'
        with self.assertRaises(SafeError):verify_article(expected,dict(expected,digest='其他内容'))
    def test_duplicate_targets_are_rejected(self):
        entry=dict(chapter_id='02',old_draft_media_id='test-old',title='16S最佳实践｜2. 范围',draft={'author':'Peter'},source_commit='a'*40)
        with self.assertRaises(SafeError):validate_plan({'entries':[entry,entry]})

if __name__=='__main__':unittest.main()
