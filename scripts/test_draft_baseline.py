import copy
import unittest
from sync_wechat_revision import SafeError, old_article_signature, verify_old_article


class DraftBaselineTests(unittest.TestCase):
    def setUp(self):
        self.article = dict(title='16S最佳实践｜1. Topic', author='Peter', digest='Summary',
                            content_source_url='https://example.org/source', thumb_media_id='test-cover',
                            content='<p>Observed result</p><img src="http://mmbiz.qpic.cn/example/0">')

    def test_display_size_and_whitespace_rewriting_is_allowed(self):
        remote = copy.deepcopy(self.article)
        remote['content'] = '<p>Observed   result</p><img src="https://mmbiz.qpic.cn/example/640">'
        verify_old_article({'old_article_signature': old_article_signature(self.article)}, remote)

    def test_edited_old_body_is_retained(self):
        remote = copy.deepcopy(self.article)
        remote['content'] = '<p>User edited result</p>'
        with self.assertRaises(SafeError):
            verify_old_article({'old_article_signature': old_article_signature(self.article)}, remote)

    def test_edited_metadata_is_retained(self):
        remote = dict(self.article, digest='User edited summary')
        with self.assertRaises(SafeError):
            verify_old_article({'old_article_signature': old_article_signature(self.article)}, remote)


if __name__ == '__main__':
    unittest.main()
