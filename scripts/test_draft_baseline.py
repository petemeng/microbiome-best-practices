import copy
import unittest
from sync_wechat_revision import SafeError, old_article_signature, verify_old_article, mapped_article_payload, body_signature


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

    def test_current_asset_ledger_rebinds_a_stale_archive(self):
        entry={'cover_media_id':'current-test-cover','body_images':[{'url':'https://mmbiz.qpic.cn/current/0'}]}
        rebound=mapped_article_payload(self.article,entry)
        self.assertEqual(rebound['thumb_media_id'],'current-test-cover')
        self.assertEqual(body_signature(rebound['content']),('Observed result',['https://mmbiz.qpic.cn/current/0']))
        self.assertEqual(self.article['thumb_media_id'],'test-cover')

    def test_ambiguous_archive_asset_count_is_rejected(self):
        with self.assertRaises(SafeError):
            mapped_article_payload(self.article,{'cover_media_id':'test-cover','body_images':[]})


if __name__ == '__main__':
    unittest.main()
