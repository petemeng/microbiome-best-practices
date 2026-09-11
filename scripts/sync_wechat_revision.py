#!/usr/bin/env python3
"""Replace explicitly mapped drafts and verify the full body before deleting old ones.

Only draft and image-material endpoints are implemented. Publishing and mass
sending are deliberately unavailable. Credentials and resumable journals remain
in private files and are never printed. Run without --execute for read-only
identity checks, or use --verify-existing to inspect already-created drafts.
"""
from __future__ import annotations
import argparse,copy,hashlib,html as stdhtml,json,os,re,threading
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from urllib.parse import urlsplit

import requests
from lxml import html
from finalize_wechat_revision import atomic_json

class SafeError(RuntimeError):
    pass

ALLOWED={'draft/get','draft/add','draft/delete','draft/batchget','draft/count',
         'media/uploadimg','material/add_material'}

class Client:
    def __init__(self,token_file):
        self.token_file=token_file
        self.session=requests.Session()
    def api(self,endpoint,payload=None,*,file=None,media_type=None):
        if endpoint not in ALLOWED:raise SafeError('Endpoint is outside draft synchronization scope')
        if self.token_file.stat().st_mode & 0o077:raise SafeError('Credential file must be private (mode 600)')
        token=json.loads(self.token_file.read_text())['access_token']
        params={'access_token':token}
        if media_type:params['type']=media_type
        url='https://api.weixin.qq.com/cgi-bin/'+endpoint
        try:
            if file:
                with file.open('rb') as stream:
                    response=self.session.post(url,params=params,files={'media':(file.name,stream,'image/jpeg')},timeout=(15,90))
            elif endpoint=='draft/count':
                response=self.session.get(url,params=params,timeout=(15,90))
            else:
                response=self.session.post(url,params=params,data=json.dumps(payload,ensure_ascii=False).encode(),
                    headers={'Content-Type':'application/json; charset=utf-8'},timeout=(15,90))
            response.raise_for_status()
            response.encoding='utf-8'
            result=response.json()
        except (requests.RequestException,ValueError):
            raise SafeError(f'{endpoint}: request did not return a confirmed JSON response; inspect the journal before retrying a write') from None
        if result.get('errcode',0)!=0:
            raise SafeError(f'{endpoint}: WeChat error code {int(result["errcode"])}')
        return result

def article_from_response(response):
    items=response.get('news_item',[])
    if len(items)!=1:raise SafeError('Expected exactly one article in the mapped draft')
    return items[0]

def body_signature(content):
    tree=html.fromstring(content)
    text=' '.join(stdhtml.unescape(tree.text_content()).split())
    images=[i.get('src') or i.get('data-src') for i in tree.xpath('//img')]
    return text,images

def image_identity(url):
    """WeChat rewrites upload URLs from HTTP /0 to HTTPS /640 for display."""
    parsed=urlsplit(url or '')
    path=parsed.path
    if path.rsplit('/',1)[-1] in ('0','640'):
        path=path.rsplit('/',1)[0]+'/{display-size}'
    return parsed.hostname,path

def verify_article(expected,actual):
    for key in ('title','author','digest','content_source_url','thumb_media_id'):
        if actual.get(key)!=expected.get(key):raise SafeError(f'Remote {key} differs from the reviewed article')
    expected_text,expected_images=body_signature(expected['content'])
    actual_text,actual_images=body_signature(actual.get('content',''))
    if expected_text!=actual_text or list(map(image_identity,expected_images))!=list(map(image_identity,actual_images)):
        raise SafeError('Remote body text or images differ from the reviewed article; old draft retained')

def cdn_url(url):
    parsed=urlsplit(url)
    return parsed.scheme in ('https','http') and (parsed.hostname or '').endswith('.qpic.cn')

def validate_plan(plan):
    keys=[e['chapter_id'] for e in plan['entries']]
    ids=[e['old_draft_media_id'] for e in plan['entries']]
    if len(set(keys))!=len(keys) or len(set(ids))!=len(ids):raise SafeError('Plan contains duplicate draft targets')
    for e in plan['entries']:
        if not e['title'].startswith(f'16S最佳实践｜{int(e["chapter_id"])}. ') or e['draft']['author']!='Peter':
            raise SafeError('Unexpected series identity')
        if not re.fullmatch('[0-9a-f]{40}',e['source_commit']):raise SafeError('Immutable source commit required')

def synchronize(entry,args):
    key=entry['chapter_id'];client=Client(args.token_file)
    journal_path=args.journal_directory/f'{args.journal_prefix}-{key}.json'
    state=json.loads(journal_path.read_text()) if journal_path.exists() else {k:copy.deepcopy(v) for k,v in entry.items() if k!='draft'}
    if state['old_draft_media_id']!=entry['old_draft_media_id'] or state['source_commit']!=entry['source_commit']:
        raise SafeError('Existing journal belongs to a different revision')
    def save():atomic_json(journal_path,state)
    if not state.get('created'):
        if args.verify_existing:raise SafeError('No created draft is available for verification')
        old=article_from_response(client.api('draft/get',{'media_id':entry['old_draft_media_id']}))
        if not old.get('title','').startswith(f'16S最佳实践｜{int(key)}. ') or old.get('author')!='Peter':
            raise SafeError('Old mapped draft identity does not match this chapter')
        if not args.execute:return {'chapter':key,'identity_verified':True,'mutated':False}
        state['old_identity_verified']=True;save()
    if not args.execute and not args.verify_existing:
        return {'chapter':key,'created_previously':True,'mutated':False}
    if args.verify_existing and not state.get('created'):
        raise SafeError('No created draft is available for verification')
    if not state.get('cover_media_id'):
        if not args.execute:raise SafeError('Missing cover in verification-only mode')
        state['cover_media_id']=client.api('material/add_material',file=Path(state['cover_file']),media_type='thumb')['media_id'];save()
    for image in state['images']:
        path=Path(image['file'])
        if hashlib.sha256(path.read_bytes()).hexdigest()!=image['sha256']:raise SafeError('Reviewed body image changed')
        if not image.get('url'):
            if not args.execute:raise SafeError('Missing body image in verification-only mode')
            image['url']=client.api('media/uploadimg',file=path)['url'];save()
        if not cdn_url(image['url']):raise SafeError('Body image is not on WeChat image CDN')
    article=copy.deepcopy(entry['draft'])
    article['thumb_media_id']=state['cover_media_id']
    article['content_source_url']=entry['source_url']
    article['show_cover_pic']=1
    for image in state['images']:
        article['content']=article['content'].replace('src="'+image['src']+'"','src="'+image['url']+'"')
    if any(not cdn_url(url or '') for url in body_signature(article['content'])[1]):raise SafeError('Unresolved article image')
    if not state.get('created'):
        if state.get('creation_pending'):
            raise SafeError('An earlier creation is unconfirmed; reconcile it before another add request')
        state['creation_pending']=True;save()
        result=client.api('draft/add',{'articles':[article]})
        state.update(draft_media_id=result['media_id'],created=True,creation_pending=False);save()
    if state['draft_media_id']==state['old_draft_media_id']:raise SafeError('New and old draft must be distinct')
    remote=article_from_response(client.api('draft/get',{'media_id':state['draft_media_id']}))
    verify_article(article,remote)
    state.update(verified=True,full_body_verified=True,
        content_sha256=hashlib.sha256(article['content'].encode()).hexdigest(),
        remote_content_sha256=hashlib.sha256(remote['content'].encode()).hexdigest(),
        verification_scope='Full remote body text and image asset identities (allowing WeChat HTTPS/display-size URL rewriting), title, author, digest, source URL and cover ID matched against the reviewed submission.')
    save()
    if args.execute and not state.get('old_draft_deleted'):
        old=article_from_response(client.api('draft/get',{'media_id':state['old_draft_media_id']}))
        if not old.get('title','').startswith(f'16S最佳实践｜{int(key)}. ') or old.get('author')!='Peter':
            raise SafeError('Old draft identity changed before deletion')
        client.api('draft/delete',{'media_id':state['old_draft_media_id']})
        state['old_draft_deleted']=True;save()
    return {'chapter':key,'created':state['created'],'full_body_verified':True,'old_deleted':state.get('old_draft_deleted',False)}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--token-file',type=Path,required=True)
    p.add_argument('--journal-directory',type=Path,required=True)
    p.add_argument('--journal-prefix',required=True)
    p.add_argument('--jobs',type=int,default=3)
    mode=p.add_mutually_exclusive_group()
    mode.add_argument('--execute',action='store_true')
    mode.add_argument('--verify-existing',action='store_true')
    args=p.parse_args()
    args.journal_directory.mkdir(mode=0o700,parents=True,exist_ok=True)
    plan=json.loads(args.plan.read_text());validate_plan(plan)
    errors=[]
    with ThreadPoolExecutor(max_workers=max(1,min(args.jobs,4))) as pool:
        futures={pool.submit(synchronize,e,args):e['chapter_id'] for e in plan['entries']}
        for future in as_completed(futures):
            key=futures[future]
            try:print(json.dumps(future.result()),flush=True)
            except Exception as error:
                message=str(error) if isinstance(error,SafeError) else type(error).__name__
                errors.append(key);print(json.dumps({'chapter':key,'status':'failed','error':message}),flush=True)
    print(json.dumps({'status':'failed' if errors else 'passed','chapters':len(futures),'failed':errors}),flush=True)
    return bool(errors)

if __name__=='__main__':raise SystemExit(main())
