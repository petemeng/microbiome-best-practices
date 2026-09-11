#!/usr/bin/env python3
"""Prepare a private, resumable revision plan; never contact WeChat or delete a draft."""
import argparse,hashlib,json,os,re
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--project-root',type=Path,default=Path(__file__).resolve().parents[1])
p.add_argument('--bundle',type=Path,required=True)
p.add_argument('--live-map',type=Path,required=True)
p.add_argument('--source-commit',required=True)
p.add_argument('--output',type=Path,required=True)
args=p.parse_args()
if not re.fullmatch('[0-9a-f]{40}',args.source_commit):raise SystemExit('A full immutable commit is required')
root=args.project_root.resolve();bundle=args.bundle.resolve()
live=json.loads(args.live_map.read_text());old_bundle=Path(live['source_bundle'])
report=json.loads((bundle/'report.json').read_text())
if report['status']!='passed':raise SystemExit('Only a passed bundle can be planned')
mapped={e['chapter_id']:e for e in live['entries']}
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
entries=[]
for item in report['items']:
    key=item['chapter_id'];old=mapped[key]
    if old['author']!='Peter' or not item['title'].startswith(f'16S最佳实践｜{int(key)}. '):raise ValueError('Unexpected series identity')
    draft=json.loads(Path(item['draft_json']).read_text())
    reuse={}
    for image in old['body_images']:
        candidate=(old_bundle/key/image['src']).resolve()
        if candidate.is_relative_to(old_bundle.resolve()) and candidate.is_file():reuse[sha(candidate)]=image['url']
    images=[]
    for image in item['embedded_images']:
        path=Path(image['local_path']).resolve()
        if not path.is_relative_to(bundle):raise ValueError('Upload outside the reviewed bundle')
        digest=sha(path)
        if digest!=image['sha256']:raise ValueError('Reviewed image changed')
        images.append({'src':image['relative_src'],'file':str(path),'sha256':digest,
            'url':reuse.get(digest),'resolution':[image['width'],image['height']]})
    cover=Path(item['cover_image']).resolve()
    if not cover.is_relative_to(bundle):raise ValueError('Cover outside the reviewed bundle')
    old_cover=old_bundle/key/'cover.jpg'
    cover_id=old['cover_media_id'] if old_cover.is_file() and sha(old_cover)==sha(cover) else None
    relative=Path(item['source_qmd']).resolve().relative_to(root).as_posix()
    source_url=f'https://github.com/petemeng/microbiome-best-practices/blob/{args.source_commit}/{relative}'
    draft['content_source_url']=source_url
    entries.append({'chapter_id':key,'title':item['title'],'old_draft_media_id':old['draft_media_id'],
        'cover_media_id':cover_id,'cover_file':str(cover),'images':images,'draft':draft,
        'source_commit':args.source_commit,'source_url':source_url,'created':False,'verified':False,'old_draft_deleted':False})
plan={'baseline_map_sha256':sha(args.live_map),'live_map':str(args.live_map.resolve()),
    'bundle':str(bundle),'source_commit':args.source_commit,'entries':entries}
args.output.parent.mkdir(parents=True,exist_ok=True)
fd=os.open(args.output,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
with os.fdopen(fd,'w') as out:json.dump(plan,out,ensure_ascii=False,indent=2);out.write('\n')
print(json.dumps({'chapters':[e['chapter_id'] for e in entries],
    'body_images_reused':sum(bool(i['url']) for e in entries for i in e['images']),
    'body_images_to_upload':sum(not i['url'] for e in entries for i in e['images']),
    'covers_reused':sum(bool(e['cover_media_id']) for e in entries)},ensure_ascii=False))
