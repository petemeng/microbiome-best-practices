#!/usr/bin/env python3
"""Merge verified draft replacements into the private canonical bundle and map.

This command does not contact WeChat. Each selected journal must record creation,
full-body readback and deletion of its exact previously mapped draft. It is safe
to resume after an interrupted local merge.
"""
import argparse,hashlib,json,os,shutil
from datetime import datetime,timezone
from pathlib import Path

def atomic_json(path,value):
    temporary=path.with_name(path.name+'.tmp')
    fd=os.open(temporary,os.O_WRONLY|os.O_CREAT|os.O_TRUNC,0o600)
    with os.fdopen(fd,'w') as stream:
        json.dump(value,stream,ensure_ascii=False,indent=2)
        stream.write('\n');stream.flush();os.fsync(stream.fileno())
    os.replace(temporary,path)

def repath(value,old,new):
    if isinstance(value,str):return value.replace(str(old),str(new))
    if isinstance(value,list):return [repath(x,old,new) for x in value]
    if isinstance(value,dict):return {k:repath(v,old,new) for k,v in value.items()}
    return value

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--journal-directory',type=Path,required=True)
    p.add_argument('--journal-prefix',required=True)
    p.add_argument('--destination',type=Path,required=True)
    args=p.parse_args()
    plan=json.loads(args.plan.read_text());map_path=Path(plan['live_map'])
    live=json.loads(map_path.read_text());mapped={e['chapter_id']:e for e in live['entries']}
    journals={}
    for e in plan['entries']:
        key=e['chapter_id']
        j=json.loads((args.journal_directory/f'{args.journal_prefix}-{key}.json').read_text())
        if not all(j.get(k) for k in ('created','verified','full_body_verified','old_draft_deleted')):
            raise ValueError(f'Incomplete remote replacement: {key}')
        if j['old_draft_media_id']!=e['old_draft_media_id'] or j['source_commit']!=e['source_commit']:
            raise ValueError(f'Journal does not match reviewed plan: {key}')
        if mapped[key]['draft_media_id'] not in (e['old_draft_media_id'],j['draft_media_id']):
            raise ValueError(f'Canonical mapping changed independently: {key}')
        journals[key]=j
    old_bundle=Path(live['source_bundle']).resolve();destination=args.destination.resolve()
    if not destination.exists():
        shutil.copytree(old_bundle,destination)
        os.chmod(destination,0o700)
    backup=args.journal_directory/f'{args.journal_prefix}-map-before.json'
    if not backup.exists():atomic_json(backup,live)
    bundle=Path(plan['bundle']);new_report=json.loads((bundle/'report.json').read_text())
    report=repath(json.loads((destination/'report.json').read_text()),old_bundle,destination)
    items={i['chapter_id']:i for i in report['items']}
    revised={i['chapter_id']:i for i in new_report['items']}
    for e in plan['entries']:
        key=e['chapter_id'];j=journals[key]
        shutil.copytree(bundle/key,destination/key,dirs_exist_ok=True)
        body=e['draft']['content']
        for image in j['images']:
            file=Path(image['file'])
            if hashlib.sha256(file.read_bytes()).hexdigest()!=image['sha256']:
                raise ValueError(f'Reviewed asset changed: {key}')
            body=body.replace('src="'+image['src']+'"','src="'+image['url']+'"')
        if 'src="images/' in body:raise ValueError(f'Unresolved body image: {key}')
        draft=dict(e['draft'],content=body,thumb_media_id=j['cover_media_id'])
        atomic_json(destination/key/'draft.json',draft)
        mapped[key].update(title=e['title'],author='Peter',old_draft_media_id=e['old_draft_media_id'],
            draft_media_id=j['draft_media_id'],cover_media_id=j['cover_media_id'],
            cover_kind='representative_figure_only',
            body_images=[{k:i[k] for k in ('src','url','resolution')} for i in j['images']],
            body_image_count=len(j['images']),created=True,verified=True,old_draft_deleted=True,
            full_body_verified=bool(j.get('full_body_verified')),
            source_commit=e['source_commit'],content_sha256=hashlib.sha256(body.encode()).hexdigest(),
            verification_scope=j['verification_scope'])
        items[key]=repath(revised[key],bundle,destination)
    now=datetime.now(timezone.utc).isoformat()
    report.update(generated_at=now,items=[items[k] for k in sorted(items)],
        item_count=len(items),selected_chapters=[int(k) for k in sorted(items)],status='passed')
    report.setdefault('revision_batches',[]).append({'chapters':sorted(journals),'source_commit':plan['source_commit'],'completed_at':now})
    atomic_json(destination/'report.json',report)
    live.update(generated_at=now,source_bundle=str(destination),entries=[mapped[k] for k in sorted(mapped)])
    live['github'].update(commit=plan['source_commit'],pushed=True)
    live['wechat'].update(updated_chapters=sorted(journals),drafts_created=len(journals),
        old_drafts_deleted=len(journals),existing_series_drafts_retained=55-len(journals),
        inventory_verified=False,published=False,mass_sent=False)
    atomic_json(map_path,live)
    print(json.dumps({'updated_chapters':sorted(journals),'canonical_entries':len(mapped),
        'full_remote_body_readback':all(j.get('full_body_verified') for j in journals.values()),
        'inventory_recheck_required':True}))

if __name__=='__main__':main()
