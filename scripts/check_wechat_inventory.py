#!/usr/bin/env python3
"""Read-only inventory validation for the mapped 55-article series."""
import argparse,json,re
from pathlib import Path
from sync_wechat_revision import Client,SafeError
from finalize_wechat_revision import atomic_json

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--token-file',type=Path,required=True)
    p.add_argument('--live-map',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--baseline',type=Path)
    p.add_argument('--record-live-map',action='store_true')
    args=p.parse_args();client=Client(args.token_file)
    total=client.api('draft/count')['total_count'];items=[]
    for offset in range(0,total,20):
        items.extend(client.api('draft/batchget',{'offset':offset,'count':20,'no_content':1}).get('item',[]))
    live=json.loads(args.live_map.read_text());mapped={e['chapter_id']:e for e in live['entries']}
    series=[];unrelated=[]
    for item in items:
        articles=item['content']['news_item']
        matches=[(a,re.match(r'^16S最佳实践｜([1-9][0-9]?)\. ',a['title'])) for a in articles]
        matches=[(a,m) for a,m in matches if m]
        if matches:
            if len(articles)!=1 or len(matches)!=1:raise SafeError('Series draft contains multiple articles')
            a,m=matches[0]
            series.append({'chapter_id':f'{int(m.group(1)):02}','media_id':item['media_id'],'title':a['title'],'author':a.get('author')})
        else:unrelated.append(item['media_id'])
    ids={i['media_id'] for i in items};expected={e['draft_media_id'] for e in mapped.values()}
    old={e.get('old_draft_media_id') for e in mapped.values()}-expected
    checks={'total':total,'listed':len(items),'unique_drafts':len(ids),'series':len(series),
        'unique_ordinals':len({s['chapter_id'] for s in series}),
        'missing_mapped_chapters':[k for k,e in mapped.items() if e['draft_media_id'] not in ids],
        'all_mapped_old_drafts_absent':not bool(old & ids),'unrelated_drafts':len(unrelated),
        'mapped_identity_matches':all(s['chapter_id'] in mapped and s['media_id']==mapped[s['chapter_id']]['draft_media_id']
            and s['title']==mapped[s['chapter_id']]['title'] and s['author']=='Peter' for s in series)}
    if args.baseline:
        baseline=json.loads(args.baseline.read_text())
        checks['unrelated_drafts_retained']=set(baseline['unrelated_ids'])==set(unrelated)
    passed=(total==len(items)==len(ids) and len(series)==55 and {s['chapter_id'] for s in series}==set(mapped)
        and not checks['missing_mapped_chapters'] and checks['all_mapped_old_drafts_absent']
        and checks['mapped_identity_matches'] and checks.get('unrelated_drafts_retained',True))
    report={'status':'passed' if passed else 'failed','checks':checks,'unrelated_ids':sorted(unrelated)}
    atomic_json(args.output,report)
    if args.record_live_map and passed:
        live['wechat'].update(final_draft_count=total,series_draft_count=55,unrelated_draft_count=len(unrelated),
            inventory_verified=True,inventory_checks=checks,published=False,mass_sent=False)
        atomic_json(args.live_map,live)
    print(json.dumps({'status':report['status'],'checks':checks},ensure_ascii=False))
    return not passed

if __name__=='__main__':
    try:raise SystemExit(main())
    except SafeError as error:print(str(error));raise SystemExit(1)
