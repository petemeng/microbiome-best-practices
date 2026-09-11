#!/usr/bin/env python3
"""Recheck existing fresh-run artifacts without rerunning or relabelling a computation.

Original run reports are kept untouched. Promotion requires a successful R exit,
the exact executed snapshot, unchanged reader code, and every non-export check.
"""
import argparse,hashlib,json,re
from pathlib import Path
from reader_reproducibility import frontmatter,reader_blocks,script_text

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--work-root',type=Path,required=True)
p.add_argument('--project-root',type=Path,default=Path(__file__).resolve().parents[1])
p.add_argument('--output',type=Path,required=True)
args=p.parse_args()
export_names={'figure_recreated','vector_export_recreated','figures_exported'}
results=[]
for report in sorted(args.work_root.glob('[0-9]*.json')):
    original=json.loads(report.read_text());n=original['chapter']
    source=next((args.project_root/'chapters').glob(f'{n:02}-*.qmd'))
    text=source.read_text();expected_script=script_text(text,reader_blocks(text))
    snapshot=args.work_root/'scripts'/f'{n:02}.R'
    checks=[dict(c) for c in original['checks'] if c['name'] not in export_names]
    def check(name,ok,**fields):checks.append({'name':name,'passed':bool(ok),**fields})
    check('executed_snapshot_matches_current_code',snapshot.is_file() and snapshot.read_text()==expected_script)
    check('snapshot_matches_execution_record',snapshot.is_file() and hashlib.sha256(snapshot.read_bytes()).hexdigest()==original.get('script_sha256'))
    check('successful_R_exit_recorded',any(c['name']=='fresh_R_process_completed' and c['passed'] for c in checks))
    work=args.work_root/f'{n:02}'
    expected=set(re.findall(r'\]\(\.\./(figures/[^)]+\.png)\)',text))
    for name in re.findall(r'"(figures/[^"\n]+)"',expected_script):
        path=Path(name)
        if name.endswith('/') or path.suffix not in ('','.png','.pdf','.svg','.tiff'):continue
        expected.add(str(path.with_suffix('.png')) if path.suffix else name+'.png')
    for name in sorted(expected):
        path=work/name
        check('figure_recreated',path.is_file() and path.stat().st_size>1000,file=name)
        check('vector_export_recreated',path.with_suffix('.pdf').is_file() or path.with_suffix('.svg').is_file(),file=name)
    check('figures_exported',any((work/'figures').rglob('*.png')) if (work/'figures').exists() else False)
    result={'chapter':n,'status':'passed' if all(c['passed'] for c in checks) else 'failed','checks':checks,
        'original_run_report':str(report),'original_status':original['status'],'script_sha256':original.get('script_sha256'),
        'scope':'Rechecked nested figure exports from the same recorded fresh R run; no new execution is claimed.'}
    results.append(result)
output={'status':'passed' if results and all(r['status']=='passed' for r in results) else 'failed','results':results}
args.output.write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':output['status'],'chapters':len(results),'failed':[r['chapter'] for r in results if r['status']=='failed']}))
raise SystemExit(0 if output['status']=='passed' else 1)
