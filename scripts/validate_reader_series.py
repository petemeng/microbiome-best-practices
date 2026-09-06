#!/usr/bin/env python3
"""Execute each chapter's public R script in a fresh directory, without local data.

Use --article-dir to additionally check the packaged WeChat code. This is a
reader-reproduction check, not a replacement for scientific-method validators.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

from reader_reproducibility import frontmatter, reader_blocks, script_text, validate_reader_contract


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_one(root, number, work_root, timeout, article_dir):
    source = next((root / 'chapters').glob(f'{number:02}-*.qmd'))
    text = source.read_text()
    meta = frontmatter(text)
    config = meta.get('reader-reproduction', {})
    result = {'chapter': number, 'status': 'failed', 'checks': []}
    def check(name, passed, **detail):
        result['checks'].append({'name': name, 'passed': bool(passed), **detail})
    try:
        blocks = reader_blocks(text)
        script = root / config['script']
        starting_sha = sha256(script)
        check('script_matches_public_source', script.read_text() == script_text(text, blocks))
        if article_dir:
            candidates = list(article_dir.glob(f'**/{number:02}*/article.html'))
            if len(candidates) != 1:
                raise ValueError('Exactly one packaged article is required')
            validate_reader_contract(source, candidates[0].read_text(), root)
            check('wechat_code_complete_and_ordered', True)
        work = work_root / f'{number:02}'
        work.mkdir(parents=True, exist_ok=False)
        check('started_with_empty_working_directory', not any(work.iterdir()))
        snapshots = work_root / 'scripts'
        snapshots.mkdir(exist_ok=True)
        snapshot = snapshots / f'{number:02}.R'
        snapshot.write_text(script.read_text())
        env = os.environ.copy()
        env['R_LIBS_USER'] = str(root / '.r-lib')
        for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
            env[key] = '1'
        started = time.monotonic()
        log_path = work_root / f'{number:02}.log'
        with log_path.open('w') as log:
            completed = subprocess.run(['Rscript', '--vanilla', str(snapshot)], cwd=work,
                env=env, stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
        result['elapsed_seconds'] = round(time.monotonic()-started, 1)
        check('fresh_R_process_completed', completed.returncode == 0, exit_code=completed.returncode)
        if completed.returncode:
            # Keep full diagnostics in the local log, not in public material.
            result['error'] = log_path.read_text(errors='replace')[-2500:].replace(str(root), '[project]').replace(str(work_root), '[validation]')
        downloaded = sorted((work / 'data').rglob('*')) if (work/'data').exists() else []
        for path in downloaded:
            if path.is_file():
                relative = path.relative_to(work)
                original = root / relative
                check('downloaded_input_matches_reference', original.is_file() and sha256(path)==sha256(original), file=str(relative))
        check('real_inputs_downloaded', any(p.is_file() for p in downloaded))
        expected = set(re.findall(r'\]\(\.\./(figures/[^)]+\.png)\)', text))
        for name in re.findall(r'"(figures/[^"\n]+)"', snapshot.read_text()):
            if name.endswith('/') or Path(name).suffix not in ('', '.png', '.pdf', '.svg', '.tiff'):
                continue
            expected.add(str(Path(name).with_suffix('.png')) if Path(name).suffix else name+'.png')
        expected = sorted(expected)
        for png in expected:
            path = work / png
            check('figure_recreated', path.is_file() and path.stat().st_size>1000, file=png)
            check('vector_export_recreated', path.with_suffix('.pdf').is_file() or path.with_suffix('.svg').is_file(), file=png)
        check('figures_exported', any((work/'figures').rglob('*.png')) if (work/'figures').exists() else False)
        result['script_sha256'] = starting_sha
        check('source_unchanged_during_run', source.read_text() == text and sha256(script) == starting_sha)
        result['executed_code_blocks'] = len(blocks)
        result['scope'] = 'fresh public R execution, input identity, figure exports; method correctness has separate tests'
    except Exception as exc:
        result['error'] = str(exc).replace(str(root), '[project]').replace(str(work_root), '[validation]')
        check('validation_completed', False)
    result['status'] = 'passed' if result['checks'] and all(c['passed'] for c in result['checks']) else 'failed'
    (work_root / f'{number:02}.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--chapters', type=int, nargs='+', required=True)
    p.add_argument('--work-root', type=Path, required=True)
    p.add_argument('--article-dir', type=Path)
    p.add_argument('--jobs', type=int, default=2)
    p.add_argument('--timeout', type=int, default=3600)
    args = p.parse_args()
    args.work_root.mkdir(parents=True, exist_ok=True)
    results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        jobs=[pool.submit(run_one,args.project_root.resolve(),n,args.work_root.resolve(),args.timeout,args.article_dir) for n in args.chapters]
        for job in concurrent.futures.as_completed(jobs):
            result=job.result();results.append(result)
            print(json.dumps({'chapter':result['chapter'],'status':result['status'],'error':result.get('error')},ensure_ascii=False),flush=True)
            (args.work_root/'report.json').write_text(json.dumps({'status':'passed' if len(results)==len(jobs) and all(r['status']=='passed' for r in results) else 'incomplete_or_failed','results':sorted(results,key=lambda r:r['chapter'])},ensure_ascii=False,indent=2)+'\n')
    raise SystemExit(0 if all(r['status']=='passed' for r in results) else 1)


if __name__=='__main__':
    main()
