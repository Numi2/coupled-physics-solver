import csv,json
from pathlib import Path

def prepare(directory):
    d=Path(directory)
    def read(name):
        with (d/name).open() as f:return list(csv.DictReader(f))
    tissue={};thread={}
    for row in read('tissue-geometry.csv'):
        tissue.setdefault(int(row['step']),{})[int(row['node'])]=[float(row[a+'_m'])*1000 for a in 'xyz']
    for row in read('thread-geometry.csv'):
        thread.setdefault(int(row['step']),{}).setdefault(int(row['stitch']),{})[int(row['marker'])]=[float(row[a+'_m'])*1000 for a in 'xyz']
    stats=read('thread.csv');lips=read('lip-geometry.csv');frames=[]
    for s in read('traction-samples.csv'):
        step=int(s['step']);tt=[t for t in stats if int(t['step'])==step]
        ll=[p for p in lips if int(p['step'])==step]
        wound=[[[float(p[side+'_'+a+'_m'])*1000 for a in 'xyz'] for p in ll] for side in ['lower','upper']]
        frames.append(dict(lips=wound,tissue=[v for k,v in sorted(tissue[step].items())],thread=[[v for k,v in sorted(points.items())] for _,points in sorted(thread.get(step,{}).items())],gap=float(s['mean_gap_m'])*1000,pull=max([float(t['pulled_m'])*1000 for t in tt] or [0]),tension=max([float(t['tension_n'])*1000 for t in tt] or [0]),time=float(s['time_s'])))
    triangles=[int(row[a]) for row in read('surface.csv') for a in 'abc']
    assert triangles and all(0<=v<len(frames[0]['tissue']) for v in triangles)
    out=d/'replay.json';out.write_text(json.dumps(dict(triangles=triangles,frames=frames),separators=(',',':')));return out
if __name__=='__main__':
    import sys
    print(prepare(sys.argv[1]))
