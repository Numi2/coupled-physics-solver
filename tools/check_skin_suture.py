"""Independent float64 observations of exported native states, never physics."""
import csv, json, math
from pathlib import Path
from wound_geometry import check_wound

def rows(path):
    with Path(path).open() as f: return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
def add(a,b): return [x+y for x,y in zip(a,b)]
def sub(a,b): return [x-y for x,y in zip(a,b)]
def scale(a,s): return [x*s for x in a]
def norm(a): return math.sqrt(sum(x*x for x in a))
def check(directory, closure=False):
    d=Path(directory); samples=rows(d/'traction-samples.csv'); lips=rows(d/'lip-geometry.csv')
    assert samples and (d/'exit.txt').read_text().strip()=='0', 'native execution incomplete'
    command=json.loads((d/'command.json').read_text())
    planned=int(command[command.index('--steps')+1]) if '--steps' in command else 60
    assert len(samples)==planned, 'missing planned accepted frames'
    assert len(samples)==round(samples[-1]['time_s']/samples[0]['time_s']), 'missing accepted frames'
    max_gap_error=0
    for s in samples:
        points=[p for p in lips if p['step']==s['step']]
        assert len(points)==12
        gaps=[norm([p[f'upper_{a}_m']-p[f'lower_{a}_m'] for a in 'xyz']) for p in points]
        max_gap_error=max(max_gap_error,abs(sum(gaps)/len(gaps)-s['mean_gap_m']))
        assert all(p['upper_y_m']>=p['lower_y_m'] for p in points), 'crossed lips'
        assert s['minimum_J']>.2 and s['active_tetrahedra']==2304
        assert s['equilibrium_residual']<=1e-4 and s['volume_residual']<=1e-4 and s['pressure_residual']<=1e-4
        assert abs(s['mass_kg']/samples[0]['mass_kg']-1)<=1e-5
    assert max_gap_error<1e-10
    thread=rows(d/'thread.csv'); report={'frames':len(samples),'initial_gap_mm':samples[0]['mean_gap_m']*1000,'final_gap_mm':samples[-1]['mean_gap_m']*1000,'closure_fraction':1-samples[-1]['mean_gap_m']/samples[0]['mean_gap_m'],'gap_reconstruction_error_m':max_gap_error,'minimum_J':min(s['minimum_J'] for s in samples)}
    if thread:
        source=rows(d/'thread-source.csv'); anchors=rows(d/'thread-anchors.csv')
        geometry=rows(d/'tissue-geometry.csv'); forces=rows(d/'thread-forces.csv')
        nodes={}; actual={}
        for p in geometry: nodes.setdefault(int(p['step']),{})[int(p['node'])]=[p[a+'_m'] for a in 'xyz']
        for p in forces: actual.setdefault(int(p['step']),{})[int(p['node'])]=[p['f'+a+'_n'] for a in 'xyz']
        faces=[tuple(int(r[a]) for a in 'abc') for r in rows(d/'surface.csv')]
        report.update(check_wound(nodes,faces))
        force_error=0; tension_error=0; material_error=0; gradient_error=0
        off=all(p['tension_n']==0 for p in thread)
        for step in sorted(nodes):
            predicted={}
            for material in source:
                loop=int(material['stitch']); aa=[a for a in anchors if a['stitch']==loop]
                corners=[]
                for c in range(4):
                    selected=[a for a in aa if a['corner']==c];assert abs(sum(a['weight'] for a in selected)-1)<1e-6
                    v=[0.,0.,0.]
                    for a in selected:v=add(v,scale(nodes[step][int(a['node'])],a['weight']))
                    corners.append(v)
                edge=[sub(corners[(c+1)%4],corners[c]) for c in range(4)];L=sum(norm(e) for e in edge)
                t=[scale(e,1/norm(e)) for e in edge]
                record=next(p for p in thread if p['step']==step and p['stitch']==loop)
                rest=record['loop_rest_m'];T=0 if off else material['EA_n']*max(L/rest-1,0)
                tension_error=max(tension_error,abs(T-record['tension_n']))
                material_error=max(material_error,abs(rest+record['tail_rest_m']-material['total_rest_m']))
                for a in aa:
                    c=int(a['corner']);n=int(a['node']); f=scale(sub(t[c],t[(c+3)%4]),a['weight']*T)
                    predicted[n]=add(predicted.get(n,[0.,0.,0.]),f)
                # Independent central difference of elastic energy at one corner.
                if T>1e-5:
                    delta=1e-8;energies=[]
                    for sign in [-1,1]:
                        cc=[p[:] for p in corners];cc[0][1]+=sign*delta
                        ll=sum(norm(sub(cc[(c+1)%4],cc[c])) for c in range(4))
                        energies.append(.5*material['EA_n']*max(ll-rest,0)**2/rest)
                    derivative=(energies[1]-energies[0])/(2*delta)
                    gradient_error=max(gradient_error,abs(derivative+T*(t[0][1]-t[3][1])))
            for n in set(predicted)|set(actual.get(step,{})):
                force_error=max(force_error,norm(sub(predicted.get(n,[0.,0.,0.]),actual.get(step,{}).get(n,[0.,0.,0.]))))
            balance=[0.,0.,0.]
            for f in actual.get(step,{}).values(): balance=add(balance,f)
            assert norm(balance)<1e-6,'thread reaction does not balance'
        assert force_error<1e-6 and tension_error<1e-6 and material_error<1e-8 and gradient_error<1e-8
        markers=rows(d/'thread-geometry.csv');first={int(p['stitch']):p for p in thread if p['step']==0};last={int(p['stitch']):p for p in thread if p['step']==samples[-1]['step']}
        exiting=sum(first[int(p['stitch'])]['loop_rest_m']>=p['material_m']>last[int(p['stitch'])]['loop_rest_m'] for p in markers if p['step']==samples[-1]['step'])
        report.update(thread_force_error_n=force_error,thread_tension_error_n=tension_error,energy_gradient_error_n=gradient_error,material_length_error_m=material_error,markers_pulled_through=exiting,pull_mm=max(p['pulled_m'] for p in thread)*1000,peak_tension_mn=max(p['tension_n'] for p in thread)*1000)
        if closure:
            assert not off and exiting>=3 and report['pull_mm']>=1.19,'no actual material passage'
            assert report['closure_fraction']>=.7,'wound did not close sufficiently'
            assert report['final_gap_mm']<=.18,'remaining wound gap exceeds bench target'
        elif off:
            assert max(abs(s['mean_gap_m']-samples[0]['mean_gap_m']) for s in samples)<1e-6,'no-thread control moved'
    (d/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':
    import sys
    print(json.dumps(check(sys.argv[1], '--closure' in sys.argv),indent=2))
