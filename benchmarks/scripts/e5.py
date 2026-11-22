"""E5 timing: run from a directory that holds the e3/e4/ivk variant. Fixed stratified sample, single process,
cold start per case (all caches cleared), decisions compared to the stored truth and to the E4 A2 decision."""
import json,glob,random,sys,time,os,numpy as np
sys.path.insert(0,os.getcwd())
import e4,e3
src=sys.argv[1]; reps=int(sys.argv[2]); label=sys.argv[3]
R=[json.loads(l) for f in sorted(glob.glob(src+'/e4_G*.jsonl')) for l in open(f)]
random.seed(7)
S=[]
for g in ('G1','G2','G3','G4'):
    pool=[r for r in R if r['group']==g and r['B6']['dec']!='INCONCLUSIVE']; random.shuffle(pool); S+=pool[:100]
def clr():
    e3._rootcache.clear(); e3._modecache.clear()
    if hasattr(e3,'_sccache'): e3._sccache.clear()
allt=[]; decs=None
for rep in range(reps):
    ts=[]; ds=[]
    for r in S:
        box=tuple(map(tuple,r['box'])); clr(); t0=time.perf_counter()
        d,w=e4.b6_decision_first(box,r['Fs'],r['Fp'],r['L']); ts.append(time.perf_counter()-t0); ds.append((d,w))
    allt.append(ts); decs=ds
T=np.min(np.array(allt),axis=0)        # per-case best of reps (removes scheduler noise)
wrong=sum(1 for r,(d,w) in zip(S,decs) if d!=r['truth']); incon=sum(1 for d,w in decs if d=='INCONCLUSIVE')
ws=[w for d,w in decs if w]
print(f"{label:34s} n={len(S)} median {np.median(T)*1e3:6.2f} ms  mean {np.mean(T)*1e3:6.2f}  P90 {np.percentile(T,90)*1e3:6.2f}  max {T.max()*1e3:6.1f} | wrong {wrong} inconclusive {incon} | median width {np.median(ws):.1e}")
json.dump(dict(label=label,t=T.tolist(),dec=[d for d,w in decs],w=[w for d,w in decs]),open(f"e5_{label.replace(' ','_').replace('+','p')}.json","w"))
