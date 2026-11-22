"""Theorem-dependence ablation on the 1,600-case benchmark. All variants use the same float spectral evaluator."""
import json, sys, time, numpy as np, e4
rng=np.random.default_rng(20261122)
cases=[json.loads(l) for g in ('G1','G2','G3','G4') for l in open(f'e4_{g}.jsonl')]
N=int(sys.argv[1]) if len(sys.argv)>1 else len(cases); cases=cases[:N]
def dec_corners(pU,pL,L): return 'SAFE_ALL' if pU<L else 'EXCEEDS_ALL' if pL>=L else 'MIXED'
def dec_samples(peaks,L):
    a=np.array(peaks)
    return 'SAFE_ALL' if np.all(a<L) else 'EXCEEDS_ALL' if np.all(a>=L) else 'MIXED'
def endpeak(B1,B2,Fs,Fp): return max(e4.endpoints_float(B1,B2,Fs,Fp))
V=[('end time only, corners',1),('uniform grid 5 pts, corners',5),('uniform grid 20 pts, corners',20),
   ('uniform grid 50 pts, corners',50),('uniform grid 200 pts, corners',200),('Theorem 3.5 endpoints, corners (corner rule)','T'),
   ('Theorem 3.5 endpoints, 20 random box samples','S20'),('Theorem 3.5 endpoints, 100 random box samples','S100')]
res={v[0]:[] for v in V}; tim={v[0]:[] for v in V}
for c in cases:
    (b1l,b1h),(b2l,b2h)=c['box']; Fs,Fp,L=c['Fs'],c['Fp'],c['L']; tr=c['truth']
    for name,spec in V:
        t0=time.perf_counter()
        if spec=='T': d=dec_corners(endpeak(b1h,b2h,Fs,Fp),endpeak(b1l,b2l,Fs,Fp),L)
        elif isinstance(spec,int):
            d=dec_corners(e4.full_history_peak(b1h,b2h,Fs,Fp,nt=spec),e4.full_history_peak(b1l,b2l,Fs,Fp,nt=spec),L)
        else:
            k=int(spec[1:]); pk=[endpeak(rng.uniform(b1l,b1h),rng.uniform(b2l,b2h),Fs,Fp) for _ in range(k)]; d=dec_samples(pk,L)
        tim[name].append(time.perf_counter()-t0); res[name].append((d,tr,c['sep']))
out={}
for name,_ in V:
    r=res[name]; dis=sum(d!=t for d,t,s in r); unsafe=sum(d=='SAFE_ALL' and t!='SAFE_ALL' for d,t,s in r)
    far=[x for x in r if x[2]>=1e-14]; disf=sum(d!=t for d,t,s in far)
    out[name]=dict(n=len(r),disagree=dis,unsafe_called_safe=unsafe,n_gap_ge_1e14=len(far),disagree_gap_ge_1e14=disf,median_ms=1e3*float(np.median(tim[name])))
    print(f"{name:50s} dis={dis:5d} unsafe={unsafe:5d} | gap>=1e-14: {disf}/{len(far)} | {out[name]['median_ms']:.2f} ms")
json.dump(out,open('ablate.json','w'),indent=1)
