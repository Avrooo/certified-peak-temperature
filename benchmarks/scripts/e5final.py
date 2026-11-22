import json,glob,sys,time,signal,math,numpy as np
from fractions import Fraction as F
import stacks
FR=stacks.load_stack("frozen","frozen",256,50,52,0,0,0)
V7=stacks.load_stack("v7",".",96,24,24,1,1,1)
FB=stacks.load_stack("fb",".",256,50,52,1,1,1)       # fallback kernel: 256 bit, validated optimisations on
R=[json.loads(l) for f in sorted(glob.glob('e4_G*.jsonl')) for l in open(f)]
class TO(Exception): pass
def _h(s,f): raise TO()
signal.signal(signal.SIGALRM,_h)
def run(st,r,limit=20):
    box=tuple(map(tuple,r['box'])); stacks.clear(st); 
    if hasattr(st.e3,'FAILS'): st.e3.FAILS[0]=0
    signal.alarm(limit); t0=time.perf_counter()
    try: d,w=st.e4.b6_decision_first(box,r['Fs'],r['Fp'],r['L']); err=None
    except TO: d,w,err='TIMEOUT',None,'timeout'
    except Exception as ex: d,w,err='ERROR',None,type(ex).__name__
    finally: signal.alarm(0)
    return dict(dec=d,w=w,t=time.perf_counter()-t0,err=err,fails=getattr(st.e3,'FAILS',[0])[0],P=st.ivk.P)
P_ADEQUATE=72     # 2^-72*1e6 ~ 2e-16 < tolerance floor 1e-13: below this the rounding floor can dominate a requested tolerance
def with_fallback(r,v7res):
    invalid=v7res.get('err')=='ValueError'          # invalid input: precision cannot help, report it
    if not invalid and (v7res['dec'] in ('ERROR','TIMEOUT') or (v7res['dec']=='INCONCLUSIVE' and (v7res['fails']>0 or v7res.get('P',256)<P_ADEQUATE))):
        fb=run(FB,r); return dict(dec=fb['dec'],t=v7res['t']+fb['t'],fallback=True,w=fb['w'])
    return dict(dec=v7res['dec'],t=v7res['t'],fallback=False,w=v7res['w'])
if __name__=="__main__":
    out=[]
    for i,r in enumerate(R):
        fr=run(FR,r); v7=run(V7,r); fb=with_fallback(r,v7)
        out.append(dict(seed=r['seed'],group=r['group'],kind=r['kind'],cls=r['cls'],truth=r['truth'],sep=r['sep'],fr=fr,v7=v7,fb=fb))
        if i%200==0: print(i,flush=True)
    json.dump(out,open('e5final_main.json','w'),default=float); print("done")
