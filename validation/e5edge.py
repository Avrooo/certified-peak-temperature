import json,glob,signal,time,numpy as np
import stacks
import e5final as E
V7=E.V7; FB=E.FB; FR=E.FR
# --- (A) zero modes inside the benchmark ---
C=E.R; zero=[r for r in C if r['box'][1][0]==0.0 or r['box'][0][0]==0.0]
D={d['seed']:d for d in json.load(open('e5final_main.json'))}
print("(A) zero-mode cases (B=0 at a corner):",len(zero),"| wrong V7:",sum(D[r['seed']]['v7']['dec'] not in('INCONCLUSIVE',r['truth']) for r in zero),"| V7 differs from frozen:",sum(D[r['seed']]['v7']['dec']!=D[r['seed']]['fr']['dec'] for r in zero))
# --- (B) arithmetic failure / degenerate inputs: outcome must be a decision consistent with truth, INCONCLUSIVE, or a clean failure -- never a wrong claim ---
base=dict(Fs=0.05,Fp=0.2)
edge=[("tau=0 (Fs==Fp)",((3,3.2),(1,1.1)),0.05,0.05,0.4),("Fs=0 invalid",((3,3.2),(1,1.1)),0.0,0.2,0.4),("Fp<Fs invalid",((3,3.2),(1,1.1)),0.2,0.05,0.4),
 ("huge B 1e9",((1e9,2e9),(1e8,2e8)),0.05,0.2,0.9),("tiny B 1e-12",((1e-12,2e-12),(1e-12,2e-12)),0.05,0.2,0.01),("B1=B2=0",((0.0,0.0),(0.0,0.0)),0.05,0.2,0.5),
 ("B1=0,B2>0",((0.0,0.0),(2,3)),0.05,0.2,0.5),("B2=0,B1>0",((2,3),(0.0,0.0)),0.05,0.2,0.05),("L=0",((3,3.2),(1,1.1)),0.05,0.2,0.0),("L=1",((3,3.2),(1,1.1)),0.05,0.2,1.0),
 ("long pulse Fp=1e3",((3,3.2),(1,1.1)),0.05,1000.0,0.9999),("point box B1=B2",((3,3),(3,3)),0.05,0.2,0.5)]
print("(B) arithmetic-failure / degenerate inputs (frozen | V7 | V7+fallback):")
for name,box,Fs,Fp,L in edge:
    out=[]
    for st in (FR,V7):
        stacks.clear(st)
        if hasattr(st.e3,'FAILS'): st.e3.FAILS[0]=0
        signal.alarm(30)
        try: d,w=st.e4.b6_decision_first(box,Fs,Fp,L); out.append(d)
        except E.TO: out.append("TIMEOUT")
        except Exception as ex: out.append("ERR:"+type(ex).__name__)
        finally: signal.alarm(0)
    rv={'dec':out[1],'fails':V7.e3.FAILS[0]}; fb=E.with_fallback(dict(box=box,Fs=Fs,Fp=Fp,L=L),dict(dec=out[1],fails=V7.e3.FAILS[0],t=0,w=None))
    print(f"   {name:20s} frozen={out[0]:13s} V7={out[1]:13s} V7fails={V7.e3.FAILS[0]} fallback_used={fb['fallback']} -> {fb['dec']}")
# --- (C) precision exhaustion: crippled kernels must degrade to INCONCLUSIVE/failure, never a wrong decision ---
sample=[]; import random; random.seed(7)
for g in ('G1','G2','G3','G4'):
    pool=[r for r in C if r['group']==g and r['sep']>=1e-14]; random.shuffle(pool); sample+=pool[:60]
print("(C) precision exhaustion on",len(sample),"cases (separation >=1e-14):")
for P,K,N in ((40,10,10),(56,14,14),(64,16,16),(80,20,20),(96,24,24)):
    st=stacks.load_stack(f"p{P}",".",P,K,N,1,1,1); c=w=i=er=fb_ok=0; wrong_fb=0
    for r in sample:
        res=E.run(st,r,limit=20)
        d=res['dec']
        if d in('ERROR','TIMEOUT'): er+=1
        elif d=='INCONCLUSIVE': i+=1
        elif d==r['truth']: c+=1
        else: w+=1
        fbres=E.with_fallback(r,res) if (res['dec'] in('ERROR','TIMEOUT') or (res['dec']=='INCONCLUSIVE')) else None
        if fbres and (res['dec'] in('ERROR','TIMEOUT') or res['fails']>0): 
            fb_ok+= fbres['dec']==r['truth']
    print(f"   P={P:3d} bits (K=N={K}): correct {c:3d} wrong {w} inconclusive {i:3d} errors {er} | fallback(256) after failure -> correct {fb_ok}")
