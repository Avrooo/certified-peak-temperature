import json,glob,time,signal,numpy as np,e4_floor as E, e3
R=[json.loads(l) for f in sorted(glob.glob('e4_G*.jsonl')) for l in open(f)]
class TO(Exception): pass
def _h(s,f): raise TO()
signal.signal(signal.SIGALRM,_h)
import sys
floors=[1e-3,1e-5,1e-7,1e-9,1e-11,1e-13]; out={}
for fl in floors:
    E.FLOOR=fl; rows=[]; t0=time.time()
    for r in R:
        box=tuple(map(tuple,r['box'])); E.clr(); signal.alarm(30)
        try: d,w=E.b6_decision_first(box,r['Fs'],r['Fp'],r['L'])
        except TO: d,w='TIMEOUT',None
        except Exception as ex: d,w='ERROR:'+type(ex).__name__,None
        finally: signal.alarm(0)
        rows.append((d,w,r['truth'],r['sep']))
    out[str(fl)]=rows
    dec=sum(1 for d,w,t,s in rows if d in('SAFE_ALL','EXCEEDS_ALL','MIXED')); dis=sum(1 for d,w,t,s in rows if d in('SAFE_ALL','EXCEEDS_ALL','MIXED') and d!=t)
    print(fl,'decided',dec,'disagree',dis,'inconclusive',sum(d=='INCONCLUSIVE' for d,*_ in rows),'other',sum(d.startswith(('TIME','ERR')) for d,*_ in rows),f'{time.time()-t0:.0f}s',flush=True)
    json.dump(out,open('floorrun.json','w'))
