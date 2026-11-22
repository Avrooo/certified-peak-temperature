"""Independent finite-volume check, 16 cases. Cell-centred FV, backward Euler, dt=0.25*h^2*? (dt ~ h^2), exact switch step."""
import json,glob,math,numpy as np
from scipy.linalg import solve_banded
import e3
from fractions import Fraction
R=[json.loads(l) for f in sorted(glob.glob('e4_G*.jsonl')) for l in open(f)]
def fv(B1,B2,Fs,Fp,N,c=1.0):
    h=1.0/N; dt0=c*h*h
    n1=max(1,math.ceil(Fs/dt0)); n2=max(1,math.ceil((Fp-Fs)/dt0))
    u=np.ones(N); th=[]; 
    def ab(B,dt):
        r=dt/h**2; a=np.zeros((3,N)); a[1,:]=1+2*r; a[0,1:]=-r; a[2,:-1]=-r
        a[1,0]=1+r                      # Neumann at x=0
        a[1,N-1]=1+r+r*(2*B*h/(2+B*h))*1.0/(1.0)*0+ r*( B*h/(1+B*h/2))  # Robin face: flux=B*u_b, u_b=u_N/(1+Bh/2)
        return a
    best=0; vals={}
    for B,T,n,key in ((B1,Fs,n1,'s'),(B2,Fp-Fs,n2,'p')):
        dt=T/n; a=ab(B,dt)
        for i in range(n):
            u=solve_banded((1,1),a,u); ub=u[-1]/(1+B*h/2); best=max(best,1-ub)
        vals[key]=1-u[-1]/(1+B*h/2)
    return vals['s'],vals['p'],best
out=[]
sel=[]
for g in ['G1','G2','G3','G4']:
    pool=[r for r in R if r['group']==g and r['Fs']>=2e-2 and r['Fp']<=3 and max(r['box'][0][0],r['box'][1][0])<=60]
    sel+=pool[:4]
for r in sel:
    B1,B2=r['box'][0][0],r['box'][1][0]; Fs,Fp=r['Fs'],r['Fp']
    (lo,hi),ths,thp,M=e3.theta_peak(B1,B2,Fs,Fp,1e-10)
    Ns=[50,100,200,400]; res=[fv(B1,B2,Fs,Fp,N) for N in Ns]
    ps=[max(a,b) for a,b,_ in res]; hist=[m for *_,m in res]
    d=[ps[i]-ps[i+1] for i in range(3)]
    order=[math.log2(abs(d[i])/abs(d[i+1])) if d[i+1]!=0 and d[i]!=0 else None for i in range(2)]
    S=2.0**96; cert=(float(lo)/S,float(hi)/S)
    out.append(dict(seed=r['seed'],group=r['group'],B1=B1,B2=B2,Fs=Fs,Fp=Fp,cert=cert,fv=[float(x) for x in ps],hist=[float(x) for x in hist],order=[None if o is None else float(o) for o in order],
        inside=bool(cert[0]-1e-12<=ps[-1]<=cert[1]+1e-12), err=float(ps[-1]-0.5*(cert[0]+cert[1])), peak_at_switch_or_end=bool(abs(hist[-1]-ps[-1])<1e-12)))
    print(r['seed'],r['group'],f"B=({B1:.3g},{B2:.3g}) Fs={Fs:.3g} Fp={Fp:.3g} cert=[{cert[0]:.10f},{cert[1]:.10f}] fv400={ps[-1]:.10f} err={out[-1]['err']:.2e} ord={order} histmax-ends={hist[-1]-ps[-1]:.1e}",flush=True)
json.dump(out,open('fv16.json','w'),indent=1)
