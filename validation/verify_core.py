"""Independent numerical verification of the proof chain (Lemmas 3.2-3.4, Thm 3.5). double precision, K modes."""
import numpy as np, math
from e4 import closed
rng=np.random.default_rng(20261122)
def modes(B,K): z=closed.eigz(B,K); return z,closed.Nn(z)
def G(B,x,y,t,K=400):
    z,N=modes(B,K); return (np.cos(np.outer(x,z))[:,None,:]*np.cos(np.outer(y,z))[None,:,:]*np.exp(-z**2*t)/N).sum(-1)
# --- Lemma 3.2: TP2 minors of G_B
bad=0;minrel=1;n=0
for _ in range(4000):
    B=10**rng.uniform(-2,2) if rng.random()>.1 else 0.0
    t=10**rng.uniform(-2.5,0.5); x=np.sort(rng.random(2)); y=np.sort(rng.random(2))
    g=G(B,x,y,t,K=int(min(600,40+8/math.sqrt(t))))
    D=g[0,0]*g[1,1]-g[0,1]*g[1,0]; sc=g[0,0]*g[1,1]+g[0,1]*g[1,0]
    if sc<=1e-280: continue
    n+=1; minrel=min(minrel,D/sc)
    if D/sc<-1e-9: bad+=1
print("TP2 minors tested",n,"violations",bad,"min rel",minrel)
# --- Lemma 3.3: R(y,t) nondecreasing in t
bad=0;n=0
for _ in range(600):
    B=10**rng.uniform(-2,2) if rng.random()>.1 else 0.0
    y=rng.random(); ts=np.logspace(-3,0.5,60)
    R=np.array([G(B,[y],[1.0],t,K=int(min(700,40+8/math.sqrt(t))))[0,0]/G(B,[1.0],[1.0],t,K=int(min(700,40+8/math.sqrt(t))))[0,0] for t in ts])
    n+=1
    if np.any(np.diff(R)<-1e-10*np.maximum(1,R[1:])): bad+=1
print("ratio-monotone profiles tested",n,"violations",bad)
# --- Lemma 3.4 identity and Thm 3.5(iii): sign pattern of theta_s' after the switch
def two(B1,B2,Fs,K=700):
    z,N1=modes(B1,K); w,N2=modes(B2,K); c=np.sinc(z/np.pi)/N1
    Z,W=np.meshgrid(z,w,indexing='ij'); d=Z**2-W**2
    with np.errstate(all='ignore'): I=(B1-B2)*np.cos(Z)*np.cos(W)/d
    diag=0.5*(np.sinc((Z-W)/np.pi)+np.sinc((Z+W)/np.pi)); I=np.where(np.abs(Z-W)<1,diag,I)
    A=((c*np.exp(-z**2*Fs))[:,None]*I).sum(0)/N2
    return z,w,c,N1,N2,A
errs=[];pat={}
for _ in range(1500):
    B1=10**rng.uniform(-1,2.3); B2=0.0 if rng.random()<.15 else 10**rng.uniform(-2,B1 and math.log10(B1)+0.5)
    Fs=10**rng.uniform(-2,0.3); K=int(min(900,60+10/math.sqrt(Fs)))
    z,w,c,N1,N2,A=two(B1,B2,Fs,K=K)
    taus=np.logspace(-4,1.5,400)
    dth=np.array([ (A*np.cos(w)*w**2*np.exp(-w**2*t)).sum() for t in taus])
    # identity check at a few taus
    if _%50==0:
        y=np.linspace(0,1,4001); u0=(c*np.exp(-z**2*Fs))@np.cos(np.outer(z,y)); p=(c*z**2*np.exp(-z**2*Fs))@np.cos(np.outer(z,y))
        u01=u0[-1]
        for t in (0.01,0.1,0.5):
            Gy=(np.cos(w)/N2*np.exp(-w**2*t))@np.cos(np.outer(w,y)); G11=(np.cos(w)**2/N2*np.exp(-w**2*t)).sum()
            f=p*Gy; rhs=np.trapezoid(f,y)-(B1-B2)*u01*G11
            lhs=(A*np.cos(w)*w**2*np.exp(-w**2*t)).sum()
            errs.append(abs(lhs-rhs)/(abs(lhs)+abs(rhs)+1e-30))
    s=np.sign(dth[np.abs(dth)>1e-12*np.max(np.abs(dth))]); ch=int(np.sum(s[1:]!=s[:-1]))
    key=("B1<=B2" if B1<=B2 else "B1>B2",ch,"first sign "+("-" if s[0]<0 else "+"))
    pat[key]=pat.get(key,0)+1
print("identity max rel err",max(errs),len(errs))
for k,v in sorted(pat.items()): print(k,v)
