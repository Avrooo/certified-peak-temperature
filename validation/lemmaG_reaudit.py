"""Independent re-audit of Lemma 4.1 (=G) after corrections: every tail piece tested separately, plus G1-G3 pointwise."""
import numpy as np, io, contextlib
from math import pi
with contextlib.redirect_stdout(io.StringIO()):
    from closed import eigz,Nn
rng=np.random.default_rng(20261121)
Nmin=0.5*(1-1/(2*pi))
def S(t,M): return np.exp(-M*M*pi*pi*t)/(1-np.exp(-(2*M+1)*pi*pi*t))
def Sexact(t,M,K=2000): k=np.arange(M,K); return np.sum(np.exp(-(k*pi)**2*t))
viol={k:0 for k in("G1","G2","G3a","G3b","S","offA","diag","RB","total")}; worst={k:0.0 for k in viol}; n=0
def rec(k,ok,ratio,bound=1.0):
    global viol,worst
    if bound<1e-10: ratio=0.0
    if not ok: viol[k]+=1
    worst[k]=max(worst[k],ratio)
cases=[]
for _ in range(1500):
    kind=rng.choice(["gen","eq","near","b1zero","b2zero","big","short"])
    B1=10**rng.uniform(-2,2.5); B2=10**rng.uniform(-2,2.5)
    ts=10**rng.uniform(-2.3,-0.3); tau=10**rng.uniform(-2.5,-0.2)
    if kind=="eq": B2=B1
    if kind=="near": B2=B1*(1+10**rng.uniform(-9,-3))
    if kind=="b1zero": B1=0.0
    if kind=="b2zero": B2=0.0
    if kind=="big": B1=10**rng.uniform(2,3.5); B2=10**rng.uniform(1,3.5)
    if kind=="short": ts=10**rng.uniform(-3.3,-2.3)
    cases.append((B1,B2,ts,tau))
K=int(max(60,min(700,np.sqrt(40/(pi*pi*3e-4)))))
for B1,B2,ts,tau in cases:
    tp=ts+tau; K=int(min(700,max(60,np.sqrt(38/(pi**2*min(ts,tau))))))
    z=eigz(B1,K); w=eigz(B2,K); N1=Nn(z); N2=Nn(w)
    c=np.sinc(z/pi)/N1
    if B1==0.0: c[1:]=0.0     # exact: sin(m*pi)=0 (G5); double-precision sinc(m*pi) is only ~1e-17
    rec("G1",bool((N1>=Nmin-1e-15).all() and (N2>=Nmin-1e-15).all()),0)
    m=np.arange(1,K); 
    if B1>0:
        r=np.abs(c[1:])/(B1/(Nmin*(m*pi)**2)); rec("G2",bool((r<=1+1e-12).all()),float(r.max()))
    else: rec("G2",bool((np.abs(c[1:])<1e-14).all()),0)
    Z,W=np.meshgrid(z,w,indexing='ij'); d=Z**2-W**2
    with np.errstate(all='ignore'): I=(B1-B2)*np.cos(Z)*np.cos(W)/d
    diag=0.5*(np.sinc((Z-W)/pi)+np.sinc((Z+W)/pi)); idx=np.abs(Z-W)<1.0
    I=np.where(idx,diag,I)
    Id=0.5*(np.sinc((Z-W)/pi)+np.sinc((Z+W)/pi))
    rec("G3a",bool((np.abs(Id)<=1+1e-12).all()),float(np.abs(Id).max()))
    # G3b: off-diagonal (m!=n) bound with f_M for several M; and separation >= pi/2
    ii,jj=np.meshgrid(np.arange(K),np.arange(K),indexing='ij'); off=ii!=jj
    ok=True; wr=0
    if (np.abs(Z-W)[off]>=pi/2-1e-12).all()==False: ok=False
    for M in (1,2,3,5,8,13):
        f=min(1.0,2*abs(B1-B2)/(pi**2*M)); mask=off&(np.maximum(ii,jj)>=M)
        if mask.any():
            v=np.abs(Id[mask]); wr=max(wr,float(v.max()/max(f,1e-300)) if f>0 else float(v.max()))
            if (v>f*(1+1e-9)+1e-13).any(): ok=False
    rec("G3b",ok,wr)
    T=(c*np.exp(-z**2*ts))[:,None]*Id*(np.cos(w)/N2*np.exp(-w**2*tau))[None,:]
    # NOTE: use the exact overlap for off-diagonal (closed form agrees with Id numerically); Id here is exact integral for all pairs
    for M in (1,2,4,8,16,32):
        if M>=K-5: continue
        f=min(1.0,2*abs(B1-B2)/(pi**2*M))
        for k_,t_ in (("S",ts),("S",tau),("S",tp)):
            rec("S",bool(Sexact(t_,M)<=S(t_,M)*(1+1e-12)+1e-300),float(Sexact(t_,M)/max(S(t_,M),1e-300)))
        mm,nn=np.meshgrid(np.arange(K),np.arange(K),indexing='ij')
        RA_off=np.abs(T[(mm>=M)&(mm!=nn)]).sum(); RA_diag=np.abs(T[(mm>=M)&(mm==nn)]).sum(); RB=np.abs(T[(mm<M)&(nn>=M)]).sum()
        Eoff=B1*f/(Nmin**2*(M*pi)**2)*S(ts,M)*(1+S(tau,1)); Ed=B1/(Nmin**2*(M*pi)**2)*S(tp,M); EB=f/Nmin**2*(1+B1/6)*S(tau,M)
        rec("offA",bool(RA_off<=Eoff*(1+1e-9)+1e-15),float(RA_off/max(Eoff,1e-300)),Eoff)
        rec("diag",bool(RA_diag<=Ed*(1+1e-9)+1e-15),float(RA_diag/max(Ed,1e-300)),Ed)
        rec("RB",bool(RB<=EB*(1+1e-9)+1e-15),float(RB/max(EB,1e-300)),EB)
        # total error vs full sum (K large)
        um=T[:M,:M].sum(); full=T.sum(); err=abs(full-um); tot=Eoff+Ed+EB
        rec("total",bool(err<=tot*(1+1e-9)+1e-15),float(err/max(tot,1e-300)),tot)
    n+=1
print("cases",n,"violations",viol); print("max ratio (actual/bound; <=1 is valid)",{k:round(v,4) for k,v in worst.items()})
