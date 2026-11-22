import math, time, os
from fractions import Fraction
import numpy as np
from scipy.optimize import brentq
import ivk as k
from ivk import fl,add,sub,mul,div,divi,muli,sq,neg,cn,PI,ONE,ZERO,exp_neg,sinc,cosT,sinT,hull

NMIN=cn(42,100)                      # Lemma G1: N>=0.4204>0.42
class CertFail(Exception): pass
def _g(m,d,B):                       # (m pi + d) sin d - B cos d  (d,B point intervals)
    return sub(mul(add(muli(PI,m),d),sinT(d)),mul(B,cosT(d)))
_rootcache={}
def root_delta(m,B):
    """validated enclosure of delta=zeta_m-m*pi in (0,pi/2) for zeta tan zeta=B (B float>0).
    sign(g)=sign(h) on (0,pi/2) (cos>0) and h=(m pi+d)tan d is strictly increasing => unique root in any bracket with g(a)<0<g(b)."""
    key=(m,B)
    if key in _rootcache: return _rootcache[key]
    Bi=fl(B); seed=brentq(lambda d:(m*math.pi+d)*math.tan(d)-B,1e-300,math.pi/2-1e-13,xtol=1e-17,rtol=8.9e-16)
    w=2e-15*(1+seed)
    for _ in range(12):
        a=seed-w
        if a<=0: a=seed*0.5
        b=min(seed+w,math.pi/2-1e-12)
        A,Bb=fl(a),fl(b)
        ga=_g(m,A,Bi); gb=_g(m,Bb,Bi)
        if ga[1]<0 and gb[0]>0 and Bb[1]<k.HALFPI[0]:
            _rootcache[key]=(A[0],Bb[1]); return _rootcache[key]
        w*=6
    raise CertFail(f"root bracket m={m} B={B}")
class Mode:
    __slots__=("z","d","cz","sz","N","w2","zero")
def make_mode(B,m):
    md=Mode()
    if B==0:
        md.d=(0,0); md.z=muli(PI,m); md.cz=cn((-1)**m); md.sz=ZERO
        md.N=ONE if m==0 else cn(1)   # sinc(2 m pi)=0 for m>=1 => N=1/2 ; fixed below
        if m>0: md.N=cn(1,2)
        md.zero=(m==0); return md
    d=root_delta(m,B); md.d=d; md.zero=False
    md.z=add(muli(PI,m),d); s=sinT(d); c=cosT(d); sg=(-1)**m
    md.cz=muli(c,sg); md.sz=muli(s,sg)
    md.N=divi(add(ONE,div(mul(s,c),md.z)),2)    # N=1/2(1+ sin2z/(2z)) = 1/2(1+ s c / z)
    return md
_modecache={}
def modes(B,M):
    out=[]
    for m in range(M):
        key=(B,m)
        if key not in _modecache: _modecache[key]=make_mode(B,m)
        out.append(_modecache[key])
    return out
def coef(md):                         # c_m = sinc(z)/N
    return div(ONE,md.N) if md.zero else div(div(md.sz,md.z),md.N)
# ---------- Lemma G tail (all interval/upward) ----------
def S_tail(t,M):                      # upper bound on sum_{j>=M} exp(-(j pi)^2 t), t interval>0
    pi2=sq(PI); a=exp_neg(mul(mul(cn(M*M),pi2),t)); r=exp_neg(mul(mul(cn(2*M+1),pi2),t))
    den=sub(ONE,r)
    if den[0]<=0: raise CertFail("geometric ratio not <1")
    return div(a,den)[1]
def EM_two(B1,B2,Fs,Fp,M):
    Fs_,Fp_=fl(Fs),fl(Fp); tau=sub(Fp_,Fs_); B1i=fl(B1)
    pi2=sq(PI); dB=sub(fl(B1),fl(B2)); dB=(0,max(abs(dB[0]),abs(dB[1])))
    f=div(muli(dB,2),mul(pi2,cn(M))); fM=min(f[1],S_)
    Mp2=mul(cn(M*M),pi2); base=mul(B1i,div(ONE,mul(sq(NMIN),Mp2)))[1]   # B1/(Nmin^2 (M pi)^2)
    S_s=S_tail(Fs_,M); S_t1=S_tail(tau,1); S_p=S_tail(Fp_,M); S_tM=S_tail(tau,M)
    up=lambda x,y:mul((0,x),(0,y))[1]
    offA=up(up(base,fM),up(S_s,S_+S_t1))        # S_ = ONE (scaled): 1+S(tau,1)
    diag=up(base,S_p)
    cm=add(ONE,divi(B1i,6))[1]
    nm2=div(ONE,sq(NMIN))[1]
    B_=up(up(up(fM,nm2),cm),S_tM)
    return offA+diag+B_
S_=k.S
def EM_one(B1,Fs,M):
    Fs_=fl(Fs); B1i=fl(B1); Mp2=mul(cn(M*M),sq(PI))
    return mul((0,mul(B1i,div(ONE,mul(NMIN,Mp2)))[1]),(0,S_tail(Fs_,M)))[1]
def EM_float(B1,B2,Fs,Fp,M,one=False):
    tau=Fp-Fs; Nm=0.42; f=min(1,2*abs(B1-B2)/(math.pi**2*M))
    S=lambda t,M:math.exp(-(M*math.pi)**2*t)/(1-math.exp(-(2*M+1)*math.pi**2*t))
    if one: return B1/(Nm*(M*math.pi)**2)*S(Fs,M)
    b=B1/(Nm**2*(M*math.pi)**2)
    return b*f*S(Fs,M)*(1+S(tau,1))+b*S(Fp,M)+f/Nm**2*(1+B1/6)*S(tau,M)
def pick_M(fun,eps,Mmax=2000):
    M=1
    while M<Mmax and fun(M)>eps: M+=1 if M<20 else max(1,M//8)
    return M
# ---------- certified endpoint evaluations ----------
def u_single(B1,Fs,eps,fixedM=None):
    if fixedM is None:
        M=pick_M(lambda M:EM_float(B1,0,Fs,Fs,M,one=True),eps/4)
        while True:
            E=EM_one(B1,Fs,M)
            if E<=k.fl(eps)[1]: break
            M=int(M*1.3)+1
    else: M=fixedM; E=EM_one(B1,Fs,M)
    md=modes(B1,M); Fi=fl(Fs); s=ZERO
    for q in md: s=add(s,mul(mul(coef(q),q.cz),exp_neg(mul(sq(q.z),Fi))))
    return (s[0]-E,s[1]+E),M
def u_two(B1,B2,Fs,Fp,eps,fixedM=None):
    if fixedM is None:
        M=pick_M(lambda M:EM_float(B1,B2,Fs,Fp,M),eps/4)
        while True:
            E=EM_two(B1,B2,Fs,Fp,M)
            if E<=k.fl(eps)[1]: break
            M=int(M*1.3)+1
    else: M=fixedM; E=EM_two(B1,B2,Fs,Fp,M)
    m1=modes(B1,M); m2=modes(B2,M); Fi=fl(Fs); tau=sub(fl(Fp),Fi); dB=sub(fl(B1),fl(B2))
    a=[mul(coef(q),exp_neg(mul(sq(q.z),Fi))) for q in m1]
    dn=[mul(div(q.cz,q.N),exp_neg(mul(sq(q.z),tau))) for q in m2]
    tot=ZERO
    for i,q in enumerate(m1):
        row=ZERO
        for j,r in enumerate(m2):
            if i!=j:
                den=mul(sub(q.z,r.z),add(q.z,r.z))
                I=div(mul(mul(dB,q.cz),r.cz),den)
            else:
                zm=sub(q.d,r.d); zp=add(q.z,r.z)
                if i==0: sp=sinc(zp)
                else: sp=div(sinT(add(q.d,r.d)),zp)
                I=divi(add(sinc(zm),sp),2)
            row=add(row,mul(I,dn[j]))
        tot=add(tot,mul(a[i],row))
    return (tot[0]-E,tot[1]+E),M
def theta_peak(B1,B2,Fs,Fp,eps):
    u1,M1=u_single(B1,Fs,eps); u2,M2=u_two(B1,B2,Fs,Fp,eps)
    ths=sub(ONE,u1); thp=sub(ONE,u2)
    return (max(ths[0],thp[0]),max(ths[1],thp[1])),ths,thp,max(M1,M2)
