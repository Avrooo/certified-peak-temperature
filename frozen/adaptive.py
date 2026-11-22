import math, time, numpy as np
from math import exp, erfc, sqrt, pi, sin, cos, atan
from fast import bounds as seg_bounds, Lsemi, erfcx, dg
SQPI=sqrt(pi)
# ---------- short-pulse branch: image series with rigorous remainder ----------
def ierfc_ub(z):            # ierfc(z) <= e^{-z^2}/sqrt(pi) * (sqrt(z^2+2)-z)/(z+sqrt(z^2+2))   (from erfcx(z)>2/(sqrt(pi)(z+sqrt(z^2+2))))
    r=sqrt(z*z+2); return exp(-z*z)/SQPI*(r-z)/(z+r)
def img_partial(F,M):
    """(lower, upper) for Delta g(F)=4 sqrt(F) sum_{m>=1} ierfc(m/sqrt F): first M terms exact (skipped if z^2>40), rest bounded"""
    r=sqrt(F); s=0.0; m=1
    while m<=M:
        z=m/r
        if z*z>40: break
        s+=exp(-z*z)*(1/SQPI-z*erfcx(z)); m+=1
    z=m/r; q=exp(-(2*m+1)/F)
    tail=ierfc_ub(z)/(1-q)
    return 4*r*s, 4*r*(s+tail)
# ---------- long-pulse branch: modal series, all terms positive ----------
def root(k,Bi):
    """root of z tan z = Bi in (k pi, k pi + pi/2): psi(z)=z-k*pi-atan(Bi/z) is increasing, so safeguarded Newton stays bracketed"""
    base=k*pi; a=base; b=base+pi/2
    z=min(max(sqrt(Bi),1e-9),b-1e-9) if k==0 else a+atan(Bi/(k*pi))
    z=min(max(z,a+1e-12),b-1e-12)
    for _ in range(8):
        psi=z-base-atan(Bi/z)
        if psi>0: b=z
        else: a=z
        d=1+Bi/(z*z+Bi*Bi)
        zn=z-psi/d
        if abs(psi)<1e-15: break
        z=zn if a<=zn<=b else 0.5*(a+b)
    return z
def mode(k,Bi,F):
    z=root(k,Bi); return 4*sin(z)*cos(z)/(2*z+sin(2*z))*exp(-z*z*F)      # a_k>0 ; S=1-sum a_k
def modal_tail(M,F):         # a_k <= e^{-k^2 pi^2 F}/(k pi) for k>=1; sum over k>=M
    return exp(-M*M*pi*pi*F)/(M*pi)/(1-exp(-(2*M+1)*pi*pi*F))
FSW=0.25
def decide(Bi,F,th,stats=None):
    if F<FSW:
        L=Lsemi(Bi,F)
        if L>=th: return 'exceed','L'
        for M in range(0,5):
            lo_d,up_d=img_partial(F,M); up=min(1.0,L+Bi*up_d)
            if up<th: return 'safe','img%d'%M
        for N in (2,4,8):
            lo,upN=seg_bounds(Bi,F,N)
            if upN<th: return 'safe','seg%d'%N
            if lo>=th: return 'exceed','seg%d'%N
        return decide_modal(Bi,F,th)
    return decide_modal(Bi,F,th)
def decide_modal(Bi,F,th):
    s=0.0
    for k in range(0,int(sqrt(45/F)/pi)+4):
        s+=mode(k,Bi,F); up=1-s
        if up<th: return 'safe','modal%d'%(k+1)
        lo=up-modal_tail(k+1,F)
        if lo>=th: return 'exceed','modal%d'%(k+1)
    return 'undecided','modal'
def eig_scalar(Bi,F):
    M=int(sqrt(30/F)/pi)+3
    return 1-sum(mode(k,Bi,F) for k in range(M))
if __name__=="__main__":
    from pilot_lib import roots, step
    for Bi in (0.005,0.1,0.5):
        z,C=roots(Bi,10)
        assert max(abs(root(k,Bi)-z[k]) for k in range(6))<1e-11,Bi
    rng=np.random.default_rng(7); M=300
    Bu=10**rng.uniform(np.log10(.005),np.log10(.5),M); Fu=10**rng.uniform(np.log10(.005),np.log10(5),M)
    S=np.array([step(b,f,*roots(b,300),1.0)[0] for b,f in zip(Bu,Fu)])
    # validity of bracket stages
    bad=0
    for b,f,s in zip(Bu,Fu,S):
        if f>=FSW:
            acc=0.0
            for k in range(6):
                acc+=mode(k,b,f); up=1-acc; lo=up-modal_tail(k+1,f)
                bad+= (s>up+1e-12) or (s<lo-1e-12)
        else:
            L=Lsemi(b,f)
            for Mi in range(5):
                lo_d,up_d=img_partial(f,Mi); bad+= s>min(1.0,L+b*up_d)+1e-12
    print("bracket stage violations:",bad)
    def tm(fn,reps=9):
        fn();fn();ts=[]
        for _ in range(reps):
            t0=time.perf_counter();fn();ts.append(time.perf_counter()-t0)
        return np.median(ts)/M*1e6
    print("eigen (modal, fixed M): %.1f us"%tm(lambda:[eig_scalar(b,f) for b,f in zip(Bu,Fu)]))
    print("\n r=th/S   decided  wrong  mean_us   stages used")
    for r in (0.2,0.5,0.9,0.98,0.999,1.001,1.02,1.1,2,5):
        th=r*S; res=[decide(b,f,t) for b,f,t in zip(Bu,Fu,th)]
        truth=['exceed' if s>=t else 'safe' for s,t in zip(S,th)]
        wrong=sum(x[0]!='undecided' and x[0]!=tr for x,tr in zip(res,truth)); dec=np.mean([x[0]!='undecided' for x in res])
        from collections import Counter; c=Counter(x[1] for x in res)
        print(" %6.3f  %5.1f%%   %d   %6.1f   %s"%(r,100*dec,wrong,tm(lambda:[decide(b,f,t) for b,f,t in zip(Bu,Fu,th)]),dict(c.most_common(5))))
    for lo,hi,name in [(0.005,0.1,"Fo<0.1"),(0.1,0.25,"0.1-0.25"),(0.25,1,"0.25-1"),(1,5.01,"Fo>1")]:
        m=(Fu>=lo)&(Fu<hi); B2,F2,S2=Bu[m],Fu[m],S[m]
        th=0.98*S2
        def a(): [decide(b,f,t) for b,f,t in zip(B2,F2,th)]
        def e(): [eig_scalar(b,f) for b,f in zip(B2,F2)]
        def tt(fn):
            fn();fn();ts=[]
            for _ in range(9):
                t0=time.perf_counter();fn();ts.append(time.perf_counter()-t0)
            return np.median(ts)/len(B2)*1e6
        print("  r=0.98 %-9s n=%3d adaptive %.1f us  eigen %.1f us"%(name,m.sum(),tt(a),tt(e)))
