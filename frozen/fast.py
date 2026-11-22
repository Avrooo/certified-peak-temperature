import math, time, numpy as np
from math import exp, erfc, sqrt, pi
SQPI=sqrt(pi); ISQPI=1/SQPI
def erfcx(x): return exp(x*x)*erfc(x)
def dg(s):
    """image part of the constant-flux surface response: g(s)-2sqrt(s/pi), exact up to tails <1e-15"""
    if s<=0: return 0.0
    if s<0.2:
        r=sqrt(s); tot=0.0; m=1
        while True:
            z=m/r
            if z*z>40: break
            tot+=exp(-z*z)*(ISQPI-z*erfcx(z)); m+=1
        return 4*r*tot
    t=0.0
    for n in range(1,7): t+=exp(-n*n*pi*pi*s)/(n*n)
    return s+1/3-(2/(pi*pi))*t-2*sqrt(s/pi)
def Lsemi(Bi,F): return 1-erfcx(Bi*sqrt(F))
def bounds(Bi,F,N):
    """returns (lower, upper). N=0: closed form. N>0: N uniform segments on the smooth image part"""
    L=Lsemi(Bi,F); dgF=dg(F)
    U0=min(1.0,L+Bi*dgF)
    if N==0: return L,U0
    t=[F*j/N for j in range(N+1)]
    d=[dg(F-tj) for tj in t]; d[-1]=0.0
    inc=[d[j]-d[j+1] for j in range(N)]                       # int of image kernel over segment j
    U=L+Bi*sum(erfcx(Bi*sqrt(t[j]))*inc[j] for j in range(N)) # flux <= Bi*erfcx at left end
    # lower: D>=int q_low k_img - Bi int Dup k_semi
    qlow=[Bi*(1-min(1.0,Lsemi(Bi,t[j+1])+Bi*dg(t[j+1]))) for j in range(N)]
    A=sum(qlow[j]*inc[j] for j in range(N))
    kap=lambda s:2*sqrt(max(s,0)/pi)
    B=Bi*sum(Bi*dg(t[j+1])*(kap(F-t[j])-kap(F-t[j+1])) for j in range(N))
    lo=L+max(0.0,A-B)
    return lo,min(U,U0)
if __name__=="__main__":
    import sys; sys.path.insert(0,'.')
    from pilot_lib import roots, step
    Bi,F=0.1,0.1; z,C=roots(Bi,300); S=step(Bi,F,z,C,1.0)[0]
    print("worked example exact %.9f"%S)
    for N in (0,1,2,4,8):
        lo,up=bounds(Bi,F,N); print(" N=%d  lo %.9f  up %.9f  rel width %.2e"%(N,lo,up,(up-lo)/S))
    rng=np.random.default_rng(1); M=400
    Bu=10**rng.uniform(np.log10(.005),np.log10(.5),M); Fu=10**rng.uniform(np.log10(.005),np.log10(5),M)
    Sx=np.array([step(b,f,*roots(b,300),1.0)[0] for b,f in zip(Bu,Fu)])
    print("\nunseen 400 cases (seed 1):  N  violations  median W_r   90th pct W_r   max W_r   us/query")
    for N in (0,1,2,4,8,16):
        res=np.array([bounds(b,f,N) for b,f in zip(Bu,Fu)])
        viol=np.sum((Sx<res[:,0]-1e-12)|(Sx>res[:,1]+1e-12)); W=(res[:,1]-res[:,0])/Sx
        t0=time.perf_counter()
        for _ in range(5):
            for b,f in zip(Bu,Fu): bounds(b,f,N)
        us=(time.perf_counter()-t0)/5/M*1e6
        print("   %2d   %d      %.2e     %.2e      %.3f    %6.1f"%(N,viol,np.median(W),np.percentile(W,90),W.max(),us))
    # by pulse length
    res=np.array([bounds(b,f,4) for b,f in zip(Bu,Fu)]); W=(res[:,1]-res[:,0])/Sx
    for name,m in [("Fo<0.1",Fu<0.1),("0.1-0.3",(Fu>=0.1)&(Fu<0.3)),("0.3-1",(Fu>=0.3)&(Fu<1)),("Fo>=1",Fu>=1)]:
        print("  N=4, %-8s median W_r %.2e  max %.3f"%(name,np.median(W[m]),W[m].max()))
