"""Independent validation of kernel enclosures against exact-rational Taylor evaluation (different algorithm:
no argument halving, no fixed-point, Fractions with 120+ terms and Lagrange remainder), incl. near zeros, large roots."""
import math, random
from fractions import Fraction as F
import ivk as k
import e3
PIq=F(k.PI[0]+k.PI[1],2*k.S); PIe=F(k.PI[1]-k.PI[0],2*k.S)   # center; use independent pi below
def pi_exact(terms=400):                      # Gauss-Legendre-free: Machin with big rationals, tail rigorous
    def at(x,n): return sum(F((-1)**i,(2*i+1)*x**(2*i+1)) for i in range(n)),F(1,(2*n+1)*x**(2*n+1))
    a,ta=at(5,150); b,tb=at(239,60); return 16*a-4*b,16*ta+4*tb
PIR,PIT=pi_exact()
def q_exp_neg(t,N=None):
    N=int(3*t)+300   # exp(-t), t rational >=0, direct Taylor, alternating? use exp(t) series and reciprocal for rigor
    s=F(0); term=F(1)
    for i in range(N): s+=term; term=term*t/(i+1)
    rem=term*2     # tail <= 2*next term when t/(N)<1/2 ; valid if N>2t -> use N large
    assert t<F(N,2)
    return F(1)/(s+rem),F(1)/s          # [lo,hi] for exp(-t)=1/exp(t), exp(t) in [s,s+rem]
def q_cos(x,N=120):
    s=F(0); term=F(1)
    for i in range(N): s+=term; term=-term*x*x/((2*i+1)*(2*i+2))
    r=abs(term)*2; return s-r,s+r
def q_sin(x,N=120):
    s=F(0); term=x
    for i in range(N): s+=term; term=-term*x*x/((2*i+2)*(2*i+3))
    r=abs(term)*2; return s-r,s+r
def ov(iv_,lo,hi): return F(iv_[0],k.S)<=hi and lo<=F(iv_[1],k.S)       # overlap (both contain truth)
def sub_(iv_,lo,hi): return F(iv_[0],k.S)<=lo and hi<=F(iv_[1],k.S)       # kernel encloses the narrow rational enclosure
random.seed(7); fails=0; n=0; notsub=0
pts=[0.0,1e-9,0.5,1.0,math.pi/2-1e-9,math.pi/2,math.pi/2+1e-9,math.pi-1e-9,math.pi,3.2]+[random.uniform(0,3.2) for _ in range(150)]
for x in pts:
    X=k.fl(x); xf=F(x)
    for kern,ref in ((k.cosT,q_cos),(k.sinT,q_sin)):
        lo,hi=ref(xf); n+=1; r=kern(X); fails+=not ov(r,lo,hi); notsub+=not sub_(r,lo,hi)
    lo,hi=q_sin(xf); lo2,hi2=(F(1),F(1)) if x==0 else (lo/xf,hi/xf)
    r=k.sinc(X); n+=1; fails+=not ov(r,lo2,hi2); notsub+=not sub_(r,lo2,hi2)
for y in [0.0,1e-12,0.3,1.0,7.7,30.0,99.9,100.0,150.0,400.0]+[random.uniform(0,120) for _ in range(120)]:
    r=k.exp_neg(k.fl(y)); lo,hi=q_exp_neg(F(y)); n+=1
    fails+=not ov(r,lo,hi); notsub+=not sub_(r,lo,hi) if y<100 else 0
print(f"transcendentals: {n} checks, overlap failures {fails}, kernel-not-enclosing-narrow-rational {notsub}")
# pi
print("pi: kernel encloses independent rational pi:",F(k.PI[0],k.S)<=PIR-PIT and PIR+PIT<=F(k.PI[1],k.S))
# roots incl large spectral indices & large B: validate bracket against exact rational sign test with independent Taylor
bad=0; tot=0; maxw=0
for B in (1e-9,0.05,1.0,7.3,40.0,500.0):
    for m in (0,1,2,10,50,300,600):
        d=e3.root_delta(m,B); a,b=F(d[0],k.S),F(d[1],k.S); maxw=max(maxw,float(b-a)); tot+=1
        def g(dl):                          # exact-rational g via independent Taylor (take midpoints of its tiny enclosure)
            sl,sh=q_sin(dl); cl,ch=q_cos(dl); sm=(sl+sh)/2; cm=(cl+ch)/2
            return (m*PIR+dl)*sm-F(B)*cm
        ok=g(a)<0<g(b) and b<PIR/2 and a>0
        bad+=not ok
print(f"root brackets: {tot} validated independently (rational Taylor sign test), failures {bad}, max bracket width {maxw:.1e}")
