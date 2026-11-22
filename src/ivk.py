"""Minimal validated interval kernel: fixed-point big ints (scale 2^P), outward rounding.
Every op returns an interval guaranteed to contain the exact result for all inputs in the input intervals."""
from fractions import Fraction
import math
import os
P=int(os.environ.get("IVK_P","96")); S=1<<P
KT=int(os.environ.get("IVK_K","24")); EXPN=int(os.environ.get("IVK_N","24")); OPT={"exp1":os.environ.get("IVK_EXP1","1")=="1"}
def fl(x):                      # exact float/Fraction/int -> point interval
    f=Fraction(x)*S; return (math.floor(f), math.ceil(f))
def cn(a,b=1): return fl(Fraction(a,b))
def add(a,b): return (a[0]+b[0],a[1]+b[1])
def sub(a,b): return (a[0]-b[1],a[1]-b[0])
def neg(a): return (-a[1],-a[0])
def mul(a,b):
    a0,a1=a; b0,b1=b
    if a0>=0 and b0>=0: return ((a0*b0)>>P, -((-(a1*b1))>>P))
    if a1<=0 and b1<=0: return ((a1*b1)>>P, -((-(a0*b0))>>P))
    p=(a0*b0,a0*b1,a1*b0,a1*b1); return (min(p)>>P, -((-max(p))>>P))
def muli(a,k): 
    p=(a[0]*k,a[1]*k); return (min(p),max(p))
def divi(a,k):                  # k>0 int
    return (a[0]//k, -((-a[1])//k))
def div(a,b):
    if not (b[0]>0 or b[1]<0): raise ZeroDivisionError("division by interval containing 0")
    lo=[(x<<P)//y for x in a for y in b]; hi=[-((-(x<<P))//y) for x in a for y in b]
    return (min(lo),max(hi))
def sq(a):
    if a[0]>=0: return (a[0]*a[0]>>P, -((-a[1]*a[1])>>P))
    if a[1]<=0: return (a[1]*a[1]>>P, -((-a[0]*a[0])>>P))
    m=max(-a[0],a[1]); return (0,-((-m*m)>>P))
def absmax(a): return max(abs(a[0]),abs(a[1]))
def widen(a,e): return (a[0]-e,a[1]+e)
def hull(a,b): return (min(a[0],b[0]),max(a[1],b[1]))
def tof(a): return (a[0]/S,a[1]/S)
ONE=cn(1); ZERO=(0,0)
# ---- pi (Machin) ----
def _atan_inv(x,terms):
    s=Fraction(0)
    for k in range(terms): s+=Fraction((-1)**k,(2*k+1)*x**(2*k+1))
    tail=Fraction(1,(2*terms+1)*x**(2*terms+1)); return s,tail
_a,ta=_atan_inv(5,100); _b,tb=_atan_inv(239,40)
_pi=16*_a-4*_b; _pe=16*ta+4*tb
PI=(math.floor((_pi-_pe)*S), math.ceil((_pi+_pe)*S))
HALFPI=divi(PI,2)
# ---- exp(-t) for exact point t>=0 (scaled int) ----
def _exp_neg_pt(t):
    j=0
    while (t>>j)>S//2: j+=1
    r=(t>>j, -((-t)>>j))
    s=ONE; term=ONE; n=EXPN
    for k in range(1,n+2):
        term=divi(mul(term,neg(r)),k)
        if k<=n: s=add(s,term)
    e=absmax(term); s=(s[0]-e,s[1]+e)
    for _ in range(j): s=mul(s,s)
    return s
_E100=None
def exp_neg(y):                 # y interval, y[0]>=0 ; returns enclosure of exp(-y)
    global _E100
    if y[0]<0: raise ValueError("exp_neg needs y>=0")
    def pt(t):
        if t>=100*S:
            global _E100
            if _E100 is None: _E100=_exp_neg_pt(100*S)
            return (0,_E100[1])
        return _exp_neg_pt(t)
    if OPT["exp1"]:
        d=y[1]-y[0]; l,h=pt(y[0]) if y[0]>0 else (S,S)
        lo=0 if d>=S else (l*(S-d))>>P
        return (max(0,lo), min(S,h))
    return (max(0,pt(y[1])[0]), min(S,pt(y[0])[1]) if y[0]>0 else S)
# ---- sinc, cos by Taylor with explicit remainder (valid for any interval, intended |x|<=4) ----
def sinc(x,K=None):
    K=K or KT
    x2=sq(x); s=ONE; term=ONE
    for k in range(1,K+2):
        term=divi(mul(term,neg(x2)),(2*k)*(2*k+1))
        if k<=K: s=add(s,term)
    e=absmax(term); return (s[0]-e,s[1]+e)
def cosT(x,K=None):
    K=K or KT
    x2=sq(x); s=ONE; term=ONE
    for k in range(1,K+2):
        term=divi(mul(term,neg(x2)),(2*k-1)*(2*k))
        if k<=K: s=add(s,term)
    e=absmax(term); return (s[0]-e,s[1]+e)
def sinT(x): return mul(x,sinc(x))
