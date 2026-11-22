"""Level-2 audit harness: differential tests of every kernel operation against exact Fractions, plus MUTATION tests
(deliberately broken rounding) to show the tests can detect enclosure violations."""
import random, math, types, os, sys
from fractions import Fraction as F
src=open('ivk.py').read()
def load(s):
    m=types.ModuleType('ivm'); m.__dict__['__name__']='ivm'; exec(compile(s,'ivm','exec'),m.__dict__); return m
def rnd_iv(k,rng,scale=1.0):
    a=F(rng.randint(-10**9,10**9),rng.randint(1,10**9))*F(scale)*rng.choice([1,1,1e-6,1e3])
    b=a+abs(F(rng.randint(0,10**6),rng.randint(1,10**9)))*rng.choice([0,1,1e-9,1e-3])
    iv=(math.floor(a*k.S),math.ceil(b*k.S)); return (a,b),iv
def contains(iv,k,lo,hi): return F(iv[0],k.S)<=lo and hi<=F(iv[1],k.S)
def run(k,n=4000,seed=1):
    rng=random.Random(seed); bad={}
    def chk(name,ok): bad[name]=bad.get(name,0)+(not ok)
    for _ in range(n):
        (a,b),A=rnd_iv(k,rng); (c,d),B=rnd_iv(k,rng)
        ea,eb=F(A[0],k.S),F(A[1],k.S); ec,ed=F(B[0],k.S),F(B[1],k.S)   # exact values the kernel intervals denote
        chk('add',contains(k.add(A,B),k,ea+ec,eb+ed)); chk('sub',contains(k.sub(A,B),k,ea-ed,eb-ec))
        p=[ea*ec,ea*ed,eb*ec,eb*ed]; chk('mul',contains(k.mul(A,B),k,min(p),max(p)))
        sqr=[ea*ea,eb*eb]; lo=0 if ea<=0<=eb else min(sqr); chk('sq',contains(k.sq(A),k,lo,max(sqr)))
        kk=rng.randint(1,10**6); chk('divi',contains(k.divi(A,kk),k,ea/kk,eb/kk))
        chk('muli',contains(k.muli(A,-kk),k,-eb*kk,-ea*kk))
        if ec>0 or ed<0:
            q=[ea/ec,ea/ed,eb/ec,eb/ed]; chk('div',contains(k.div(A,B),k,min(q),max(q)))
        x=F(rng.randint(-10**12,10**12),rng.randint(1,10**12)); X=k.fl(x)
        chk('fl',F(X[0],k.S)<=x<=F(X[1],k.S))
    return bad
if __name__=="__main__":
    k=load(src); print("REAL kernel   :",run(k))
    mutants={
     "mul hi uses floor":("-((-max(p))>>P)","(max(p)>>P)"),
     "mul lo uses ceil":("return (min(p)>>P, -((-max(p))>>P))","return (-((-min(p))>>P), -((-max(p))>>P))"),
     "div hi uses floor":("hi=[-((-(x<<P))//y) for x in a for y in b]","hi=[((x<<P)//y) for x in a for y in b]"),
     "divi hi uses floor":("return (a[0]//k, -((-a[1])//k))","return (a[0]//k, a[1]//k)"),
     "sq straddle lo wrong":("return (0,-((-m*m)>>P))","return (m*m>>P,-((-m*m)>>P))"),
     "fl uses round":("return (math.floor(f), math.ceil(f))","return (round(f), round(f))"),
    }
    for name,(old,new) in mutants.items():
        assert old in src,name
        km=load(src.replace(old,new,1)); r=run(km); print(f"MUTANT {name:24s} detected:",any(v>0 for v in r.values()), {a:b for a,b in r.items() if b})
