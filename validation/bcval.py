import subprocess,random,math,types
from fractions import Fraction as F
from decimal import Decimal,getcontext
getcontext().prec=130
src=open('ivk.py').read()
def load(s):
    m=types.ModuleType('ivm'); m.__dict__['__name__']='ivm'; exec(compile(s,'ivm','exec'),m.__dict__); return m
def bc(exprs):
    prog="scale=120\n"+"\n".join(exprs)+"\n"
    out=subprocess.run(['bc','-l'],input=prog,capture_output=True,text=True).stdout.replace("\\\n","").split()
    return [Decimal(x) for x in out]
def run(k,tag):
    random.seed(3); xs=[0.0,1e-9,0.7,1.0,math.pi/2,math.pi-1e-9,3.1]+[random.uniform(0,3.2) for _ in range(60)]
    ys=[0.0,1e-9,0.5,5.0,33.3,99.0]+[random.uniform(0,95) for _ in range(40)]
    def dec(x): return Decimal(F(x).numerator)/Decimal(F(x).denominator)
    ex=[]; 
    for x in xs: d=repr(Decimal(x)) if False else format(Decimal(x),'f'); ex+= [f"c({d})",f"s({d})"]
    for y in ys: d=format(Decimal(y),'f'); ex+=[f"e(-({d}))"]
    vals=bc(ex); bad=0; n=0; i=0
    tol=Decimal(10)**-100
    def inside(iv,v): return Decimal(iv[0])/Decimal(k.S)-tol<=v<=Decimal(iv[1])/Decimal(k.S)+tol
    for x in xs:
        X=k.fl(x); c,s=vals[i],vals[i+1]; i+=2
        bad+= not inside(k.cosT(X),c); bad+= not inside(k.sinT(X),s); n+=2
        if x!=0: bad+= not inside(k.sinc(X),s/Decimal(x)); n+=1
    for y in ys:
        v=vals[i]; i+=1; bad+= not inside(k.exp_neg(k.fl(y)),v); n+=1
    # pi
    pi=bc(["4*a(1)"])[0]; bad+= not inside(k.PI,pi); n+=1
    print(f"{tag}: {n} bc-oracle comparisons, enclosure failures {bad}")
    return bad
getcontext().prec=130
k=load(src); run(k,"REAL kernel (vs bc, 120-digit independent oracle)")
# mutants: weakened Taylor (K small, remainder dropped)
m1=src.replace("e=absmax(term); return (s[0]-e,s[1]+e)","return s").replace('KT=int(os.environ.get("IVK_K","24"))','KT=5')
run(load(m1),"MUTANT (K=5, no Taylor remainder)")
m2=src.replace("e=absmax(term); s=(s[0]-e,s[1]+e)","pass").replace('EXPN=int(os.environ.get("IVK_N","24"))','EXPN=12')
run(load(m2),"MUTANT (exp: n=12, no remainder)")
