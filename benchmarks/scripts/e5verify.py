import json,glob,random,time
from fractions import Fraction as F
from decimal import Decimal as D
import stacks,e4ref as R
FR=stacks.load_stack("frozen","frozen",256,50,52,0,0,0); V7=stacks.load_stack("v7",".",96,24,24,1,1,1)
C=[json.loads(l) for f in sorted(glob.glob('e4_G*.jsonl')) for l in open(f)]
random.seed(3); near=[r for r in C if r['sep']<1e-14]; oth=random.sample([r for r in C if r['sep']>=1e-14],200); S=near+oth
n=0; v7_miss=fr_miss=nonov=0; wr=[]
for r in S:
    (b1l,b1h),(b2l,b2h)=r['box']
    for B1,B2 in ((b1h,b2h),(b1l,b2l)):
        ref=F(R.peak(B1,B2,r['Fs'],r['Fp'])); ivs=[]
        for st in (FR,V7):
            stacks.clear(st); pk,_,_,M=st.e3.theta_peak(B1,B2,r['Fs'],r['Fp'],1e-13); ivs.append((F(pk[0],st.ivk.S),F(pk[1],st.ivk.S)))
        (fl_,fh),(vl,vh)=ivs; n+=1
        fr_miss+= not(fl_<=ref<=fh); v7_miss+= not(vl<=ref<=vh); nonov+= max(fl_,vl)>min(fh,vh)
        wr.append(float((vh-vl)/(fh-fl_)))
wr.sort()
print(f"corner evaluations vs 60-digit reference: {n} (334 near-ties + 200 others, eps=1e-13)")
print("reference outside interval: frozen",fr_miss,"| V7",v7_miss,"| frozen/V7 intervals disjoint:",nonov)
print("width ratio V7/frozen: min %.3f median %.3f max %.3f"%(wr[0],wr[len(wr)//2],wr[-1]))
