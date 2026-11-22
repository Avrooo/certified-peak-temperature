import numpy as np, io, contextlib, time
with contextlib.redirect_stdout(io.StringIO()):
    from twostage import theta_s_two as th
rng=np.random.default_rng(2026)
def pk(B1,B2,Fs,Fp,K=80): return max(th(B1,B2,Fs,Fs,K),th(B1,B2,Fs,Fp,K))
# ---------- E1: broad conditions, end-only vs theorem-aware vs dense-history
N=400; rows=[]
for i in range(N):
    B1=10**rng.uniform(-2,1.3); B2=0.0 if i%10==0 else 10**rng.uniform(-2,1.3)
    Fs=10**rng.uniform(-3,-0.3); Fp=Fs*10**rng.uniform(0.02,1.3)
    end=th(B1,B2,Fs,Fp) if B2>0 else th(B1,1e-12,Fs,Fp); sw=th(B1,max(B2,1e-12),Fs,Fs)
    P=max(end,sw)
    ts=np.linspace(Fs*0.5,Fp,150); dense=max(th(B1,max(B2,1e-12),Fs,t) for t in ts); dense=max(dense,sw)
    rows.append((B1,B2,end,sw,P,dense))
R=np.array(rows)
miss=R[:,4]-R[:,2]
print("E1  cases:",N," B1>B2:",int((R[:,0]>R[:,1]).sum()))
for thr in (1e-2,1e-4,1e-6): print(f"   end-only underestimates peak by >{thr:g}: {np.mean(miss>thr)*100:.1f}%")
print("   max theorem-aware minus dense-history (should be >=0 or ~0):",float(np.max(R[:,5]-R[:,4])))
# ---------- E2: uncertain convection, box +-20%
def box(B1,B2,r=0.2): return B1*(1-r),B1*(1+r),B2*(1-r),B2*(1+r)
for delta in (1e-2,1e-4,1e-6):
    fs={"M1 nominal end-only":0,"M2 nominal peak":0,"M3 20 random samples":0,"M3b 200 random samples":0,"M4 corner rule":0}; n=0
    for i in range(300):
        B1=10**rng.uniform(-1.5,1.0); B2=10**rng.uniform(-2,0.7)
        Fs=10**rng.uniform(-2.5,-0.5); Fp=Fs*10**rng.uniform(0.05,1.0)
        l1,h1,l2,h2=box(B1,B2)
        Tup=pk(h1,h2,Fs,Fp)
        L=Tup-delta          # truth: some admissible pair exceeds L, so "safe for all" is FALSE
        if L<=0: continue
        n+=1
        if th(B1,B2,Fs,Fp)<L: fs["M1 nominal end-only"]+=1
        if pk(B1,B2,Fs,Fp)<L: fs["M2 nominal peak"]+=1
        for name,m in (("M3 20 random samples",20),("M3b 200 random samples",200)):
            ok=all(pk(rng.uniform(l1,h1),rng.uniform(l2,h2),Fs,Fp)<L for _ in range(m))
            if ok: fs[name]+=1
        if Tup<L: fs["M4 corner rule"]+=1
    print(f"E2  margin {delta:g}: false-'safe-for-all' declarations out of {n}:",{k:v for k,v in fs.items()})
# ---------- timing (scalar python, K=80)
args=[(10**rng.uniform(-1,1),10**rng.uniform(-2,0.5),10**rng.uniform(-2.5,-0.5)) for _ in range(40)]
t0=time.perf_counter()
for B1,B2,Fs in args: pk(B1*1.2,B2*1.2,Fs,Fs*5)
tc=(time.perf_counter()-t0)/len(args)
t0=time.perf_counter()
for B1,B2,Fs in args:
    for t in np.linspace(Fs,Fs*5,200): th(B1*1.2,B2*1.2,Fs,t)
th_=(time.perf_counter()-t0)/len(args)
print(f"timing per box decision: corner rule (2 corners x 2 times)= {2*tc*1e3:.2f} ms ; one corner with full 200-point history = {th_*1e3:.1f} ms")
