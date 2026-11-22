import io,contextlib,math,time,json,sys,numpy as np
from decimal import Decimal as D
from scipy.linalg import solve_banded
import e3,ivk as k,e4ref as R
with contextlib.redirect_stdout(io.StringIO()): import closed
SCHED=[1e-3,1e-5,1e-7,1e-9,1e-11,1e-13]
def clr(): e3._rootcache.clear(); e3._modecache.clear()
# ---------------- float building blocks ----------------
def Kfor(B1,B2,Fs,Fp,tol=1e-14,one=False):
    return min(600,e3.pick_M(lambda M:e3.EM_float(B1,B2,Fs,Fp,M,one=one),tol)+2)
def single_float(B1,Fs):
    K=Kfor(B1,0,Fs,Fs,one=True); z=closed.eigz(B1,K); c=np.sinc(z/math.pi)/closed.Nn(z)
    return 1-float((c*np.cos(z)*np.exp(-z**2*Fs)).sum())
def two_float(B1,B2,Fs,Fp):
    K=Kfor(B1,B2,Fs,Fp); return 1-float(closed.double(B1,B2,Fs,Fp,K).sum())
def endpoints_float(B1,B2,Fs,Fp): return single_float(B1,Fs),two_float(B1,B2,Fs,Fp)
# B2: full history on uniform time grid (no knowledge of Fs)
def full_history_peak(B1,B2,Fs,Fp,K=300,nt=200):
    z=closed.eigz(B1,K); w=closed.eigz(B2,K); N1=closed.Nn(z); N2=closed.Nn(w); c=np.sinc(z/math.pi)/N1
    Z,W=np.meshgrid(z,w,indexing='ij'); d=Z**2-W**2
    with np.errstate(all='ignore'): I=(B1-B2)*np.cos(Z)*np.cos(W)/d
    diag=0.5*(np.sinc((Z-W)/math.pi)+np.sinc((Z+W)/math.pi)); I=np.where(np.abs(Z-W)<1,diag,I)
    v=((c*np.exp(-z**2*Fs))[:,None]*I).sum(0)*np.cos(w)/N2
    best=-1
    for i in range(1,nt+1):
        t=Fp*i/nt
        if t<=Fs: u=(c*np.cos(z)*np.exp(-z**2*t)).sum()
        else: u=(v*np.exp(-w**2*(t-Fs))).sum()
        best=max(best,1-float(u))
    return best
# B1: finite difference (backward Euler), deficit u=1-theta, node-centred with Robin ghost node
def fd_peak(B1,B2,Fs,Fp,N=100,steps=1200):
    h=1.0/N; u=np.ones(N+1); best=0.0
    def band(B,dt):
        ab=np.zeros((3,N+1)); r=dt/h**2
        ab[1,:]=1+2*r; ab[0,1:]=-r; ab[2,:-1]=-r; ab[0,1]=-2*r         # x=0 Neumann: u0_t=(2/h^2)(u1-u0)
        ab[2,N-1]=-2*r; ab[1,N]=1+2*r+2*B*dt/h                           # x=1: u_t=(2/h^2)(uN-1 - uN)-(2B/h)uN
        return ab
    n1=steps//2; n2=steps-n1
    for B,T,n in ((B1,Fs,n1),(B2,Fp-Fs,n2)):
        dt=T/n; ab=band(B,dt)
        for _ in range(n):
            u=solve_banded((1,1),ab,u); best=max(best,1-u[-1])
    return best
# ---------------- decision rule from approximate corner peaks ----------------
def decide_float(pU,pL,L):
    if pU<L: return 'SAFE_ALL'
    if pL>=L: return 'EXCEEDS_ALL'
    return 'MIXED'
# ---------------- interval methods ----------------
def ivl_peak_status(pk,L):
    Lq=k.fl(L)
    if pk[1]<Lq[0]: return 'below'
    if pk[0]>=Lq[1]: return 'above'
    return None
def combine(Us,Ls):
    if Us=='below': return 'SAFE_ALL'
    if Ls=='above': return 'EXCEEDS_ALL'
    if Us=='above' and Ls=='below': return 'MIXED'
    return 'INCONCLUSIVE'
def b4_fixed(box,Fs,Fp,L,M=64):
    (b1l,b1h),(b2l,b2h)=box; w=[]
    def corner(B1,B2):
        try:
            u1,_=e3.u_single(B1,Fs,0,fixedM=M); u2,_=e3.u_two(B1,B2,Fs,Fp,0,fixedM=M)
        except e3.CertFail: return None,None
        ths=k.sub(k.ONE,u1); thp=k.sub(k.ONE,u2); pk=(max(ths[0],thp[0]),max(ths[1],thp[1]))
        w.append((pk[1]-pk[0])/k.S); return ivl_peak_status(pk,L),pk
    Us,_=corner(b1h,b2h)
    if Us=='below': return 'SAFE_ALL',max(w) if w else None
    Ls,_=corner(b1l,b2l); return combine(Us,Ls),max(w) if w else None
def b5_adaptive(box,Fs,Fp,L):
    (b1l,b1h),(b2l,b2h)=box; w=[]
    def corner(B1,B2):
        for eps in SCHED:
            try: pk,_,_,M=e3.theta_peak(B1,B2,Fs,Fp,eps)
            except e3.CertFail: return None
            w.append((pk[1]-pk[0])/k.S); st=ivl_peak_status(pk,L)
            if st: return st
        return None
    Us=corner(b1h,b2h)
    if Us=='below': return 'SAFE_ALL',w[-1]
    Ls=corner(b1l,b2l); return combine(Us,Ls),(max(w) if w else None)
def b6_decision_first(box,Fs,Fp,L):
    (b1l,b1h),(b2l,b2h)=box; w=[]
    def endpoint(which,B1,B2,approx):
        eps_list=[e for e in SCHED if e<=max(1e-13,min(1e-3,abs(approx-L)/4))] or [SCHED[-1]]
        eps_list=sorted(set([max(1e-13,min(1e-3,abs(approx-L)/4))]+[e for e in SCHED if e<max(1e-13,min(1e-3,abs(approx-L)/4))]),reverse=True)
        for eps in eps_list:
            try:
                u,_=(e3.u_single(B1,Fs,eps) if which=='s' else e3.u_two(B1,B2,Fs,Fp,eps))
            except e3.CertFail: return None
            th=k.sub(k.ONE,u); w.append((th[1]-th[0])/k.S); Lq=k.fl(L)
            if th[1]<Lq[0]: return 'below'
            if th[0]>=Lq[1]: return 'above'
        return None
    def corner(B1,B2,ap):
        order=sorted(('s','p'),key=lambda q:-ap[q]); res=[]
        for q in order:
            r=endpoint(q,B1,B2,ap[q]); res.append(r)
            if r=='above': return 'above'
            if r=='below' and len(res)==2: return 'below' if res[0]=='below' else None
        return None
    sU=endpoints_float(b1h,b2h,Fs,Fp); sL=endpoints_float(b1l,b2l,Fs,Fp)
    apU=dict(s=sU[0],p=sU[1]); apL=dict(s=sL[0],p=sL[1]); pU=max(sU); pL=max(sL)
    if pU<L:   # likely safe: certify upper first
        Us=corner(b1h,b2h,apU)
        if Us=='below': return 'SAFE_ALL',max(w)
        Ls=corner(b1l,b2l,apL); return combine(Us,Ls),max(w)
    if pL>=L:  # likely exceeds: certify lower first
        Ls=corner(b1l,b2l,apL)
        if Ls=='above': return 'EXCEEDS_ALL',max(w)
        Us=corner(b1h,b2h,apU); return combine(Us,Ls),max(w)
    Us=corner(b1h,b2h,apU); Ls=corner(b1l,b2l,apL); return combine(Us,Ls),(max(w) if w else None)
# ---------------- case generation ----------------
def gen_box(group,rng):
    if group=='G2':
        kind=rng.choice(['short','bigB','near','zero'])
        Fs=10**rng.uniform(-3.5,-2.7) if kind=='short' else 10**rng.uniform(-2.5,-0.5)
        b1=10**rng.uniform(0,1.2) if kind!='bigB' else 10**rng.uniform(2,2.7); b2=10**rng.uniform(0,1.2) if kind!='bigB' else 10**rng.uniform(1.3,2.7)
        w1=1+10**rng.uniform(-2,-.5); w2=w1
        if kind=='near': b2=b1*(1+10**rng.uniform(-9,-5)); w1=w2=1+1e-7
        box=((b1,b1*w1),(b2,b2*w2))
        if kind=='zero':
            box=((b1,b1*w1),(0.0,0.0)) if rng.random()<.5 else ((0.0,0.0),(0.5,0.5*w1))
    elif group=='G4':
        kind='indep'; Fs=10**rng.uniform(-3,0.0); b1=10**rng.uniform(-0.7,2.3); b2=10**rng.uniform(-0.7,2.3)
        box=((b1,b1*(1+rng.uniform(.02,.8))),(b2,b2*(1+rng.uniform(.02,.8))))
    else:
        kind='ord'; Fs=10**rng.uniform(-3,-0.4); b1=10**rng.uniform(-1.3,1.7); b2=10**rng.uniform(-1.3,1.7)
        box=((b1,b1*(1+rng.uniform(.05,.5))),(b2,b2*(1+rng.uniform(.05,.5))))
    ratio=10**rng.uniform(0.05,1.3) if group!='G4' else 10**rng.uniform(0.05,2.0)
    return kind,box,Fs,Fs*ratio
def make_case(group,seed):
    rng=np.random.default_rng(seed); 
    while True:
        kind,box,Fs,Fp=gen_box(group,rng)
        (b1l,b1h),(b2l,b2h)=box
        pU=R.peak(b1h,b2h,Fs,Fp); pL=R.peak(b1l,b2l,Fs,Fp)
        u=rng.random(); cls='safe' if u<.25 else 'exceed' if u<.5 else 'mixed' if u<.8 else 'neareq'
        if group=='G3': s=D(10)**D(float(rng.uniform(-13,-8)))
        elif group=='G4': s=D(10)**D(float(rng.uniform(-10,-2)))
        else: s=D(float(rng.choice([1e-2,1e-3,1e-4,1e-6])))
        if cls=='safe': Ld=pU+s
        elif cls=='exceed': Ld=pL-s
        elif cls=='mixed':
            if pU-pL<=2*s: continue
            Ld=pL+s
        else:
            e=D(float(rng.uniform(-3e-15,3e-15))); Ld=(pU if rng.random()<.5 else pL)+e
        L=float(Ld); Lx=D(L)
        truth='SAFE_ALL' if pU<Lx else 'EXCEEDS_ALL' if pL>=Lx else 'MIXED'
        sep=float(min(abs(Lx-pU),abs(Lx-pL)))
        return dict(seed=seed,group=group,kind=kind,cls=cls,box=box,Fs=Fs,Fp=Fp,L=L,truth=truth,sep=sep,pU=float(pU),pL=float(pL))
def run_case(args):
    group,seed=args; c=make_case(group,seed); box=c['box']; (b1l,b1h),(b2l,b2h)=box; Fs,Fp,L=c['Fs'],c['Fp'],c['L']; out=dict(c)
    def timed(name,f):
        clr(); t0=time.perf_counter()
        try: r=f()
        except Exception as ex: r=('ERROR:'+repr(ex)[:60],None)
        out[name]=dict(dec=r[0],w=r[1],t=time.perf_counter()-t0)
    timed('B1',lambda:(decide_float(fd_peak(b1h,b2h,Fs,Fp),fd_peak(b1l,b2l,Fs,Fp),L),None))
    timed('B2',lambda:(decide_float(full_history_peak(b1h,b2h,Fs,Fp),full_history_peak(b1l,b2l,Fs,Fp),L),None))
    timed('B3',lambda:(decide_float(max(endpoints_float(b1h,b2h,Fs,Fp)),max(endpoints_float(b1l,b2l,Fs,Fp)),L),None))
    timed('B4',lambda:b4_fixed(box,Fs,Fp,L)); timed('B5',lambda:b5_adaptive(box,Fs,Fp,L)); timed('B6',lambda:b6_decision_first(box,Fs,Fp,L))
    return out
if __name__=="__main__":
    group=sys.argv[1]; n=int(sys.argv[2]); seed0=int(sys.argv[3]); tag=sys.argv[4]
    with open(f"e4_{tag}.jsonl","a") as f:
        for i in range(n):
            r=run_case((group,seed0+i)); f.write(json.dumps(r,default=float)+"\n"); f.flush()
