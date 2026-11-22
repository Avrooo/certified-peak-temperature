import numpy as np, io, contextlib
from math import pi
with contextlib.redirect_stdout(io.StringIO()):
    from twostage import eig, theta_s_two
def eigz(B,K):
    if B==0.0:
        z=np.arange(K)*pi; return z   # exact: zeta_n=n*pi, zeta_0=0
    return eig(B,K)[0]
def Nn(z): return 0.5*(1+np.sinc(2*z/pi))   # sinc form: N(0)=1 exactly
def double(B1,B2,Fs,Fp,K):
    """exact (up to K) double series for u(1,t_p) with CLOSED-FORM overlap; returns terms matrix"""
    z=eigz(B1,K); w=eigz(B2,K); tau=Fp-Fs
    N1=Nn(z); N2=Nn(w)
    c=np.sinc(z/pi)/N1                     # coefficient of the constant 1 in stage-1 basis
    Z,W=np.meshgrid(z,w,indexing='ij')    # rows m (stage1), cols n (stage2)
    d=Z**2-W**2
    with np.errstate(all='ignore'):
        I=(B1-B2)*np.cos(Z)*np.cos(W)/d   # closed form for z_m != zeta_n
    diag=0.5*(np.sinc((Z-W)/pi)+np.sinc((Z+W)/pi))
    near=np.abs(Z-W)<1.0                  # same-index pairs: use stable sinc form
    I=np.where(near,diag,I)
    T=(c*np.exp(-z**2*Fs))[:,None]*I*(np.cos(w)/N2*np.exp(-w**2*tau))[None,:]
    return T
def tail_bound(B1,B2,Fs,Fp,M):
    """analytic bound on |sum of all terms with m>=M or n>=M| (indices 0..M-1 kept)"""
    tau=Fp-Fs; dB=abs(B1-B2); Nmin=0.39
    def S(t,M):  # sum_{k>=M} exp(-(k pi)^2 t)
        a=np.exp(-(M*pi)**2*t); r=np.exp(-(2*M+1)*pi**2*t); return a/(1-r) if r<1 else np.inf
    # |c_m e^{-z^2 Fs}| <= cm_bound ; m>=1: B1/(Nmin z^2); m=0: 1/Nmin
    # |cos(w_n)/N2 e^{-w^2 tau}| <= 1/Nmin e^{-(n pi)^2 tau}
    # off-diagonal: 1/|z^2-w^2| <= 2/(pi*max(z,w)); diagonal |I|<=1
    # row-wise bound, region A: m>=M (all n): sum_m |c_m|e^{-z_m^2 Fs} * sum_n (1/Nmin)e^{-(n pi)^2 tau} * min(1, dB*2/(pi*z_m)) 
    # (crude: use the larger of the two cases via min(1,.) on every pair, since |I|<=1 always)
    Sm_all_tau=1+S(tau,1)                           # sum_{n>=0} e^{-(n pi)^2 tau}
    amax_m=lambda m0: (B1/(Nmin*(m0*pi)**2))        # |c_m| for m>=m0>=1
    # region m>=M:
    f=min(1.0,dB*2/(pi*M*pi))                      # valid for OFF-diagonal pairs only (|z-w|>=pi/2)
    tA=(1/Nmin)*Sm_all_tau*amax_m(M)*f*S(Fs,M)
    tA+=B1/(Nmin**2*(M*pi)**2)*S(Fp,M)             # diagonal m=n>=M: |I|<=1, no smallness in dB
    # region n>=M, m<M: sum_m |c_m| e^{-z_m^2 Fs}*(1/Nmin)*S_n * min(1,2 dB/(pi*M*pi))
    cm_sum=1/Nmin+B1/(Nmin*6.0)                    # sum_m |c_m| <= 1/Nmin + (B1/Nmin) sum 1/(m pi)^2 = ... ; sum_{m>=1} 1/(m pi)^2 = 1/6
    tB=(1/Nmin)*S(tau,M)*cm_sum*min(1.0,dB*2/(pi*M*pi))
    return tA+tB
def adversarial():
    rng=np.random.default_rng(5); bad=0; n=0; worst=0
    for B1 in (0.05,1.0,8.0,40.0,200.0):
      for rel in (0.0,1e-9,1e-6,1e-3,0.05,0.5,-1.0):
        B2=0.0 if rel==-1.0 else B1*(1+rel)
        for Fs in (3e-4,3e-3,3e-2,0.3):
          for r in (1.01,1.5,5,50):
            Fp=Fs*r; Kr=600; T=double(B1,B2,Fs,Fp,Kr); uref=T.sum()
            for tol in (1e-2,1e-4,1e-6):
              M=next((M for M in range(2,Kr//2) if tail_bound(B1,B2,Fs,Fp,M)<tol),None)
              n+=1
              if M is None: continue
              act=abs(T[:M,:M].sum()-uref); b=tail_bound(B1,B2,Fs,Fp,M)
              if act>b+1e-13: bad+=1
              worst=max(worst,M)
    print('adversarial: checks',n,'violations',bad,'max M',worst)
if __name__=="__main__":
    adversarial()
    rng=np.random.default_rng(11)
    bad=0; rows=[]
    for i in range(60):
        B1=10**rng.uniform(-1.5,1.2); B2=0.0 if i%8==0 else 10**rng.uniform(-1.5,1.2)
        Fs=10**rng.uniform(-2.2,-0.4); Fp=Fs*10**rng.uniform(0.05,1.2)
        ref_th=th=theta_s_two(B1,max(B2,1e-12),Fs,Fp,K=120)
        Kref=300; T=double(B1,B2,Fs,Fp,Kref); uref=T.sum()
        cond=np.abs(T).sum()/max(abs(uref),1e-300)
        err0=abs((1-uref)-th)
        out={}
        for tol in (1e-2,1e-4,1e-6):
            M=next((M for M in range(2,Kref) if tail_bound(B1,B2,Fs,Fp,M)<tol),None)
            if M is None: out[tol]=None; continue
            uM=T[:M,:M].sum(); actual=abs(uM-uref); out[tol]=(M,actual,actual<=tail_bound(B1,B2,Fs,Fp,M)+1e-12)
            if not out[tol][2]: bad+=1
        rows.append((B1,B2,Fs,Fp,err0,cond,out))
    print("closed-form double series vs theta_s_two(K=120): max abs diff",max(r[4] for r in rows))
    print("cancellation ratio (sum|T|/|sum T|): median %.2f max %.2f"%(np.median([r[5] for r in rows]),max(r[5] for r in rows)))
    print("bound violations (actual truncation error > bound):",bad)
    for tol in (1e-2,1e-4,1e-6):
        Ms=[r[6][tol][0] for r in rows if r[6][tol]]; fail=sum(1 for r in rows if r[6][tol] is None)
        print(f"tol {tol:g}: modes M needed (m,n each): median {np.median(Ms):.0f}, max {max(Ms)}, no M<{300} found: {fail}")
