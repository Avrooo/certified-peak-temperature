import numpy as np, time
from math import sqrt,pi
from scipy.special import erfcx as sp_erfcx
from adaptive import root
from fast import dg
def gfun(s):
    return 0.0 if s<=0 else 2*sqrt(s/pi)+dg(s)
# ---------- exact reference: spectral propagator across the Bi switch ----------
def eig(Bi,K):
    z=np.array([root(k,Bi) for k in range(K)]); c=4*np.sin(z)/(2*z+np.sin(2*z)); return z,c
def theta_s_two(Bi1,Bi2,Fs,F,K=80):
    """surface temperature at F (>=Fs) for Bi=Bi1 on [0,Fs), Bi2 afterwards; ambient 1, initial 0"""
    z1,c1=eig(Bi1,K); u=c1*np.exp(-z1**2*Fs)                    # u=1-theta coefficients in basis cos(z1 x)
    if F<=Fs: return 1-np.sum(c1*np.cos(z1)*np.exp(-z1**2*F))
    z2,_=eig(Bi2,K)
    Z1,Z2=np.meshgrid(z1,z2)
    I=0.5*(np.sinc((Z2-Z1)/pi)+np.sinc((Z2+Z1)/pi))              # int_0^1 cos(z2 x)cos(z1 x) dx
    Nm=0.5*(1+np.sin(2*z2)/(2*z2))
    cp=(I@u)/Nm
    return 1-np.sum(cp*np.cos(z2)*np.exp(-z2**2*(F-Fs)))
def theta_s_path(Bi1,Bi2,Fs,Fp,n=400,K=80):
    ts=np.linspace(1e-4,Fp,n); return ts,np.array([theta_s_two(Bi1,Bi2,Fs,t,K) for t in ts])
# ---------- certified bracket (Bi1<=Bi2): comparison envelopes + step-flux Volterra bound ----------
def Llow(Bi,t): return 1-np.exp((Bi*sqrt(t))**2)*__import__('math').erfc(Bi*sqrt(t)) if t>0 else 0.0
def Uup(Bi,t):
    return 0.0 if t<=0 else min(1.0,Llow(Bi,t)+Bi*dg(t))
def bracket(Bi1,Bi2,Fs,Fp,N1,N2):
    ts=list(np.linspace(0,Fs,N1+1))+list(np.linspace(Fs,Fp,N2+1)[1:]); Bs=[Bi1]*N1+[Bi2]*N2
    gg=[gfun(Fp-t) for t in ts]; up=lo=0.0
    for j in range(len(Bs)):
        w=gg[j]-gg[j+1]
        up+=Bs[j]*(1-Llow(Bi1,ts[j]))*w                               # theta(t)>=S_{Bi1}(t)>=Llow
        lo+=Bs[j]*(1-Uup(Bi2,ts[j+1]))*w                              # theta(t)<=S_{Bi2}(t)<=Uup
    return max(lo,0.0),min(up,1.0)
# ---------- fast Volterra product-integration solver (the 'fast solver' competitor) ----------
def volterra(Bi1,Bi2,Fs,Fp,N):
    h=Fp/N; lag=np.array([gfun(k*h) for k in range(N+1)]); w=lag[1:]-lag[:-1]   # w[m]=int over one cell at lag m..m+1
    # kernel weight for flux in cell i at time t_n: g(t_n-t_{i-1})-g(t_n-t_i)=lag[n-i+1]-lag[n-i]=w[n-i]
    th=np.zeros(N+1); q=np.zeros(N+1)
    for n in range(1,N+1):
        Bn=Bi1 if n*h<=Fs+1e-15 else Bi2
        A=np.dot(q[1:n],w[n-np.arange(1,n)]) if n>1 else 0.0
        th[n]=(A+Bn*w[0])/(1+Bn*w[0]); q[n]=Bn*(1-th[n])
    return th[-1]
if __name__=="__main__":
    # validity of reference: constant Bi and K convergence
    from pilot_lib import roots, step
    z,C=roots(0.1,300); print("ref sanity (Bi1=Bi2=0.1,F=0.5): spectral %.10f vs eigen %.10f"%(theta_s_two(.1,.1,.2,.5),step(.1,.5,z,C,1.0)[0]))
    print("K convergence (0.03->0.3, Fs=.1, Fp=.5): K=40 %.10f  K=80 %.10f  K=160 %.10f"%tuple(theta_s_two(.03,.3,.1,.5,K) for K in (40,80,160)))
    print("Volterra product integration error (same case, exact %.8f):"%theta_s_two(.03,.3,.1,.5), [ "N=%d err %.1e"%(N,abs(volterra(.03,.3,.1,.5,N)-theta_s_two(.03,.3,.1,.5))) for N in (64,256,1024)])
