import math
from decimal import Decimal as D, getcontext
from scipy.optimize import brentq
import e3
getcontext().prec=60
EPS=D(10)**-70
def _atan_inv(x):
    s=D(0); t=D(1)/x; x2=x*x; k=0
    while abs(t)>EPS: s+=t/(2*k+1)*(1 if k%2==0 else -1) if False else (t/(2*k+1) if k%2==0 else -t/(2*k+1)); t=t/x2; k+=1
    return s
PI=16*_atan_inv(D(5))-4*_atan_inv(D(239))
def sincos(x):
    x2=x*x; s=D(0); c=D(0); ts=x; tc=D(1); k=0
    while abs(ts)>EPS or abs(tc)>EPS:
        s+=ts; c+=tc
        ts=-ts*x2/((2*k+2)*(2*k+3)); tc=-tc*x2/((2*k+1)*(2*k+2)); k+=1
    return s,c
_rc={}
def root(m,B):
    key=(m,B)
    if key in _rc: return _rc[key]
    if B==0: _rc[key]=D(0); return D(0)
    seed=brentq(lambda d:(m*math.pi+d)*math.tan(d)-float(B),1e-300,math.pi/2-1e-13,xtol=1e-17,rtol=8.9e-16)
    d=D(seed); Bd=D(B); mp=m*PI
    for _ in range(7):
        s,c=sincos(d); g=(mp+d)*s-Bd*c; gp=(1+Bd)*s+(mp+d)*c; d=d-g/gp
    s,c=sincos(d); assert abs((mp+d)*s-Bd*c)<D(10)**-45
    _rc[key]=d; return d
def mode(B,m):
    d=root(m,B); z=m*PI+d; s,c=sincos(d); sg=1 if m%2==0 else -1
    if z==0: return dict(d=d,z=z,cz=D(1),N=D(1),c=D(1))
    N=(1+s*c/z)/2; return dict(d=d,z=z,cz=sg*c,N=N,c=(sg*s/z)/N)
def sinc(x):
    if x==0: return D(1)
    s,_=sincos(abs(x)); return s/abs(x)
def theta_endpoints(B1,B2,Fs,Fp,tol=1e-32):
    B1d,B2d,Fsd,Fpd=D(B1),D(B2),D(Fs),D(Fp); tau=Fpd-Fsd
    M=e3.pick_M(lambda M:e3.EM_float(B1,B2,Fs,Fp,M),tol,Mmax=900)+2
    M1=e3.pick_M(lambda M:e3.EM_float(B1,0,Fs,Fs,M,one=True),tol,Mmax=900)+2
    m1=[mode(B1d,m) for m in range(max(M,M1))]; m2=[mode(B2d,n) for n in range(M)]
    u1=sum(q['c']*q['cz']*(-(q['z']**2)*Fsd).exp() for q in m1[:M1])
    a=[q['c']*(-(q['z']**2)*Fsd).exp() for q in m1[:M]]
    dn=[r['cz']/r['N']*(-(r['z']**2)*tau).exp() for r in m2]
    dB=B1d-B2d; u2=D(0)
    for i in range(M):
        q=m1[i]; row=D(0)
        for j in range(M):
            r=m2[j]
            if i!=j: I=dB*q['cz']*r['cz']/((q['z']-r['z'])*(q['z']+r['z']))
            else:
                zp=q['z']+r['z']; sp=sinc(zp) if i==0 else (sincos(q['d']+r['d'])[0]/zp)
                I=(sinc(q['d']-r['d'])+sp)/2
            row+=I*dn[j]
        u2+=a[i]*row
    return 1-u1,1-u2            # theta_s(Fs), theta_s(Fp)
def peak(B1,B2,Fs,Fp): a,b=theta_endpoints(B1,B2,Fs,Fp); return max(a,b)
