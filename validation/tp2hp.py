import numpy as np, math
from decimal import Decimal as D, getcontext
import e4ref as R
getcontext().prec=70
TWO=2*R.PI
def cosr(a):
    k=(a/TWO).to_integral_value(); r=a-k*TWO; return R.sincos(r)[1]
def Gh(B,x,y,t,K):
    s=D(0); Bd=D(B); td=D(t)
    for m in range(K):
        md=R.mode(B,m); z=md['z']; 
        s+=cosr(z*D(x))*cosr(z*D(y))*(-(z*z)*td).exp()/md['N']
    return s
cases=[(82.3,0.00373,(0.6499,0.809),(0.0107,0.2137)),(0,0.00323,(0.8104,0.9223),(0.1422,0.7235)),(0.0719,0.00361,(0.8258,0.8681),(0.0029,0.3836)),(4.62,0.00359,(0.117,0.9146),(0.0153,0.1518)),(25.3,0.00357,(0.1408,0.1588),(0.259,0.9614)),(2.08,0.00336,(0.7214,0.7967),(0.0396,0.3005))]
for B,t,x,y,*_ in [(c[0],c[1],c[2],c[3]) for c in cases]:
    K=260
    g=[[Gh(B,xx,yy,t,K) for yy in y] for xx in x]
    Dd=g[0][0]*g[1][1]-g[0][1]*g[1][0]
    print("B=%g t=%g  G11=%.3e G22=%.3e  D=%.6e  D>=0: %s"%(B,t,g[0][0],g[1][1],Dd,Dd>=0))
