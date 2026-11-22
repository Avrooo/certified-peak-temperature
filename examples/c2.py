import json,math,numpy as np
import solver,e4,e3
import ivk
a=0.2/(1200*1800); Lm=0.004
Fs=a*30/Lm**2; Fp=a*180/Lm**2
B1n,B2n=60*Lm/0.2,15*Lm/0.2
L=(65-25)/(100-25)
box=((B1n*0.8,B1n*1.2),(B2n*0.8,B2n*1.2))
print("Fs,Fp,L",Fs,Fp,L,box)
# certified decision
dec=solver.decide(box,Fs,Fp,L)
print("V7 decision",dec)
# float corner peaks, ends, time grid
def surf(B1,B2,t):
    return None
def peaks(B1,B2):
    ths,thp=e4.endpoints_float(B1,B2,Fs,Fp)
    h20=e4.full_history_peak(B1,B2,Fs,Fp,K=300,nt=20); h200=e4.full_history_peak(B1,B2,Fs,Fp,K=300,nt=200)
    return dict(th_s=ths,th_p=thp,end=thp,peak=max(ths,thp),grid20=h20,grid200=h200)
lo=peaks(*(box[0][0],box[1][0])); hi=peaks(box[0][1],box[1][1])
print("lower",lo); print("upper",hi)
# certified enclosures at both corners
for nm,(b1,b2) in (("lower",(box[0][0],box[1][0])),("upper",(box[0][1],box[1][1]))):
    (l,h),ths,thp,M=e3.theta_peak(b1,b2,Fs,Fp,1e-10); S=2.0**96
    print(nm,"cert peak",float(l)/S,float(h)/S,"theta_s(Fs)",[float(x)/S for x in ths],"theta_p",[float(x)/S for x in thp],"M",M)
