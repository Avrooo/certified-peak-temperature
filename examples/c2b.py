import e4,e3,solver
a=0.2/(1200*1800); Lm=0.004
def run(h1,h2,ts,tp,lim,T0=25,Tinf=100,unc=0.2,show=True):
    B1,B2=h1*Lm/0.2,h2*Lm/0.2; Fs,Fp=a*ts/Lm**2,a*tp/Lm**2; L=(lim-T0)/(Tinf-T0)
    box=((B1*(1-unc),B1*(1+unc)),(B2*(1-unc),B2*(1+unc)))
    out={}
    for nm,(b1,b2) in (("lo",(box[0][0],box[1][0])),("hi",(box[0][1],box[1][1]))):
        s,p=e4.endpoints_float(b1,b2,Fs,Fp); g20=e4.full_history_peak(b1,b2,Fs,Fp,300,20)
        out[nm]=(s,p,g20)
    d=solver.decide(box,Fs,Fp,L)
    print(f"h=({h1},{h2}) ts={ts} tp={tp} B=({B1:.3g},{B2:.3g}) Fs={Fs:.4g} Fp={Fp:.4g}")
    for k,v in out.items(): print(" ",k,"theta_s(ts)=%.6f theta_s(tp)=%.6f grid20=%.6f  -> degC switch %.2f end %.2f"%(v[0],v[1],v[2],T0+75*v[0],T0+75*v[1]))
    print("  V7",d[0],d[1]['width'],"L",L)
run(500,5,10,180,60)
run(400,5,15,180,60)
