import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import e4; from e4 import closed
import os; OUT=os.environ.get('OUTDIR','.')+'/'; DPI=900
plt.rcParams.update({'font.size':10,'axes.labelsize':10,'axes.titlesize':10,'xtick.labelsize':9.5,'ytick.labelsize':9.5,
 'legend.fontsize':9,'axes.spines.top':False,'axes.spines.right':False,'mathtext.default':'regular','font.family':'DejaVu Sans'})
NAVY='#1f4e79'; RED='#b23a2e'; GREY='#6b7078'; GOLD='#d9a441'

# ---------- shared spectral two-stage solver (dimensionless)
def two_stage(B1,B2,Fs,K=300):
    z=closed.eigz(B1,K); w=closed.eigz(B2,K); N1=closed.Nn(z); N2=closed.Nn(w); c=np.sinc(z/np.pi)/N1
    Z,W=np.meshgrid(z,w,indexing='ij'); d=Z**2-W**2
    with np.errstate(all='ignore'): I=(B1-B2)*np.cos(Z)*np.cos(W)/d
    diag=0.5*(np.sinc((Z-W)/np.pi)+np.sinc((Z+W)/np.pi)); I=np.where(np.abs(Z-W)<1,diag,I)
    A=((c*np.exp(-z**2*Fs))[:,None]*I).sum(0)/N2
    return z,w,c,A
def theta(x,t,Fs,z,w,c,A):
    x=np.atleast_1d(x)
    if t<=Fs: u=(c[:,None]*np.cos(np.outer(z,x))*np.exp(-z[:,None]**2*t)).sum(0)
    else: u=(A[:,None]*np.cos(np.outer(w,x))*np.exp(-w[:,None]**2*(t-Fs))).sum(0)
    return 1-u

# ---------- Figure 1
Fs=0.03; ts_=np.linspace(2e-4,0.6,900)
def surf(B1,B2):
    z,w,c,A=two_stage(B1,B2,Fs); return np.array([theta(1.0,t,Fs,z,w,c,A)[0] for t in ts_])
fig,ax=plt.subplots(1,2,figsize=(6.3,2.9))
a=ax[0]; y1=surf(8,0.3); y2=surf(8,0.0)
a.plot(ts_,y1,color=NAVY,lw=1.6,label=r'$B_1=8\to B_2=0.3$'); a.plot(ts_,y2,color=RED,lw=1.6,label=r'$B_1=8\to B_2=0$')
m=ts_>Fs; i=np.argmin(np.where(m,y1,9)); a.plot(ts_[i],y1[i],'o',color=NAVY,ms=5); print('turning point t=',ts_[i])
a.axvline(Fs,color='gray',ls=':',lw=1); a.text(Fs+0.01,0.03,r'switch $t_s$',color='gray',fontsize=9)
a.set_xlabel('Fourier number $t$'); a.set_ylabel(r'surface temperature $\theta_s(t)$'); a.set_title('(a) drop in $B$',fontsize=10)
a.set_ylim(0,0.72); a.legend(loc='upper right',frameon=False,fontsize=8.5)
b=ax[1]; b.plot(ts_,surf(1,4),color=NAVY,lw=1.6,label=r'$B_1=1\to B_2=4$'); b.plot(ts_,surf(2,2),color=GREY,lw=1.6,label=r'$B_1=B_2=2$')
b.axvline(Fs,color='gray',ls=':',lw=1); b.set_xlabel('Fourier number $t$'); b.set_title(r'(b) $B$ not decreasing',fontsize=10)
b.set_ylim(0,0.95); b.legend(loc='lower right',frameon=False,fontsize=8.5)
plt.tight_layout(); plt.savefig(OUT+'fig2.png',dpi=DPI); plt.close()

# ---------- Figure 2 (from stored per-case results)
d=json.load(open('e5final_main.json'))
edges=[1e-2,1e-4,1e-6,1e-8,1e-10,1e-12,1e-14]
def b(s):
    for k,e in enumerate(edges):
        if s>=e: return k
    return len(edges)
n=len(edges)+1; ok=np.zeros(n,int); inc=np.zeros(n,int); bad=np.zeros(n,int)
for r in d:
    k=b(r['sep']); dec=r['v7']['dec']
    if dec=='INCONCLUSIVE': inc[k]+=1
    elif dec==r['truth']: ok[k]+=1
    else: bad[k]+=1
print(ok.sum(),inc.sum(),bad.sum(),ok,inc)
labs=[r'$\geq10^{-2}$',r'$10^{-4}$–$10^{-2}$',r'$10^{-6}$–$10^{-4}$',r'$10^{-8}$–$10^{-6}$',r'$10^{-10}$–$10^{-8}$',r'$10^{-12}$–$10^{-10}$',r'$10^{-14}$–$10^{-12}$',r'$<10^{-14}$']
fig,a=plt.subplots(figsize=(5.6,3.3)); x=np.arange(n)
a.bar(x,ok,color=NAVY,label='reference-consistent'); a.bar(x,inc,bottom=ok,color=GOLD,label='inconclusive'); a.bar(x,bad,bottom=ok+inc,color=RED,label='reference disagreement')
for k in range(n): a.text(k,ok[k]+inc[k]+bad[k]+6,f'{inc[k]}',ha='center',fontsize=9)
a.text(0.01,0.98,'numbers above bars: inconclusive cases',transform=a.transAxes,fontsize=8.5,va='top',color='#444')
a.set_xticks(x); a.set_xticklabels(labs,rotation=35,ha='right',fontsize=9); a.set_ylabel('number of cases')
a.set_xlabel(r'gap between threshold $\theta_{lim}$ and nearest corner peak'); a.set_ylim(0,470)
a.legend(loc='upper center',bbox_to_anchor=(0.55,0.93),frameon=False,fontsize=9)
plt.tight_layout(); plt.savefig(OUT+'fig4.png',dpi=DPI); plt.close()

# ---------- Figures 3 and 4
al=0.2/(1200*1800); Lm=0.004
def hist(B1,B2,Fs,Fp,nt=600):
    z,w,c,A=two_stage(B1,B2,Fs); ts=np.unique(np.concatenate([np.linspace(Fp/nt,Fp,nt),[Fs]]))
    return ts*Lm**2/al, 25+75*np.array([theta(1.0,t,Fs,z,w,c,A)[0] for t in ts])
sc=[(60,15,30,180,65,'A: 60 to 15 W/m²K'),(500,5,10,180,60,'B: 500 to 5 W/m²K')]
fig,ax=plt.subplots(1,2,figsize=(6.3,3.0))
for A_,(h1,h2,ts,tp,lim,ttl) in zip(ax,sc):
    B1,B2=h1*Lm/.2,h2*Lm/.2; Fs,Fp=al*ts/Lm**2,al*tp/Lm**2
    for f,st,lb in ((0.8,'--','lower corner'),(1.2,'-','upper corner')):
        t,T=hist(B1*f,B2*f,Fs,Fp); A_.plot(t,T,st,color='k',lw=1.4,label=lb)
    A_.axvline(ts,color='gray',lw=1,ls=':'); A_.axhline(lim,color='r',lw=1.2,label=f'threshold {lim} °C')
    A_.set_title(ttl,fontsize=10); A_.set_xlabel('time (s)'); A_.legend(fontsize=8.5,frameon=False,loc='upper right' if h1==500 else 'lower right')
ax[0].set_ylabel('surface temperature (°C)')
plt.tight_layout(); plt.savefig(OUT+'fig_c2.png',dpi=DPI); plt.close()

x=np.linspace(0,1,400); chk=[]
fig,ax=plt.subplots(1,2,figsize=(6.3,3.2))
for A_,(h1,h2,ts,tp,lim,ttl) in zip(ax,sc):
    B1,B2=h1*Lm/.2,h2*Lm/.2; Fs,Fp=al*ts/Lm**2,al*tp/Lm**2; z,w,c,Ac=two_stage(B1,B2,Fs)
    times=[0.25*ts,0.5*ts,ts,ts+(3 if ts==10 else 6),ts*2.5,(4*ts if ts==30 else 6*ts),tp]
    for t in times:
        th=theta(x,t*al/Lm**2,Fs,z,w,c,Ac); T=25+75*th
        chk.append((th.argmax()==len(x)-1,np.all(np.diff(th)>-1e-12)))
        A_.plot(x*4,T,lw=1.5,ls='-' if t<=ts else '--',label=f'{t:g} s')
    A_.set_title(ttl+f', switch {ts} s',fontsize=10); A_.set_xlabel('depth from insulated face (mm)')
    lo,hi=A_.get_ylim(); A_.set_ylim(lo,hi+0.28*(hi-lo))
    A_.legend(loc='upper left',fontsize=8,ncol=4,frameon=False,columnspacing=0.7,handlelength=1.4,handletextpad=0.4)
ax[0].set_ylabel('temperature (°C)')
plt.tight_layout(); plt.savefig(OUT+'fig_c2_profiles.png',dpi=DPI); plt.close()
print('profiles ok',all(a and b for a,b in chk),len(chk))
