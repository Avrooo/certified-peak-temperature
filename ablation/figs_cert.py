import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.labelsize':10,'axes.titlesize':10,'xtick.labelsize':9.5,'ytick.labelsize':9.5,'legend.fontsize':9,
 'axes.spines.top':False,'axes.spines.right':False,'mathtext.default':'regular','font.family':'DejaVu Sans'})
NAVY='#1f4e79'; GOLD='#d9a441'; GREY='#6b7078'; RED='#b23a2e'
d=json.load(open('floorrun.json')); floors=['0.001','1e-05','1e-07','1e-09','1e-11','1e-13']
edges=[1e-2,1e-4,1e-6,1e-8,1e-10,1e-12,1e-14]
def b(s):
    for k,e in enumerate(edges):
        if s>=e: return k
    return len(edges)
n=len(edges)+1
tab={}; 
for f in floors:
    rows=d[f]; dec=np.zeros(n,int); tot=np.zeros(n,int); 
    for dd,w,t,s in rows:
        k=b(s); tot[k]+=1; dec[k]+= dd in('SAFE_ALL','EXCEEDS_ALL','MIXED')
    tab[f]=(dec,tot)
print({f:int(tab[f][0].sum()) for f in floors}); print(tab['1e-13'][1])
for f in floors: print(f,np.round(tab[f][0]/tab[f][1],3))
fig,ax=plt.subplots(1,2,figsize=(6.3,3.1))
a=ax[0]; x=np.arange(len(floors)); dec=[tab[f][0].sum() for f in floors]; inc=[1600-v for v in dec]
a.bar(x,dec,color=NAVY,label='decided'); a.bar(x,inc,bottom=dec,color=GOLD,label='inconclusive')
for i,v in enumerate(dec): a.text(i,v/2,f'{v}',ha='center',color='white',fontsize=8.5)
for i,v in enumerate(inc): a.text(i,dec[i]+v/2,f'{v}',ha='center',fontsize=8.5)
a.set_xticks(x); a.set_xticklabels([r'$10^{-3}$',r'$10^{-5}$',r'$10^{-7}$',r'$10^{-9}$',r'$10^{-11}$',r'$10^{-13}$'],fontsize=8.5); a.set_xlabel('finest tolerance allowed'); a.set_ylabel('cases (of 1,600)'); a.set_title('(a) decisions by tolerance floor',fontsize=10)
a.set_ylim(0,1900); a.legend(loc='upper center',ncol=2,frameon=False,fontsize=8.5)
b_=ax[1]; xs=np.arange(n)
for f,c,m,lab in (('1e-05',GREY,'s',r'$10^{-5}$'),('1e-09',GOLD,'^',r'$10^{-9}$'),('1e-13',NAVY,'o',r'$10^{-13}$')):
    dec,tot=tab[f]; b_.plot(xs,100*dec/tot,marker=m,color=c,lw=1.5,ms=4.5,label=lab)
labs=[r'$\geq10^{-2}$',r'$10^{-4}$',r'$10^{-6}$',r'$10^{-8}$',r'$10^{-10}$',r'$10^{-12}$',r'$10^{-14}$',r'$<10^{-14}$']
b_.set_xticks(xs); b_.set_xticklabels(labs,rotation=40,ha='right',fontsize=8.5); b_.set_xlabel('gap (lower end of bin)'); b_.set_ylabel('cases decided (%)'); b_.set_title('(b) decided share by gap',fontsize=10)
b_.set_ylim(-3,108); b_.legend(title='tolerance floor',title_fontsize=8.5,frameon=False,fontsize=8.5,loc='lower left')
plt.tight_layout(); plt.savefig(__import__('os').environ.get('OUTDIR','.')+'/fig_cert.png',dpi=900); plt.close()
json.dump({f:[tab[f][0].tolist(),tab[f][1].tolist()] for f in floors},open('floor_by_gap.json','w'))
