"""Load independent copies of the (ivk, e3, e4) stack with their own precision/flags, so a fallback can use a different kernel."""
import os,sys,types,io,contextlib
def load_stack(name,srcdir,P,K,N,exp1,sc,pairs):
    env={"IVK_P":str(P),"IVK_K":str(K),"IVK_N":str(N),"IVK_EXP1":"1" if exp1 else "0","E3_SC":"1" if sc else "0","E3_PAIRS":"1" if pairs else "0"}
    old={k:os.environ.get(k) for k in env}; os.environ.update(env)
    try:
        ivk=types.ModuleType(name+"_ivk"); src=open(f"{srcdir}/ivk.py").read()
        exec(compile(src,name+"_ivk","exec"),ivk.__dict__)
        e3=types.ModuleType(name+"_e3"); src=open(f"{srcdir}/e3.py").read()
        names="fl,add,sub,mul,div,divi,muli,sq,neg,cn,PI,ONE,ZERO,exp_neg,sinc,cosT,sinT,hull"
        src=src.replace("import ivk as k\n","k=__IVK__\n").replace("from ivk import "+names,"fl,add,sub,mul,div,divi,muli,sq,neg,cn,PI,ONE,ZERO,exp_neg,sinc,cosT,sinT,hull=[getattr(k,n) for n in '"+names+"'.split(',')]")
        e3.__dict__["__IVK__"]=ivk; exec(compile(src,name+"_e3","exec"),e3.__dict__)
        e4=types.ModuleType(name+"_e4"); src=open(f"{srcdir}/e4.py").read()
        src=src.replace("import e3,ivk as k,e4ref as R","k=__IVK__; e3=__E3__; import e4ref as R")
        e4.__dict__["__IVK__"]=ivk; e4.__dict__["__E3__"]=e3
        with contextlib.redirect_stdout(io.StringIO()): exec(compile(src,name+"_e4","exec"),e4.__dict__)
    finally:
        for k_,v in old.items():
            if v is None: os.environ.pop(k_,None)
            else: os.environ[k_]=v
    return types.SimpleNamespace(name=name,ivk=ivk,e3=e3,e4=e4,P=P)
def clear(st):
    st.e3._rootcache.clear(); st.e3._modecache.clear()
    if hasattr(st.e3,'_sccache'): st.e3._sccache.clear()
