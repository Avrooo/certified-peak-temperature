"""Certified two-stage peak decision solver (default configuration = V7).
Default: 96-bit interval kernel, K=N=24 Taylor terms, sin/cos endpoint enclosure, pair-loop refactor.
Fallback: if the 96-bit run raises an unexpected error / times out / is flagged by a certification failure,
retry once with a 256-bit kernel. Invalid input (ValueError) is never retried. An error that survives the retry is
re-raised -- it is NOT converted to INCONCLUSIVE.
Frozen E4 baseline lives in frozen/ and *_frozenE4.py (unchanged)."""
import e4, e3
_fb=None
def _fallback_stack():
    global _fb
    if _fb is None:
        import stacks; _fb=stacks.load_stack("fb256",".",256,50,52,1,1,1)
    return _fb
def decide(box,Fs,Fp,L,fallback=True):
    """returns (decision, info) with decision in SAFE_ALL / EXCEEDS_ALL / INCONCLUSIVE; info['fallback'] says if 256-bit was used."""
    e4.validate_inputs(box,Fs,Fp,L)                       # ValueError, never retried
    e3.FAILS[0]=0
    try:
        d,w=e4.b6_decision_first(box,Fs,Fp,L)
        flagged=(d=='INCONCLUSIVE' and e3.FAILS[0]>0)
        if not flagged or not fallback: return d,dict(fallback=False,width=w,certfail=e3.FAILS[0])
    except ValueError: raise
    except Exception:
        if not fallback: raise
    st=_fallback_stack(); st.e3.FAILS[0]=0; stacks_clear=__import__('stacks').clear; stacks_clear(st)
    d,w=st.e4.b6_decision_first(box,Fs,Fp,L)             # an error here propagates
    return d,dict(fallback=True,width=w,certfail=st.e3.FAILS[0])
