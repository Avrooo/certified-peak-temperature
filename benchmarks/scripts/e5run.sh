#!/bin/sh
R=3
python3 e5.py . $R "V1 fast mul only"
IVK_EXP1=1 python3 e5.py . $R "V2 +exp single-eval"
E3_SC=1 python3 e5.py . $R "V3 +sin/cos from root ends"
E3_PAIRS=1 python3 e5.py . $R "V4 +pair-loop refactor"
IVK_P=96 IVK_K=24 IVK_N=24 python3 e5.py . $R "V5 +P=96 (adaptive precision)"
IVK_P=128 IVK_K=32 IVK_N=32 python3 e5.py . $R "V5b +P=128"
IVK_EXP1=1 E3_SC=1 E3_PAIRS=1 python3 e5.py . $R "V6 V2+V3+V4"
IVK_EXP1=1 E3_SC=1 E3_PAIRS=1 IVK_P=96 IVK_K=24 IVK_N=24 python3 e5.py . $R "V7 all (P=96)"
