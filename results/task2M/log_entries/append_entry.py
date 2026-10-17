import sys, textwrap, subprocess, datetime
es_src, en_src = sys.argv[1], sys.argv[2]
SE='/home/claude/work/repo/preregistration/original_sealed_es/PREREGISTRO_1j_1k.txt'
EN='/home/claude/work/repo/preregistration/PREREGISTRATION_LOG_EN.txt'
HS='/home/claude/work/repo/preregistration/original_sealed_es/PREREGISTRO_hashes.txt'
D=datetime.datetime.utcnow().strftime('%a %b %e %H:%M:%S UTC %Y')
W=lambda t,w: textwrap.wrap(t,w,break_long_words=False,break_on_hyphens=False)
es=[p.replace('__FECHA__',D) for p in open(es_src).read().rstrip('\n').split('\n')]
en=[p.replace('__FECHA__',D) for p in open(en_src).read().rstrip('\n').split('\n')]
assert len(es)==len(en)
les,len_=[''],['']
for a,b in zip(es,en):
    A=W(a,120); assert all(ord(c)<128 for c in a)
    for w in range(80,500):
        Bl=W(b,w)
        if len(Bl)==len(A): break
    assert len(Bl)==len(A),(b[:40])
    les+=A; len_+=Bl
open(SE,'a').write('\n'.join(les)+'\n'); open(EN,'a').write('\n'.join(len_)+'\n')
h=subprocess.check_output(['sha256sum',SE]).decode().split()[0]
open(HS,'a').write(f"{D}\n{h}  PREREGISTRO_1j_1k.txt\n")
a=open(SE).read().split('\n'); b=open(EN).read().split('\n')
print(D, h, len(a), len(b), all((x=='')==(y=='') for x,y in zip(a,b)))
