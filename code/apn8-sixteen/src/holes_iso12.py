# isomorphism classes (under GL(8,2)) of the hole sets of the maximal 12-partial spreads in p12.txt, by backtracking
# Rowen (Claude Opus 5.5), 2026-10-25. usage: python holes_iso12.py [p12.txt]  (the file argument was added for this
# folder; Rowen's copy read runs/i5/p12.txt)
import itertools, collections, sys
def span(b):
    s={0}
    for v in b: s|={x^v for x in s}
    return s
def holes(t):
    k=int(t[1]); by=[int(x,16) for x in t[2:2+4*k]]
    cov=set()
    for i in range(k): cov|=span(by[4*i:4*i+4])
    return frozenset(set(range(1,256))-cov)
def ppc(H):  # planes through each point
    pp=collections.Counter()
    Hl=sorted(H)
    seen=set()
    for a,b in itertools.combinations(Hl,2):
        if a^b not in H: continue
        for c in Hl:
            p=frozenset(span([a,b,c])-{0})
            if len(p)==7 and p<=H and p not in seen:
                seen.add(p)
                for x in p: pp[x]+=1
    return pp
def iso(H1,H2,p1,p2,count=False):
    # basis of F_2^8 from H1, rarest planes-per-point values first
    freq=collections.Counter(p1[x] for x in H1)
    order=sorted(H1,key=lambda x:(freq[p1[x]],p1[x],x))
    basis=[];S={0}
    for x in order:
        if x not in S: basis.append(x); S=span(basis)
        if len(basis)==8: break
    n=0
    def rec(j,img,S):  # img: dict point->image over span of basis[:j]
        nonlocal n
        if j==8:
            n+=1; return not count
        b=basis[j]
        used=set(img.values())
        for y in H2:
            if p2[y]!=p1[b] or y in used: continue
            ok=True; new={}
            for x,ix in img.items():
                z=x^b; iz=ix^y
                if iz in used or iz==0 and z!=0: ok=False;break
                if (z in H1)!=(iz in H2) or (z in H1 and p1[z]!=p2[iz]): ok=False;break
                new[z]=iz
            if not ok: continue
            img2=dict(img); img2.update(new)
            if len(set(img2.values()))!=len(img2): continue
            if rec(j+1,img2,None): return True
        return False
    r=rec(0,{0:0},None)
    return n if count else r
lines=[l.split() for l in open(sys.argv[1] if len(sys.argv)>1 else 'runs/i5/p12.txt') if '-w0 ' in l]
Hs=[holes(t) for t in lines]; P=[ppc(H) for H in Hs]
cls=[]
for i,t in enumerate(lines):
    for c in cls:
        if iso(Hs[c[0]],Hs[i],P[c[0]],P[i]): c.append(i);break
    else: cls.append([i])
for c in cls:
    aut=iso(Hs[c[0]],Hs[c[0]],P[c[0]],P[c[0]],count=True)
    s=sum(aut//int(lines[i][0].split('-a')[1].split('-')[0]) for i in c)
    print('hole class: |Aut H| =',aut,' configs',[lines[i][0] for i in c],' sum |Aut H|/|Aut P| =',s,'(covers of the complement)')
print(len(cls),'hole classes')
