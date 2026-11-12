import sys, itertools
sys.path.insert(0, '.')
from hu import hu
from yall import tables
def H(i, a, b): return [('E', i, b), ('F', i, a)]
memo = {}
def cost(x, y, z):
    if (x,y,z) not in memo:
        lab=[x,y,z]; w1=H(0,x,y)+H(1,x,z)+H(0,y,z); w2=H(1,y,z)+H(0,x,z)+H(1,x,y)
        memo[(x,y,z)] = hu(lab, w1, w2, 0)
    return memo[(x,y,z)]
lam=[3,3,3]
for mus, name in [([(3,3,3)], '(lam,lam)'), ([m for m in itertools.product(range(1,8),repeat=3) if sum(m)==9], 'all 3-col levels')]:
    n=0; big=0; ex=[]
    for mu in mus:
        for t in tables(lam, mu):
            n+=1
            a,b,c = t[0][2], t[1][1], t[2][0]
            if a and b and c and cost(a,b,c) > 4e5:
                big+=1; ex.append(('/'.join(','.join(map(str,r)) for r in t), (a,b,c), cost(a,b,c)))
    print(name, 'tables', n, 'with infeasible braid window', big, ex[:4])
