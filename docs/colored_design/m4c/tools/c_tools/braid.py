import sys
sys.path.insert(0, '.')
from hu import hu, flows
def H(i, a, b):   # atomic crossing of cells (a on i, b on i+1) -> (b, a): E_i^(b) F_i^(a)
    return [('E', i, b), ('F', i, a)]
def braid(x, y, z):
    lab = [x, y, z]
    w1 = H(0, x, y) + H(1, x, z) + H(0, y, z)      # (x,y,z)->(y,x,z)->(y,z,x)->(z,y,x)
    w2 = H(1, y, z) + H(0, x, z) + H(1, x, y)      # (x,y,z)->(x,z,y)->(z,x,y)->(z,y,x)
    return lab, w1, w2
for cells in [(1,1,1),(1,2,1),(2,1,2),(2,2,1),(2,2,2),(3,1,3),(1,3,1),(3,2,3),(2,3,2),(3,3,3)]:
    lab, w1, w2 = braid(*cells)
    print(cells, 'flows', flows(list(lab), w1), flows(list(lab), w2), 'unknowns d=0:', hu(lab, w1, w2, 0))
# E2 unit/counit certification windows (2 uprights), degree ab
for a, b in [(3,3),(3,6),(6,3),(4,5),(5,4)]:
    w = [('E', 0, b), ('F', 0, b)]
    print('E2 (a,b)=', (a,b), 'flows', flows([a,b], w), 'Hom(empty, w) d=ab unknowns', hu([a,b], [], w, a*b))
