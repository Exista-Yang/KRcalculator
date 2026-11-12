import sys
sys.path.insert(0, '.')
from hu import hu, flows
def H(i, a, b): return [('E', i, b), ('F', i, a)]
for cells in [(1,1,1,1),(2,1,2,1),(2,2,2,2),(3,1,3,1),(3,3,2,1)]:
    x,y,z,w = cells
    w1 = H(0,x,y) + H(2,z,w); w2 = H(2,z,w) + H(0,x,y)
    print('far comm', cells, flows(list(cells), w1), 'unknowns d=0', hu(list(cells), w1, w2, 0))
# merged-block intermediate of CORE Lemma 3.12: [I+s_a+s_c-s_d-s_e]: merge all three cells then split reversed
for cells in [(2,2,2),(3,3,3),(2,1,2)]:
    x,y,z = cells
    merged = [('E',1,z),('E',0,y+z),('F',0,y+z),('F',1,y)]  # (x,y,z)->(x,y+z,0)->(x+y+z,0,0)->(z... check
    try:
        print('3-H vs merge-split', cells, flows(list(cells), H(0,x,y)+H(1,x,z)+H(0,y,z)), flows(list(cells), merged))
    except AssertionError as e:
        print('skip', cells, e)
