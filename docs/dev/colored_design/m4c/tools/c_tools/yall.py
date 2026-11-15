import itertools, math, sys
def tables(lam, mu):
    m = len(lam)
    def rec(i, colrem):
        if i == m:
            if all(c == 0 for c in colrem): yield []
            return
        for row in itertools.product(*[range(min(lam[i], c) + 1) for c in colrem]):
            if sum(row) == lam[i]:
                for rest in rec(i + 1, [c - r for c, r in zip(colrem, row)]):
                    yield [list(row)] + rest
    yield from rec(0, list(mu))
def multinom(n, parts):
    r = math.factorial(n)
    for p in parts: r //= math.factorial(p)
    return r
def yflows(t):
    f = 1
    for row in t: f *= multinom(sum(row), [x for x in row if x])
    cells = [(i, j, t[i][j]) for i in range(len(t)) for j in range(len(t[0])) if t[i][j]]
    for (i, j, a), (i2, j2, b) in itertools.combinations(cells, 2):
        if i < i2 and j > j2: f *= math.comb(a + b, a)
    return f
if __name__ == "__main__":
    lam = [3,3,3]
    for k in (3, 4):
        best = []
        for mu in itertools.product(range(10), repeat=k):
            if sum(mu) != 9 or 0 in mu: continue
            for t in tables(lam, mu):
                best.append((yflows(t), mu, t))
        best.sort(reverse=True)
        print(k, 'columns: max Y-form flows', [(b[0], b[1], b[2]) for b in best[:3]])
        print('  count > 1e5:', sum(1 for b in best if b[0] > 1e5), 'of', len(best))
