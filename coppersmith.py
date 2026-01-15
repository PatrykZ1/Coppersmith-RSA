from __future__ import annotations
from typing import List, Tuple, Optional
import math
import sys
import gmpy2
from gmpy2 import mpz


def poly_add(a: List[int], b: List[int]):
    n = max(len(a), len(b))
    res = [0] * n
    for i in range(n):
        ai = a[i] if i < len(a) else 0
        bi = b[i] if i < len(b) else 0
        res[i] = ai + bi
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def poly_sub(a: List[int], b: List[int]):
    n = max(len(a), len(b))
    res = [0] * n
    for i in range(n):
        ai = a[i] if i < len(a) else 0
        bi = b[i] if i < len(b) else 0
        res[i] = ai - bi
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def poly_scalar_mul(a: List[int], k: int):
    if k == 0:
        return [0]
    return [int(ai * k) for ai in a]

def poly_mul(a: List[int], b: List[int]):
    if not a or not b:
        return [0]
    res = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            res[i + j] += ai * bj
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def poly_pow(a: List[int], e: int):
    res = [1]
    base = a[:]
    while e > 0:
        if e & 1:
            res = poly_mul(res, base)
        base = poly_mul(base, base)
        e >>= 1
    return res

def poly_shift(a: List[int], s: int):
    if a == [0]:
        return [0]
    return ([0] * s) + a[:]

def poly_degree(a: List[int]):
    return len(a) - 1

def poly_eval(a: List[int], x: int):
    # Horner
    res = 0
    for coeff in reversed(a):
        res = res * x + coeff
    return res

def int_root(n: int, k: int):
    if gmpy2:
        r, exact = gmpy2.iroot(mpz(n), k)
        return int(r), bool(exact)
    else:
        r = int(n ** (1.0 / k))
        while (r + 1) ** k <= n:
            r += 1
        while r ** k > n:
            r -= 1
        return r, (r ** k == n)

def vec_dot_float(a: List[float], b: List[float]):
    s = 0.0
    for x, y in zip(a, b):
        s += x * y
    return s

def vec_len_sq_float(a: List[float]):
    s = 0.0
    for x in a:
        s += x * x
    return s

# Gram-Schmidt and LLL
def gram_schmidt_float(B: List[List[float]]):
    n = len(B)
    if n == 0:
        return [], [], []
    m = len(B[0])
    Bstar = [[0.0] * m for _ in range(n)]
    mu = [[0.0] * n for _ in range(n)]
    norm_sq = [0.0] * n
    for i in range(n):
        v = B[i][:]
        for j in range(i):
            denom = norm_sq[j]
            if denom == 0.0:
                muij = 0.0
            else:
                muij = vec_dot_float(v, Bstar[j]) / denom
            mu[i][j] = muij
            for t in range(m):
                v[t] -= muij * Bstar[j][t]
        Bstar[i] = v
        norm_sq[i] = vec_len_sq_float(Bstar[i])
    return Bstar, mu, norm_sq

def lll_reduce_float(mat: List[List[int]], delta: float = 0.75):
    B = [list(map(int, row))[:] for row in mat]
    n = len(B)
    if n == 0:
        return []
    m = len(B[0])
    Bf = [list(float(x) for x in row) for row in B]
    Bstar, mu, norm_sq = gram_schmidt_float(Bf)
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = int(round(mu[k][j]))
            if q != 0:
                B[k] = [int(B[k][i] - q * B[j][i]) for i in range(m)]
                Bf[k] = [float(B[k][i]) for i in range(m)]
                Bstar, mu, norm_sq = gram_schmidt_float(Bf)
        if norm_sq[k] >= (delta - mu[k][k - 1] * mu[k][k - 1]) * norm_sq[k - 1]:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            Bf[k], Bf[k - 1] = Bf[k - 1], Bf[k]
            Bstar, mu, norm_sq = gram_schmidt_float(Bf)
            k = max(k - 1, 1)
    return B

def build_f_poly(M0: int, e: int, C: int):
    # f(x) = (M0 + x)^e - C
    coefs = [0] * (e + 1)
    for k in range(e + 1):
        binom = math.comb(e, k)
        coefs[k] = binom * pow(M0, e - k)
    coefs[0] -= C
    while len(coefs) > 1 and coefs[-1] == 0:
        coefs.pop()
    return coefs

def polys_for_coppersmith(f: List[int], N: int, s: int, t: int):
    polys: List[List[int]] = []
    f_pows = [[1]]
    for i in range(1, s + 1):
        f_pows.append(poly_mul(f_pows[-1], f))
    for i in range(0, s + 1):
        fi = f_pows[i] if i < len(f_pows) else poly_pow(f, i)
        Ni = pow(N, s - i)
        scaled_fi = poly_scalar_mul(fi, Ni)
        for j in range(0, t + 1):
            p = poly_shift(scaled_fi, j)
            polys.append(p)
    return polys

def poly_to_scaled_vector(p: List[int], X: int, max_deg: int):
    vec = [0] * (max_deg + 1)
    for k, a_k in enumerate(p):
        if k <= max_deg:
            vec[k] = int(a_k * pow(X, k))
    return vec

def build_lattice_matrix(polys: List[List[int]], X: int):
    degs = [len(p) - 1 for p in polys]
    max_deg = max(degs) if degs else 0
    mat = []
    for p in polys:
        mat.append(poly_to_scaled_vector(p, X, max_deg))
    return mat, max_deg

def vector_to_poly(vec: List[int], X: int):
    a = []
    for k, v in enumerate(vec):
        denom = pow(X, k)
        ak = int(v // denom) if denom != 0 else 0
        a.append(ak)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a

def find_integer_roots_bruteforce(poly: List[int], bound: int):
    roots = []
    brute = min(bound, 2000)
    for cand in range(-brute, brute + 1):
        if poly_eval(poly, cand) == 0:
            roots.append(cand)
    return roots

# main attack function
def coppersmith_univariate(N: int, e: int, C: int, M0: int, s: int = 2, t: int = 5, delta: float = 0.75):
    f = build_f_poly(M0, e, C)
    X = int(math.floor(N ** (1.0 / e))) + 1
    polys = polys_for_coppersmith(f, N, s, t)
    mat, max_deg = build_lattice_matrix(polys, X)
    reduced = lll_reduce_float(mat, delta=delta)
    def len_sq(v: List[int]):
        s = 0
        for x in v:
            s += int(x) * int(x)
        return s
    reduced_sorted = sorted(reduced, key=len_sq)
    for vec in reduced_sorted[:min(10, len(reduced_sorted))]:
        g = vector_to_poly(vec, X)
        roots = find_integer_roots_bruteforce(g, X)
        for r in roots:
            if abs(r) < X:
                if poly_eval(f, r) % N == 0:
                    return r
    return None

def demo_small():
    p, q = 2137, 21372157
    N = p * q
    e = 3
    M0 = 0xdeadbeef
    x_true = 42
    M = M0 + x_true
    C = pow(M, e, N)
    print("N =", N, " e =", e)
    print("M0 =", M0, " x_true =", x_true, " M =", M, " C =", C)
    x_found = coppersmith_univariate(N, e, C, M0, s=2, t=4, delta=0.99)
    print("x_found:", x_found)
    if x_found is not None:
        print("Recovered M:", M0 + x_found)
    else:
        print("Nie znaleziono x (prawdopodobnie poza boundem)")

if __name__ == "__main__":
    demo_small()
