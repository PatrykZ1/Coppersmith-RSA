#!/usr/bin/env python

from typing import List, Callable
from gmpy2 import mpz, mpfr
import gmpy2
import random
import json
from poly import *

gmpy2.get_context().precision = 2048

def int_root(n: int, k: int):
    r, exact = gmpy2.iroot(mpz(n), k)
    return int(r), bool(exact)

def vec_dot(a: List[float], b: List[float]):
    s = 0.0
    for x, y in zip(a, b):
        s += x * y
    return s

def vec_len_sq(a: List[float]):
    s = 0.0
    for x in a:
        s += x * x
    return s

# Gram-Schmidt and LLL
def gram_schmidt(B: List[List[mpz]]):
    n = len(B)
    if n == 0:
        return [], [], []
    m = len(B[0])
    Bstar = [[mpfr(0)] * m for _ in range(n)]
    mu = [[mpfr(0)] * n for _ in range(n)]
    norm_sq = [mpfr(0)] * n
    for i in range(n):
        v = [mpfr(x) for x in B[i]]
        for j in range(i):
            denom = norm_sq[j]
            if denom == 0:
                muij = mpfr(0)
            else:
                dot_prod = sum(v[t] * Bstar[j][t] for t in range(m))
                muij = dot_prod / denom
            mu[i][j] = muij
            for t in range(m):
                v[t] -= muij * Bstar[j][t]
        Bstar[i] = v
        norm_sq[i] = sum(x * x for x in v)
    return mu, norm_sq

def lll_reduce(mat: List[List[int]], delta: float = 0.75):
    B = [[mpz(x) for x in row] for row in mat]
    n = len(B)
    if n == 0:
        return []
    mu, norm_sq = gram_schmidt(B)
    k = 1
    delta_mpfr = mpfr(delta)
    while k < n:
        for j in range(k - 1, -1, -1):
            if abs(mu[k][j]) > 0.5:
                q = mpz(gmpy2.rint(mu[k][j]))
                for i in range(len(B[k])):
                    B[k][i] -= q * B[j][i]
                mu, norm_sq = gram_schmidt(B)
        if norm_sq[k] >= (delta_mpfr - mu[k][k - 1] ** 2) * norm_sq[k - 1]:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            mu, norm_sq = gram_schmidt(B)
            k = max(k - 1, 1)
    return [[int(x) for x in row] for row in B]

def build_f_poly(M0: int, e: int, C: int) -> Poly:
    f = Poly([M0, 1]) ** e
    f = f - Poly([C])
    return f

def polys_for_coppersmith(f: Poly, N: int, s: int, t: int):
    polys = []
    fi = Poly([1])
    for i in range(s + 1):
        Ni = pow(N, s - i)
        scaled_fi = fi * Ni
        for j in range(t + 1):
            polys.append(scaled_fi.shift(j))
        fi = fi * f
    return polys

def poly_to_scaled_vector(p: Poly, X: int, max_deg: int):
    vec = [0] * (max_deg + 1)
    for k, a in enumerate(p.coefficients):
        if k <= max_deg:
            vec[k] = int(a * pow(X, k))
    return vec

def build_lattice_matrix(polys: List[List[int]], X: int):
    degs = [p.degree for p in polys]
    max_deg = max(degs) if degs else 0
    mat = []
    for p in polys:
        mat.append(poly_to_scaled_vector(p, X, max_deg))
    return mat

def vector_to_poly(vec: List[int], X: int):
    coeffs = []
    for k, v in enumerate(vec):
        c, rem = gmpy2.f_divmod(mpz(v), mpz(X) ** k)
        coeffs.append(c)
    return Poly(coeffs)

def find_correct_integer_root(poly: Poly, bound: int, test_func: Callable[[int], bool]):
    for x in range(-min(bound, int(1e7)), min(bound, int(1e7)) + 1):
        if poly(x) == 0 and test_func(x):
            return x
    return None

def coppersmith_univariate(N: int, e: int, C: int, M0: int, s: int = 2, t: int = 5, delta: float = 0.75):
    f = build_f_poly(M0, e, C)
    X, _ = gmpy2.iroot(mpz(N), e)
    X = int(X) + 1
    polys = polys_for_coppersmith(f, N, s, t)
    mat = build_lattice_matrix(polys, X)
    reduced = lll_reduce(mat, delta=delta)
    def len_sq(v: List[int]):
        s = 0
        for x in v:
            s += int(x) * int(x)
        return s
    reduced_sorted = sorted(reduced, key=len_sq)
    for vec in reduced_sorted[:min(10, len(reduced_sorted))]:
        g = vector_to_poly(vec, X)
        root = find_correct_integer_root(g, X, lambda r: (f(r) % N) == 0)
        if root is not None:
            return root
    return None

def gen_rsa_manual(bits=256, e=3):
    assert bits % 2 == 0
    while True:
        p = gmpy2.next_prime(random.getrandbits(bits // 2))
        q = gmpy2.next_prime(random.getrandbits(bits // 2))
        if p != q and gmpy2.gcd((p - 1)*(q - 1), e) == 1:
            return int(p * q), int(p), int(q)

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

def decrypt(N, e, C, M0):
    x = coppersmith_univariate(N, e, C, M0, s=2, t=4, delta=0.99)
    return x

def open_file_and_decrypt(filename):
    data = {}
    with open(filename, 'r') as f:
        data = json.load(f)
    if not 'N' in data or not 'e' in data or not 'C' in data or not 'M0' in data:
        print('File does not contain proper data.')
        return 
    N = data['N']
    e = data['e']
    C = data['C']
    M0 = data['M0']
    x = decrypt(N, e, C, M0)
    M = M0 + x
    print(f'Given message prefix: {M0}\nRecovered message suffix: {x}\nComplete message: {M}')

if __name__ == "__main__":
    f1 = 'encrypted1.txt'
    open_file_and_decrypt(f1)
