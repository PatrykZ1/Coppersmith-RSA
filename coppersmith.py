# coppersmith_with_own_lll.py
# Edukacyjna implementacja ataku Coppersmitha (jednowymiarowy mały pierwiastek)
# z autonomiczną implementacją LLL (dokładna arytmetyka przez Fraction).
#
# Uwaga: kod jest przeznaczony do eksperymentów edukacyjnych na małych liczbach.
# Dla dużych N (np. 1024/2048-bit) konieczne są zoptymalizowane biblioteki (fpylll/Sage).
#
# Wymagane: sympy
# Opcjonalne (przyspieszenie/duża arytmetyka): gmpy2
#
# Python >= 3.8

from __future__ import annotations
from fractions import Fraction
from typing import List, Tuple, Optional
import math
import sys

try:
    import sympy as sp
except Exception:
    raise ImportError("sympy required: pip install sympy")

try:
    import gmpy2
    from gmpy2 import mpz
except Exception:
    gmpy2 = None
    mpz = int

# ---------------------------
# Pomocnicze funkcje arytmetyczne
# ---------------------------
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

def inner_product(a: List[int], b: List[int]) -> int:
    return sum(int(x) * int(y) for x, y in zip(a, b))

def vector_sub(a: List[int], b: List[int], scale: int = 1):
    return [int(x - scale * y) for x, y in zip(a, b)]

def vector_add(a: List[int], b: List[int], scale: int = 1):
    return [int(x + scale * y) for x, y in zip(a, b)]

# ---------------------------
# Dokładna implementacja Gram-Schmidt (Fraction)
# ---------------------------
def gram_schmidt_exact(basis: List[List[int]]):
    """
    Dla zadanej listy wektorów bazowych (wiersze macierzy) zwraca:
      B_star: listę ortogonalnych wektorów (Fraction),
      mu: macierz współczynników mu[i][j] = <b_i, b*_j> / <b*_j, b*_j>,
      Bstar_norm_sq: listę norm kwadratów ||b*_i||^2 jako Fraction.
    Wszystkie obliczenia w Fraction -> pełna dokładność racjonalna.
    """
    n = len(basis)
    if n == 0:
        return [], [], []
    m = len(basis[0])
    Bstar: List[List[Fraction]] = [[Fraction(0) for _ in range(m)] for _ in range(n)]
    mu: List[List[Fraction]] = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    norm_sq: List[Fraction] = [Fraction(0) for _ in range(n)]

    for i in range(n):
        # v = b_i as Fraction vector
        v = [Fraction(int(x)) for x in basis[i]]
        for j in range(i):
            # mu[i][j] = <b_i, b*_j> / ||b*_j||^2
            num = sum(vk * bk for vk, bk in zip(v, Bstar[j]))
            denom = norm_sq[j]
            if denom == 0:
                mu_ij = Fraction(0)
            else:
                mu_ij = Fraction(num, denom)
            mu[i][j] = mu_ij
            # v = v - mu_ij * b*_j
            v = [vi - mu_ij * bj for vi, bj in zip(v, Bstar[j])]
        Bstar[i] = v
        norm_sq[i] = sum(x * x for x in v)
    return Bstar, mu, norm_sq

# ---------------------------
# Implementacja LLL (dokładna / oparta na Fraction)
# ---------------------------
def lll_reduce_exact(basis: List[List[int]], delta: Fraction = Fraction(3, 4)):
    """
    Pełna implementacja algorytmu LLL z dokładną arytmetyką.
    delta (Fraction) powinien należeć do (1/4, 1)
    Zwrot: zredukowana baza (lista wektorów int).
    """
    n = len(basis)
    if n == 0:
        return []
    # Kopia bazy (będziemy ją modyfikować)
    B = [list(map(int, row))[:] for row in basis]

    # Główna pętla
    k = 1
    Bstar, mu, norm_sq = gram_schmidt_exact(B)
    # Przechodzimy aż k == n
    while k < n:
        # Size reduction: dla j = k-1..0
        for j in range(k - 1, -1, -1):
            # round mu[k][j] do najbliższej liczby całkowitej
            m = mu[k][j]
            # round m (Fraction) to nearest int
            if m >= 0:
                q = int(m + Fraction(1, 2))
            else:
                q = int(m - Fraction(1, 2))
            if q != 0:
                # b_k = b_k - q * b_j
                B[k] = vector_sub(B[k], B[j], q)
                # recompute GS for row k (najprościej cała GS)
                Bstar, mu, norm_sq = gram_schmidt_exact(B)
        # Lovasz condition:
        left = norm_sq[k]
        # compute mu[k][k-1]^2
        mu_k_k1 = mu[k][k - 1]
        right = (delta - mu_k_k1 * mu_k_k1) * norm_sq[k - 1]
        if left >= right:
            k += 1
        else:
            # swap b_k and b_{k-1}
            B[k], B[k - 1] = B[k - 1], B[k]
            # recompute GS
            Bstar, mu, norm_sq = gram_schmidt_exact(B)
            k = max(k - 1, 1)
    return B

# ---------------------------
# Narzędzia: wielomiany / macierz kratowa
# ---------------------------
def build_f(M0: int, e: int, C: int):
    x = sp.symbols('x')
    return sp.Poly(sp.expand((M0 + x) ** e - C), x)

def polys_for_coppersmith(f: sp.Poly, N: int, s: int, t: int):
    x = f.gen
    polys: List[sp.Poly] = []
    # klasyczna konstrukcja p_{i,j} = x^j * f(x)^i * N^{s-i}
    for i in range(0, s + 1):
        fi = sp.expand(f.as_expr() ** i)
        Ni = pow(N, s - i)
        for j in range(0, t + 1):
            expr = (x ** j) * fi * Ni
            polys.append(sp.Poly(sp.expand(expr), x))
    return polys

def poly_to_scaled_vector(p: sp.Poly, X: int, max_deg: int):
    x = p.gen
    # podstawienie x -> x * X
    expr = sp.expand(p.as_expr().subs(x, x * X))
    # wyznacz współczynniki 0..max_deg
    coeffs = [0] * (max_deg + 1)
    poly_expr = sp.Poly(expr, x)
    for i in range(0, max_deg + 1):
        c = poly_expr.coeff_monomial(x ** i)
        coeffs[i] = int(sp.expand(c))
    return coeffs

def build_lattice_matrix(polys: List[sp.Poly], X: int):
    degs = [p.degree() for p in polys]
    max_deg = max(degs) if degs else 0
    mat: List[List[int]] = []
    for p in polys:
        mat.append(poly_to_scaled_vector(p, X, max_deg))
    return mat, max_deg

def vector_to_poly(vec: List[int], X: int, x_symbol=None) -> sp.Poly:
    if x_symbol is None:
        x_symbol = sp.symbols('x')
    expr = 0
    for k, v in enumerate(vec):
        # oryginalny współczynnik a_k = v // X^k (integer division)
        denom = pow(X, k)
        if denom == 0:
            a_k = 0
        else:
            a_k = int(v // denom)
        expr += a_k * (x_symbol ** k)
    return sp.Poly(sp.expand(expr), x_symbol)

# ---------------------------
# Szukanie pierwiastków i weryfikacja
# ---------------------------
def find_integer_roots(poly: sp.Poly, bound: int):
    roots: List[int] = []
    if poly.degree() == 0:
        return roots
    # spróbuj numerycznych korzeni i filtracji do intów
    try:
        approx_roots = [complex(r) for r in sp.nroots(poly.as_expr())]
        for ar in approx_roots:
            if abs(ar.imag) < 1e-8:
                r_real = ar.real
                r_int = int(round(r_real))
                if abs(r_real - r_int) < 1e-6 and abs(r_int) < bound:
                    if poly.eval(r_int) == 0:
                        roots.append(r_int)
    except Exception:
        pass
    # dodatkowy bruteforce dla małych wartości
    brute = min(bound, 2000)
    for cand in range(-brute, brute + 1):
        if poly.eval(cand) == 0:
            roots.append(cand)
    # unikatowe
    return list(set(roots))

# ---------------------------
# Główna funkcja ataku
# ---------------------------
def coppersmith_univariate(N: int, e: int, C: int, M0: int,
                           s: int = 2, t: int = 5, delta: float = 0.75):
    x = sp.symbols('x')
    f = build_f(M0, e, C)
    # typowy bound X
    X = int(math.floor(N ** (1.0 / e))) + 1
    polys = polys_for_coppersmith(f, N, s, t)
    mat, max_deg = build_lattice_matrix(polys, X)
    # uruchom LLL (dokładna implementacja)
    # macierz powinna mieć wiersze jako wektory współczynników
    reduced = lll_reduce_exact(mat, Fraction(delta).limit_denominator())
    # sprawdź kilka najkrótszych wektorów (sortowanie po euklidesowej długości)
    def vec_len_sq(v):
        return sum(int(x) * int(x) for x in v)
    reduced_sorted = sorted(reduced, key=vec_len_sq)
    for vec in reduced_sorted[:min(10, len(reduced_sorted))]:
        g = vector_to_poly(vec, X, x)
        # uprość g
        g = sp.Poly(sp.expand(g.as_expr()), x)
        roots = find_integer_roots(g, X)
        for r in roots:
            if abs(r) < X:
                val = int(sp.Mod(f.eval(r), N))
                if val == 0:
                    return r
    return None

# ---------------------------
# Demo (mały przykład)
# ---------------------------
def demo_small():
    p, q = 2137, 21372157
    N = p * q
    e = 3
    M0 = 100*128 # znany prefiks
    x_true = 42
    M = M0 + x_true
    C = pow(M, e, N)
    print("N =", N, "e =", e)
    print("M0 =", M0, "x_true =", x_true, "M =", M, "C =", C)
    x_found = coppersmith_univariate(N, e, C, M0, s=2, t=4, delta=0.75)
    print("x_found:", x_found)
    if x_found is not None:
        print("Recovered M:", M0 + x_found)
    else:
        print("Nie znaleziono x (prawdopodobnie poza boundem)")

if __name__ == "__main__":
    demo_small()
