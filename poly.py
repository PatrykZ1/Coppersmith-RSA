#!/usr/bin/env python
        
class Poly:
    def __init__(self, coefficients=None):
        self.coefficients = coefficients[:] if coefficients else [0]
        self.fix_degree()

    def fix_degree(self):
        while len(self.coefficients) > 1 and self.coefficients[-1] == 0:
            self.coefficients.pop()
            
    def __add__(self, other):
        result = Poly()
        a = self.coefficients
        b = other.coefficients
        n = max(len(a), len(b))
        res = [0] * n
        for i in range(n):
            ai = a[i] if i < len(a) else 0
            bi = b[i] if i < len(b) else 0
            res[i] = ai + bi
        return Poly(res)

    def __sub__(self, other):
        return self + (-other)
        
    def __neg__(self):
        result = Poly(self.coefficients[:])
        for i in range(len(self.coefficients)):
            result.coefficients[i] = -result.coefficients[i]
        return result

    def __mul__(self, other):
        if isinstance(other, int):
            if other == 0:
                return Poly([0])
            return Poly([c * other for c in self.coefficients])

        if isinstance(other, Poly):
            a = self.coefficients
            b = other.coefficients
            if a == [0] or b == [0]:
                return Poly([0])

            res = [0] * (len(a) + len(b) - 1)
            for i, ai in enumerate(a):
                if ai == 0:
                    continue
                for j, bj in enumerate(b):
                    res[i + j] += ai * bj
            return Poly(res)

        return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __pow__(self, e: int):
        if e < 0:
            raise ValueError("negative powers not supported")

        res = Poly([1])
        base = Poly(self.coefficients)

        while e > 0:
            if e & 1:
                res = res * base
            base = base * base
            e >>= 1

        return res

    def shift(self, s: int):
        if self.coefficients == [0]:
            return Poly([0])
        return Poly(([0] * s) + self.coefficients)

    @property
    def degree(self):
        return len(self.coefficients) - 1

    def __call__(self, x):
        res = 0
        for c in reversed(self.coefficients):
            res = res * x + c
        return res

    def __repr__(self):
        return f"Poly({self.coefficients})"
