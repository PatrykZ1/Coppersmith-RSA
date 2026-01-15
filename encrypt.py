from Crypto.Util import number
import json

b = 32
e = 3

def encrypt(m, filename):
    p, q = number.getPrime(b), number.getPrime(b)
    N = p * q
    m_prefix = int(m * 0.1)
    m_suffix = m - m_prefix
    C = pow(m, e, N)
    with open(filename, 'w') as f:
        data = {'N' = N, 'e' = e, 'C' = C, 'M0' = m_prefix}
        j.dump(data, f)

filename = 'encrypted1.txt'
message = 123

encrypt(message, filename)
