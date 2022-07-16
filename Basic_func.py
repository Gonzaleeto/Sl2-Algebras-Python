import numpy as np
import sympy as sp
import math  as mt

def binomial(n, r):
    ''' Binomial coefficient, nCr, aka the "choose" function 
        n! / (r! * (n - r)!)
    '''
    p = 1    
    for i in range(1, min(r, n - r) + 1):
        p *= n
        p //= i
        n -= 1
    return p
