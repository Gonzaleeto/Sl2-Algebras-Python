import numpy as np
import sympy as sp
import math  as mt
import copy # Para usar deepcopy en algunos objetos

from sympy import MatrixSymbol, Identity

# Para imprimir en pantallas mas legible los resultados

from sympy import init_printing
init_printing(use_latex=True)

##############################################################################

def Completa_Tabla_Lie(T):
    dim = len(T[1,2])
    ind = sp.Array(T.keys())
    w = sp.zeros(1,dim)
    TT = copy.deepcopy(T)
    for i in range(dim):
        for j in range(i+1,dim):
            if not (i+1,j+1) in T:
                TT[i+1,j+1] = w
    for i in range(dim):
        TT[i+1,i+1] = w
        for j in range(i+1,dim):
            # La aritmetica se tiene que hacer elemento por elemento
            w2 = sp.zeros(1,dim)
            for k in range(dim):
                w2[k] = - TT[i+1,j+1][k]
            TT[j+1,i+1] = w2
    return TT

def Killing_Form(T):
    dim = len(T[1,2])
    M = sp.Matrix.zeros(dim,dim)
    TT = Completa_Tabla_Lie(T)
    for i in range(dim):
        for j in range(dim):
            s = 0
            for k in range(dim):
                for l in range(dim):
                    s = s + TT[j+1,l+1][k]*TT[i+1,k+1][l]
            M[j,i] = s
    return M


