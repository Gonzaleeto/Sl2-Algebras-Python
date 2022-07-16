import numpy as np
import sympy as sp
import math  as mt
import copy # Para usar deepcopy en algunos objetos

from sympy import MatrixSymbol, Identity

# Para imprimir en pantallas mas legible los resultados

from sympy import init_printing
init_printing(use_latex=True)

# --------------------------------------------- #

from Relaciones_productos import *
from Clebsch_Gordan import *

##############################################################################

def Corrige_Tabla(S,T):
    # Se hace deepcopy para no modificar original
    Ta = copy.deepcopy(T)
    d_s = len(S)
    for i in range(d_s):
        for j in range(d_s):
            for k in range(d_s):
                if Cond_Triangular(S[i],S[j],S[k]) == 0 and Ta[i+1,j+1][k] != 0:
                    Ta[i+1,j+1][k] = 0
                    print("posicion",i,",",j," en T anulada")
                    print("T[",i+1,",",j+1,"] --->", Ta[i+1,j+1])
    return Ta

def sl2Tabla_2_TablaCompleta(S,T):
    # Ta = copy.deepcopy(T)
    TT = {}
    d = 0
    d_s = len(S)
    SS = [d]
    for s in range(d_s):
        d = d+S[s]+1
        SS.append(d)
    for i in range(d_s):
        for j in range(d_s):
            for ii in range(S[i]+1):
                for jj in range(S[j]+1):
                    TT[SS[i]+ii+1,SS[j]+jj+1] = sp.zeros(1,d)
    for i in range(d_s):
        a = S[i]
        for j in range(d_s):
            b = S[j]
            for ii in range(a+1):
                for jj in range(b+1):
                    for k in range(d_s):
                        c = S[k]
                        for kk in range(c+1):
                            Cg = CG(a/2,b/2,a/2-ii,b/2-jj,c/2,c/2-kk)
                            TT[SS[i]+ii+1,SS[j]+jj+1][SS[k]+kk] += T[i+1,j+1][k]*Cg
    return TT

def sl2Tabla_2_TablaCompleta_MG(S,T,opt_sg):
    # Ta = copy.deepcopy(T)
    TT = {}
    d = 0
    d_s = len(S)
    SS = [d]
    for s in range(d_s):
        d = d+S[s]+1
        SS.append(d)
    for i in range(d_s):
        for j in range(d_s):
            for ii in range(S[i]+1):
                for jj in range(S[j]+1):
                    TT[SS[i]+ii+1,SS[j]+jj+1] = sp.zeros(1,d)
    for i in range(d_s):
        a = S[i]
        for j in range(d_s):
            b = S[j]
            for ii in range(a+1):
                for jj in range(b+1):
                    for k in range(d_s):
                        c = S[k]
                        if opt_sg in ['y','Y','yes','Yes','YES']:
                            for kk in range(c+1):
                                Mg = MG_int_sg(a,b,c,ii,jj,kk)
                                TT[SS[i]+ii+1,SS[j]+jj+1][SS[k]+kk] += T[i+1,j+1][k]*Mg
                        elif opt_sg in ['n','N','no','No','NO']:
                            for kk in range(c+1):
                                Mg = MG_int(a,b,c,ii,jj,kk)
                                TT[SS[i]+ii+1,SS[j]+jj+1][SS[k]+kk] += T[i+1,j+1][k]*Mg
    return TT

