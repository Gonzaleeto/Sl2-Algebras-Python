import numpy as np
import sympy as sp
import math  as mt
import copy # Para usar deepcopy en algunos objetos

from sympy import MatrixSymbol, Identity

# Para imprimir en pantallas mas legible los resultados

from sympy import init_printing
init_printing(use_latex=True)

##############################################################################

def xx(a,b,m):
    return mt.floor((a+b-m)/2)

def lambda_int(m,n,mu,r):
    R = 0
    if 0 <= r and r <= xx(m, n, mu):
        R = (-1)**r*mt.comb(m-r, xx(m,n,mu)-r)*mt.comb(n+r-xx(m,n, mu), r)
    return(R)

def Cond_Triangular (a, b, c):
    r = 1
    if mt.floor(a/2 + b/2 + c/2) != a/2 + b/2 + c/2:
        r = 0
    if c < abs(a-b):
        r = 0
    if a+b < c:
        r = 0
    return(r)

def LieMat(A:sp.Matrix, B:sp.Matrix):
    t1 = A * B
    t2 = B * A
    return(t1 - t2)

def AntiLieMat(A:sp.Matrix, B:sp.Matrix):
    t1 = A * B
    t2 = B * A
    return(t1 + t2)

def Multiplica_a(u,v,T):
    dim = len(u) # Tamaño de u
    w = sp.MutableDenseNDimArray.zeros(dim)
    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                w[[k]] = w[[k]] + u[i]*v[j]*T[i+1,j+1][k]
                
    return w

def Multiplica_v(u,v,T):
    dim = len(u) # Tamaño de u
    w = sp.zeros(1,dim)
    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                w[k] = w[k] + u[i]*v[j]*T[i+1,j+1][k]
                
    return w

def Verif_Abel(T):
    dim = len(T[1,2])
    cond = set() # Conjunto vacio
    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                s = sp.simplify(T[i+1,j+1][k]-T[j+1,i+1][k])
                if s != 0:
                    cond.add(s)
    return cond

def Verif_Anti_Abel(T):
    dim = len(T[1,2])
    cond = set() # Conjunto vacio
    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                s = sp.simplify(T[i+1,j+1][k]+T[j+1,i+1][k])
                if s != 0:
                    cond.add(s)
    return cond

def Jac(i,j,k,T):
    dim = len(T[1,2])
    cond = set()
    for l in range(dim):
        s = 0
        for r in range(dim):
            # partimos la suma solo por cuestion de legibilidad 
            s = s + T[i+1,j+1][r]*T[r+1,k+1][l] 
            s = s - T[i+1,k+1][r]*T[r+1,j+1][l] 
            s = s + T[j+1,k+1][r]*T[r+1,i+1][l]
        s = sp.simplify(s)
        if s != 0:
            cond.add(s)
    return cond

def Verif_Lie(T):
    dim = len(T[1,2])
    cond = set() # Conjunto vacio
    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                # Anticonmutatividad
                s1 = sp.simplify(T[i+1,j+1][k]+T[j+1,i+1][k])
                if s1 != 0:
                    cond.add(s1)
                # Jacobi
                cond_aux = Jac(i,j,k,T)
                cond.update(cond_aux)
                
    return cond

def Verif_Jac(T):
    dim = len(T[1,2])
    cond = set() # Conjunto vacio
    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                # Jacobi
                cond_aux = Jac(i,j,k,T)
                cond.update(cond_aux)
                
    return cond

def Asoc(i,j,k,T):
    dim = len(T[1,2])
    s = 0
    cond = set()
    for l in range(dim):
        for r in range(dim):
            # partimos la suma solo por cuestion de legibilidad 
            s = s + T[i+1,j+1][r]*T[r+1,k+1][l] 
            s = s - T[j+1,k+1][r]*T[i+1,r+1][l]
        s = sp.simplify(s)
        if s != 0:
            cond.add(s)
    return cond

def Asoc_pos(i,j,k,T):
    dim = len(T[1,2])
    s = 0
    cond = set()
    for l in range(dim):
        for r in range(dim):
            # partimos la suma solo por cuestion de legibilidad 
            s = s + T[i+1,j+1][r]*T[r+1,k+1][l] 
            s = s - T[j+1,k+1][r]*T[i+1,r+1][l]
        s = sp.simplify(s)
        if s != 0:
            tupla = sp.Array([i+1,j+1,k+1,s])
            cond.add(tupla)
    return cond

def Verif_Asoc(T):
    dim = len(T[1,2])
    cond = set() # Conjunto vacio
    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                # Asociatividad
                cond_aux = Asoc(i,j,k,T)
                cond.update(cond_aux)
    return cond

def Verif_Asoc_on_pos(T):
    dim = len(T[1,2])
    cond = set() # Conjunto vacio
    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                # Asociatividad
                cond_aux = Asoc_pos(i,j,k,T)
                cond.update(cond_aux)
    return cond

def Verif_LeftAlt_Law(T):
    dim = len(T[1,2])
    cond = set()
    for i in range(dim):
        for j in range(dim):
            for l in range(dim):
                sum = 0
                for r in range(dim):
                    # partimos la suma solo por cuestion de legibilidad 
                    sum = sum + T[i+1,i+1][r]*T[r+1,j+1][l] 
                    sum = sum - T[i+1,j+1][r]*T[i+1,r+1][l]
                sum = sp.simplify(sum)
                if sum != 0:
                    tupla = sp.Array([i+1,i+1,j+1,sum])
                    cond.add(tupla)
    return cond

def Verif_RightAlt_Law(T):
    dim = len(T[1,2])
    cond = set()
    for i in range(dim):
        for j in range(dim):
            for l in range(dim):
                sum = 0
                for r in range(dim):
                    # partimos la suma solo por cuestion de legibilidad 
                    sum = sum + T[i+1,j+1][r]*T[r+1,j+1][l] 
                    sum = sum - T[j+1,j+1][r]*T[i+1,r+1][l]
                sum = sp.simplify(sum)
                if sum != 0:
                    tupla = sp.Array([i+1,i+1,j+1,sum])
                    cond.add(tupla)
    return cond

def Malcev(i,j,k,T):
    dim = len(T[1,2])
    cond = set()
    for t in range(dim):
        sum1 = 0
        sum2 = 0
        for s in range(dim):
            for r in range(dim):
                sum1 = sum1 + T[i+1,j+1][r]*T[r+1,k+1][s]*T[s+1,i+1][t]
                sum1 = sum1 + T[j+1,k+1][r]*T[r+1,i+1][s]*T[s+1,i+1][t]
                sum1 = sum1 + T[k+1,i+1][r]*T[r+1,j+1][s]*T[s+1,i+1][t]
                
                sum2 = sum2 + T[i+1,j+1][r]*T[i+1,k+1][s]*T[r+1,s+1][t]
                sum2 = sum2 + T[i+1,k+1][r]*T[j+1,r+1][s]*T[s+1,i+1][t]
                sum2 = sum2 + T[i+1,k+1][r]*T[r+1,i+1][s]*T[s+1,j+1][t]
        ss = sp.simplify(sum1-sum2)
        if ss != 0:
            tupla = sp.Array([i+1,j+1,k+1,ss])
            cond.add(tupla)
    return cond

def Verif_Malcev(T):
    dim = len(T[1,2])
    cond = set()
    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                cond_aux = Malcev(i,j,k,T)
                cond.update(cond_aux)
    return cond

def Derivaciones_T(T,z):
    n = len(T[1,2])
    var = sp.Array([sp.Symbol("x{}{}".format(i+1,j+1)) for i in range(n) for j in range(n)])
    X = sp.Matrix(n,n,var)
    ecu = []
    for i in range(n):
        for j in range(n):
            Tij = sp.Matrix(n,1,T[i+1,j+1])
            lizq = X*Tij
            vi = sp.zeros(1,n)
            vi[i] = 1
            vj = sp.zeros(1,n)
            vj[j] = 1
            lder = Multiplica_v(X[:,i],vj,T) + Multiplica_v(vi,X[:,j],T)
            lder = sp.Matrix(n,1,lder)
            expr = lizq-lder
            for k in range(len(expr)):
                ecu.append(expr[k])
    M, b = sp.linear_eq_to_matrix(ecu, var)
    Sol = M.nullspace()
    Dsol = len(Sol)
    print("Dimension Der = ", Dsol)
    
    vz = sp.Array([sp.Symbol("z{}".format(i+1)) for i in range(Dsol)])
    M = sp.zeros(n,n)
    for i in range(Dsol):
        M = M + vz[i]*sp.Matrix(n,n,Sol[i])
    return M

"""Verificamos la compatibilidad entre un producto asociativo y un corchete
   de Lie
   Ta = Tabla de multiplicar de un producto asociativo
   Tl = Tabla de multiplicar de un corchete de Lie
"""
def Verif_comp_AsocLie(Ta,Tl):
    dim = len(Ta[1,2])
    cond = set()
    if dim == len(Tl[1,2]):
        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    for s in range(dim):
                        suma = 0
                        for r in range(dim):
                            suma = suma + Ta[j+1,k+1][r]*Tl[i+1,r+1][s] 
                            suma = suma - Tl[i+1,j+1][r]*Ta[r+1,k+1][s] 
                            suma = suma - Tl[i+1,k+1][r]*Ta[j+1,r+1][s]
                        suma = sp.simplify(suma)
                        if suma != 0:
                            cond.add(suma)
    else:
        print("Error en dimensiones de Tablas",dim,len(Tl[1,2]))
    return cond
                        
                    
                    
                    
                    
                    