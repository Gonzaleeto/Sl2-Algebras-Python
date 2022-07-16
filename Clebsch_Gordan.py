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

##############################################################################

def Limpiar(L):
    LL = copy.deepcopy(L)
    dim = len(L)
    k = -1
    if LL != []:
        for i in range(dim):
            if LL[i][2] == 0:
                k = i
                break
        if k>-1:
            LL.remove(LL[k])
            LL = Limpiar(LL)
    return LL

def Colapsar(L):
    Pesos = []
    k = -1
    dim = len(L)
    for i in range(dim):
        pair = [L[i][0],L[i][1]]
        if pair in Pesos:
            kp = Pesos.index(pair)
            k = i
            break
        Pesos.append(pair)
    if k == -1:
       LL = copy.deepcopy(L)
    else:
       pair_aux = [L[k][0],L[k][1],sp.simplify(sp.factor(L[kp][2]+L[k][2]))]
       LL = copy.deepcopy(L)
       LL[kp] = pair_aux
       LL.remove(LL[k])
       LL = Colapsar(LL)
    return LL

def E_simet(a,b,L):
    dim = len(L)
    R = []
    for i in range(dim):
        if L[i][0]+2 <= a:
            pair = [L[i][0]+2,L[i][1],(a+L[i][0]+2)/2*(a-L[i][0])/2*L[i][2]]
            R.append(pair)
        if L[i][1]+2 <= b:
            pair = [L[i][0],L[i][1]+2,(b+L[i][1]+2)/2*(b-L[i][1])/2*L[i][2]]
            R.append(pair)
    return Limpiar(Colapsar(R))

def E_CG(a,b,L):
    dim = len(L)
    R = []
    for i in range(dim):
        if L[i][0]+2 <= a:
            pair = [L[i][0]+2,L[i][1],sp.sqrt(a/2*(a/2+1)-L[i][0]/2*(L[i][0]/2+1))*L[i][2]]
            R.append(pair)
        if L[i][1]+2 <= b:
            pair = [L[i][0],L[i][1]+2,sp.sqrt(b/2*(b/2+1)-L[i][1]/2*(L[i][1]/2+1))*L[i][2]]
            R.append(pair)
    return Limpiar(Colapsar(R)) 

def F_simet(a,b,L):
    dim = len(L)
    R = []
    for i in range(dim):
        if L[i][0]-2 >= -a:
            pair = [L[i][0]-2,L[i][1],L[i][2]]
            R.append(pair)
        if L[i][1]-2 >= -b:
            pair = [L[i][0],L[i][1]-2,L[i][2]]
            R.append(pair)
    return Limpiar(Colapsar(R))       
       
def F_CG(a,b,L):
    dim = len(L)
    R = []
    for i in range(dim):
        if L[i][0]-2 >= -a:
            pair = [L[i][0]-2,L[i][1],sp.sqrt(a/2*(a/2+1)-L[i][0]/2*(L[i][0]/2-1))*L[i][2]]
            R.append(pair)
        if L[i][1]-2 >= -b:
            pair = [L[i][0],L[i][1]-2,sp.sqrt(b/2*(b/2+1)-L[i][1]/2*(L[i][1]/2-1))*L[i][2]]
            R.append(pair)
    return Limpiar(Colapsar(R))

# Coeficiente de Clebsch-Gordan estandares

def CG_0(j1,j2,m1,m2,j,m):
    R = 0
    if m == m1+m2 and j+j1-j2 >= 0 and j-j1+j2 >= 0 and -j+j1+j2 >= 0:
        if abs(m1) <= j1 and abs(m2) <= j2 and abs(m) <= j:
            # Aumentar uno mas el paso en el for por cuestiones de Python 
            for k in range(mt.floor(j1+j2-j)+1):
                # print("k=", k," j1+j-j2 = ",mt.floor(j1+j2-j))
                if k <= j1-m1 and k <= j2+m2 and k >= j2-j-m1 and k >= j1-j+m2:
                    t1 = sp.factorial(j1+j2-j-k)
                    t2 = sp.factorial(j1-m1-k)
                    t3 = sp.factorial(j2+m2-k)
                    t4 = sp.factorial(k-j2+j+m1)
                    t5 = sp.factorial(k-j1+j-m2)
                    x = sp.factorial(k)*t1*t2*t3*t4*t5
                    R = R + (-1)**(k)/x
            s1 = 2*j+1
            s2 = sp.factorial(j+j1-j2)
            s3 = sp.factorial(j-j1+j2)
            s4 = sp.factorial(-j+j1+j2)
            s5 = sp.factorial(j+j1+j2+1)
            
            R = sp.sqrt(s1*s2*s3*s4/s5)*R
            
            s1 = sp.factorial(j+m)
            s2 = sp.factorial(j-m)
            s3 = sp.factorial(j1+m1)
            s4 = sp.factorial(j1-m1)
            s5 = sp.factorial(j2+m2)
            s6 = sp.factorial(j2-m2)
            
            R = sp.sqrt(s1*s2*s3*s4*s5*s6)*R
    return R

def CG(j1,j2,m1,m2,j,m):
    if m >= 0:
        if j1 >= j2:
            R = CG_0(j1,j2,m1,m2,j,m)
        else:
            R = (-1)**(j-j1-j2)*CG_0(j2,j1,m2,m1,j,m)
    else:
        if j1 >= j2:
            R = (-1)**(j-j1-j2)*CG_0(j1,j2,-m1,-m2,j,-m)
        else:
            R = CG_0(j2,j1,-m2,-m1,j,-m)
    return R

# Coeficientes de Clebsch-Gordan enteros

def CG_int(n,m,mu,i,j,p):
    R = 0
    if sp.floor(i) == i and sp.floor(j) == j and sp.floor(p) == p:
        if 0 <= i <= n and 0 <= j <= m and 0 <= p <= mu:
            if Cond_Triangular(n, m, mu) != 0:
                k = xx(n,m,mu)
                if n+m-2*(i+j) == mu-2*p:
                    # Aumentar uno mas el paso en el for por cuestiones de Python
                    for r in range(k+1):
                        b1 = sp.binomial(i+j-k,i-r)
                        b2 = sp.binomial(n-r,k-r)
                        b3 = sp.binomial(m+r-k,r)
                        R = R + (-1)**(r+1)*b1*b2*b3
    return R

def CG_int_sg(n,m,mu,i,j,p):
    R = 0
    if sp.floor(i) == i and sp.floor(j) == j and sp.floor(p) == p:
        if 0 <= i <= n and 0 <= j <= m and 0 <= p <= mu:
            if Cond_Triangular(n, m, mu) != 0:
                k = xx(n,m,mu)
                sg = k % 2
                if n+m-2*(i+j) == mu-2*p:
                    # Aumentar uno mas el paso en el for por cuestiones de Python
                    for r in range(k+1):
                        b1 = sp.binomial(i+j-k,i-r)
                        b2 = sp.binomial(n-r,k-r)
                        b3 = sp.binomial(m+r-k,r)
                        R = R + (-1)**(r+1)*b1*b2*b3
                    R = (-1)**(sg+1)*R
    return R

# Coeficientes multiplicadores V(n)xV(m) --> V(s)

def DD_int(n,m,mu):
    k = xx(n,m,mu)
    R = sp.binomial(n+m-k+1,k)*sp.binomial(mu,n-k)
    return R

def MG_int(n,m,mu,i,j,p):
    mg_aux = 0
    if 0 <= p <= mu:
        k = xx(n,m,mu)
        cg_aux = CG_int(n,m,mu,n-i,m-j,mu-p)
        dd_aux = DD_int(n,m,mu)
        if cg_aux != 0 and dd_aux != 0:
            mg_aux = (-1)**k*cg_aux/dd_aux
    else:
        print("valor de p= ",p," fuera de rango")
    return mg_aux

def MG_int_sg(n,m,mu,i,j,p):
    mg_aux = 0
    if 0 <= p <= mu:
        k = xx(n,m,mu)
        cg_aux = CG_int_sg(n,m,mu,n-i,m-j,mu-p)
        dd_aux = DD_int(n,m,mu)
        if cg_aux != 0 and dd_aux != 0:
            mg_aux = (-1)**k*cg_aux/dd_aux
    else:
        print("valor de p= ",p," fuera de rango")
    return mg_aux
    
def show_coef_CG(n,m,mu,opt_CG,opt_sg):
    if opt_CG in ['std', 'Std', 'STD']:
        for i in range(n+1):
            for j in range(m+1):
                for p in range(mu+1):
                    if n+m-2*(i+j) == mu-2*p:
                        print("CG(",n,m,mu,",",i,j,p,")")
                        print (mu,",",mu-2*p," --> ",CG(n/2,m/2,n/2-i,m/2-j,mu/2,mu/2-p))
    elif opt_CG in ['int','Int','INT']:
        for i in range(n+1):
            for j in range(m+1):
                for p in range(mu+1):
                    if n+m-2*(i+j) == mu-2*p:
                        print("BG(",n,m,mu,",",i,j,p,")")
                        if opt_sg in ['y','Y','yes','Yes','YES']:
                            print(mu,",",mu-2*p," --> ",CG_int_sg(n,m,mu,i,j,p))
                        elif opt_sg in ['n','N','no','No','NO']:
                            print(mu,",",mu-2*p," --> ",CG_int(n,m,mu,i,j,p))
                        else:
                            print("Opciones de CG enteros incorrecta")
                            break
    else:
        print("opciones de CG incorrectas")
    return print("Fin")

def show_coef_MG(n,m,mu,opt_sg):
    for i in range(n+1):
        for j in range(m+1):
            for p in range(mu+1):
                if n+m-2*(i+j) == mu-2*p:
                    print("MG(",n,m,mu,",",i,j,p,")")
                    if opt_sg in ['y','Y','yes','Yes','YES']:
                        print(mu,",",mu-2*p," --> ",MG_int_sg(n,m,mu,i,j,p))
                    elif opt_sg in ['n','N','no','No','NO']:
                        print(mu,",",mu-2*p," --> ",MG_int(n,m,mu,i,j,p))
                    else:
                        print("Opciones de CG enteros incorrecta")
                        break
    return print("Fin")

    
            
            
            
            
            
            
            
            
            