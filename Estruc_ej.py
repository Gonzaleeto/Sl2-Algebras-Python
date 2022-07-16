# from __future__ import division
import numpy as np
import sympy as sp
import math  as mt
import copy # Para usar deepcopy en algunos objetos

from sympy import MatrixSymbol, Identity

# Para imprimir en pantallas mas legible los resultados

from sympy import init_printing
init_printing(use_latex=True)

##############################################################################

# Variables simbolicas

x = sp.symbols('x')
q = sp.symbols('q')
z = sp.symbols('z')

# Tablas de multiplicar

S1 = [2]

T1 = dict({(1, 1) : [q]})

T2 = dict({(1, 2) : [0, 0], (2, 2) : [0, 0], (2, 1) : [0, 0], 
           (1, 1) : [q, 0]})

TT1 = dict({(1, 1) : [x, 0, 3], (1, 2) : [0, 1, 0], (3, 3) : [0, 0, 0]})

TT2 = dict({(3, 1) : [0, (1/4)*q, 0], (1, 3) : [0, -(1/4)*q, 0], 
            (3, 2) : [0, 0, (1/4)*q], (3, 3) : [0, 0, 0], (1, 1) : [0, 0, 0], 
            (2, 2) : [0, 0, 0], (2, 1) : [(1/4)*q, 0, 0], 
            (1, 2) : [-(1/4)*q, 0, 0], (2, 3) : [0, 0, -(1/4)*q]})

TT3 = dict({(4, 4) : [0, 0, 0, 0], (3, 2) : [0, 0, (1/4)*q, 0], 
            (3, 1) : [0, (1/4)*q, 0, 0], (3, 3) : [0, 0, 0, 0], 
            (2, 3) : [0, 0, -(1/4)*q, 0], (1, 3) : [0, -(1/4)*q, 0, 0], 
            (1, 2) : [-(1/4)*q, 0, 0, 0], (2, 2) : [0, 0, 0, 0], 
            (3, 4) : [0, 0, 0, 0], (4, 3) : [0, 0, 0, 0], 
            (4, 2) : [0, 0, 0, 0], (1, 4) : [0, 0, 0, 0], (2, 4) : [0, 0, 0, 0], 
            (4, 1) : [0, 0, 0, 0], (2, 1) : [(1/4)*q, 0, 0, 0], 
            (1, 1) : [0, 0, 0, 0]})

TT4 = dict({(4, 5) : [0, 0, 0, 0, 0], (3, 1) : [(1/21)*q, 0, 0, 0, 0], 
              (1, 2) : [0, 0, 0, 0, 0], (3, 3) : [0, 0, -(1/21)*q, 0, 0], 
              (1, 4) : [0, (1/14)*q, 0, 0, 0], (1, 1) : [0, 0, 0, 0, 0], 
              (1, 3) : [(1/21)*q, 0, 0, 0, 0], (4, 2) : [0, 0, (1/42)*q, 0, 0], 
              (3, 2) : [0, -(1/42)*q, 0, 0, 0], (3, 4) : [0, 0, 0, -(1/42)*q, 0], 
              (2, 2) : [-(1/21)*q, 0, 0, 0, 0], 
              (4, 1) : [0, (1/14)*q, 0, 0, 0], (2, 1) : [0, 0, 0, 0, 0], 
              (1, 5) : [0, 0, (1/21)*q, 0, 0], (4, 4) : [0, 0, 0, 0, -(1/14)*q], 
              (5, 3) : [0, 0, 0, 0, (1/21)*q], (5, 1) : [0, 0, (1/21)*q, 0, 0], 
              (2, 5) : [0, 0, 0, (1/21)*q, 0], (5, 5) : [0, 0, 0, 0, 0], 
              (2, 3) : [0, -(1/42)*q, 0, 0, 0], (5, 4) : [0, 0, 0, 0, 0], 
              (5, 2) : [0, 0, 0, (1/21)*q, 0], (4, 3) : [0, 0, 0, -(1/42)*q, 0], 
              (2, 4) : [0, 0, (1/42)*q, 0, 0], (3, 5) : [0, 0, 0, 0, (1/21)*q]})

###############################################################################

q1 = sp.symbols('q1')
q2 = sp.symbols('q2')
q3 = sp.symbols('q3')
q4 = sp.symbols('q4')
q5 = sp.symbols('q5')
q6 = sp.symbols('q6')
q7 = sp.symbols('q7')
q8 = sp.symbols('q8')

S3 = [1,2]

T3 = dict({(1, 1) : [q1, q2], (1, 2) : [q3, q4], (2, 1) : [q5, q6], 
           (2, 2) : [q7, q8]})

###############################################################################

""" Si se usa esta forma, entonces la dimension del Array es inmutable.
    Otra punto negativo de lo anterior es que al usar numpy, el array
    creado tienen elementos en punto flotantes"""

V1 = sp.Array(np.zeros(6))

V1_a = sp.Array.zeros(6)

# Con esta forma se crea un Array con dimension mutable

V2 = sp.MutableDenseNDimArray.zeros(6)

# Las posiciones deben ser tuplas, por ejemplo

V2[[2]] = x
V2[[3]] = 5


# Formas de crear una matriz simbolica 

n_f = 3
n_c = 3

M_ej1 = sp.Array([[sp.Symbol("x_{},{}".format(i,j)) for j in range(n_c)] for i in range(n_f)])
M_ej2 = sp.MatrixSymbol('m',n_f,n_c)

# Forma de crear conjuntos de variables

nvar = 3

v_ej1 = sp.symbols('v0:%d'%nvar)
v_ej2 = sp.Array([sp.Symbol("x{}".format(i)) for i in range(nvar)])
v_ej3 = sp.Array([sp.Symbol("x{}{}".format(i,j)) for i in range(nvar) for j in range(nvar)])

M_ej3 = sp.Matrix(nvar,nvar,v_ej3)

# exp = Multiplica_v(M_ej3[:,1],v,TT2)

R0 = [[1,1,1]]

R1 = [[2, 2, 2], [2, 2, -2],[3,2,0]]

R2 = [[-1, -1, 1], [-1, -1, 1]]