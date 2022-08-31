# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 16:18:07 2022
@author: Magrini, L. & Gadotti, M.
"""

import sympy as sym
import math

# Definicao dos parametros utilizados para as simulacoes
global alpha, beta, p, gamma_a, gamma_s, gamma_as, mi, v, epsilon, tau, N
alpha, beta, p, gamma_a, gamma_s, gamma_as, mi, v, epsilon, tau, N = \
    [0.8, 0.7662, 0.5, 2/7, 1/7, 1/14, 0.0723, 2/100, 0.8, 4, 10**6]


# Variaveis do Sistema
S, V, E, A, I, R  = sym.symbols('S, V, E, A, I, R')

# Definindo a natureza dos parâmetros
N = sym.Symbol("N", positive = True)
mi = sym.Symbol("mi", positive = True)
epsilon = sym.Symbol("epsilon", positive = True)
alpha = sym.Symbol("alpha", positive = True)
gamma_as = sym.Symbol("gamma_as", positive = True)
gamma_s = sym.Symbol("gamma_s", positive = True)
gamma_a = sym.Symbol("gamma_a", positive = True)
lamb = sym.Symbol("lamb", positive = True)
tau = sym.Symbol("tau", positive = True)
beta = sym.Symbol("beta", positive = True)
exp = sym.Symbol("exp", positive = True)
p = sym.Symbol("p", positive = True)
x = sym.Symbol("x")

# X :: subsistema infectado 
X = sym.Matrix([beta*((A+I)/N)*S + beta*(1-epsilon)*((A+I)/N)*V - beta*sym.exp(-mi*tau)*((A+I)/N)*S - mi*E, alpha*beta*sym.exp(-mi*tau)*((A+I)/N)*S - (gamma_a+gamma_as+mi)*A, beta*(1-alpha)*sym.exp(-mi*tau)*((A+I)/N)*S - (gamma_s+mi)*I + gamma_as*A])
Y = sym.Matrix([E, A, I])
jacobiana = X.jacobian(Y)
J = jacobiana.det()

# Delta : Matriz de Transicao
Delta = sym.Matrix([[-mi, 0, 0], [0, -gamma_a-gamma_as-mi, 0], [0, gamma_as, -gamma_s-mi]])

# M = J - Delta
M = jacobiana - Delta

# Inversa da Matriz de Transicao
Delta_Inv = Delta.inv()

# Matriz da Proxima Geracao
MNG = -M*Delta_Inv

# R0 (Traco da MNG)
R0 = MNG.trace()

meuR0 = ((S*beta*sym.exp(-mi*tau))/N)*((alpha*gamma_as+alpha*(gamma_s+mi)+(1-alpha)*(gamma_a+gamma_as+mi)))/((gamma_a+gamma_as+mi)*(gamma_s+mi))

sym.simplify(R0-meuR0)==0

