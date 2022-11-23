from math import *

def phi(delta_a, delta_b):
    return delta_a / (delta_a + delta_b)

def delta_a(t, E_a, D_fo, D_fi):
    return 4 * t / (E_a * pi * (D_fo**2 - D_fi**2))

def delta_b(E_b, L_i, A_i):
    sum_factor = 0
    for i in range(len(L_i)):
        sum_factor += L_i[i] / A_i[i]
    return 1 / E_b * sum_factor