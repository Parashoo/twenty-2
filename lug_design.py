import math
import numpy as np

#forces = [float(x) for x in input("Give the forces (x y z):\n").split()]
#moment = float(input("Is there a moment? (give the moment)\n"))
#distance = float(input("Give the distance between lugs\n"))
#uts = float(input("Ultimate tensile strength of material\n"))
#F_ty = float(input("Yield strength of material\n"))

forces = [747.9, 2243.74, 7474.9]
moment = 184
distance = 0.1
F_tu = 420e6
F_ty = 350e6

F_1 = moment / distance
P_1 = math.sqrt((forces[1] / 2 + F_1)**2 + (forces[2] / 2)**2)
P_2 = math.sqrt((forces[1] / 2 - F_1)**2 + (forces[2] / 2)**2)

K_t = 0.8

def A_t(D, w, t):
    A_t = (w - D) * t
    return A_t

def A_br(D, t):
    A_br = D * t
    return A_br
 
def A_av(D, w, t):
    A1 = w / 2 - math.sin(math.radians(45)) * D / 2
    A2 = t * (w - D) / 2
    A_av = 6 / (4/A1 + 2/A2)
    return A_av

def K_ty(A_av, A_br):
    x = A_av / A_br
    K_ty = -0.3499 * x**5 + 1.4642 * x**4 - 2.1259 * x**3 + 0.8982 * x**2 + 1.1588 * x - 0.012
    return K_ty

def P_u(D, w, t):
    P_u = K_t * F_tu * A_t(D, w, t)
    return P_u

def P_ty(D, w, t):
    P_ty = K_ty(A_av(D, w, t), A_br(D, t)) * A_br(D, t) * F_ty
    return P_ty

def R_a(D, w, t):
    R_a = forces[0] / P_u(D, w, t)
    return R_a

def R_tr(D, w, t):
    R_tr = forces[1] / P_ty(D, w, t)
    return R_tr

def sum_of_ratios(D, w, t):
    sum_of_ratios = R_a(D, w, t)**1.6 + R_tr(D, w, t)**1.6
    return sum_of_ratios

def margin_of_safety(D, w, t):
    ms = 1 / (sum_of_ratios(D, w, t)**0.625) - 1
    return ms

D_arr = np.linspace(0.001, 0.021, 10)
w_arr = np.linspace(0.001, 0.021, 10)
t_arr = np.linspace(0.001, 0.006, 10)

for D in D_arr:
    for w in w_arr:
        if w > D:
            for t in t_arr:
                ms = margin_of_safety(D, w, t)
                if 1.5 < ms < 1.6:
                    print(f"D: {D:.3f}, w: {w:.3f}, t:{t:.3f} -> ms: {ms:.2f}")
