import math
import numpy as np

forces = [float(x) for x in input("Give the forces (x y z):\n").split()]
moment = float(input("Is there a moment? (give the moment)\n"))
distance = float(input("Give the distance between lugs\n"))
uts = float(input("Ultimate tensile strength of material\n"))
F_ty = float(input("Yield strength of material\n"))

F_1 = moment / distance
P_1 = math.sqrt((forces[1] / 2 + F_1)**2 + (forces[2] / 2)**2)
P_2 = math.sqrt((forces[1] / 2 - F_1)**2 + (forces[2] / 2)**2)

K_t = 0.8

def A_av(D, w, t):
    A1 = w / 2 - math.sin(math.radians(45)) * D / 2
    A2 = t * (w - D) / 2
    A_av = 6 / (4/A1 + 2/A2)
    return A_av

def A_t(D, w, t):
    return (w - D) * t

def A_br(D, t):
    return D * t

def K_ty(A_av, A_br):
    x = A_av / A_br
    K_ty = -0.3499 * x**5 + 1.4642 * x**4 - 2.1259 * x**3 + 0.8982 * x**2 + 1.1588 * x - 0.012
    return K_ty

def P_u(D, w, t):
    return K_t * uts * A_t(D, w, t)

def P_ty(D, w, t):
    return K_ty(A_av(D, w,t), A_br(D, t)) * A_br(D, t) * F_ty

def R_a(D, w, t):
    return forces[0] / P_u(D, w, t)

def R_tr(D, w, t):
    return forces[1] / P_ty(D, w, t)

def margin_of_safety(D, w, t):
    ms = 1 / np.power((np.power(R_a(D, w, t), 1.6) + np.power(R_tr(D, w, t), 1.6)), 0.625) - 1
    return ms

D_arr = np.linspace(0, 0.1, 5)
w_arr = np.linspace(0, 0.1, 5)
t_arr = np.linspace(0, 0.02, 5)

for D in D_arr:
    for w in w_arr:
        if w > D:
            for t in t_arr:
                ms = margin_of_safety(D, w, t)
                print(f"D: {D}, w: {w}, t:{t} -> ms: {ms}")
