import math

forces = [int(x) for x in input("Give the forces (x y z):\n").split()]
moment = int(input("Is there a moment? (give the moment)\n"))
distance = int(input("Give the distance between lugs\n"))

F_1 = moment / distance
P_1 = math.sqrt((forces[1] / 2 + F_1)**2 + (forces[2] / 2)**2)
P_2 = math.sqrt((forces[1] / 2 - F_1)**2 + (forces[2] / 2)**2)

K_t = 0.8

def A_av():
    return 1

def K_ty(A_av, A_br):
    A_ratio = A_av / A_br
    if A_ratio <= 0.6:
        K_ty = 26/20 * A_ratio
    else:
        K_ty = 30/44 * A_ratio + 26/20 * 0.6
    return K_ty
