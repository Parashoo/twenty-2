import math

forces = [int(x) for x in input("Give the forces (x y z):\n").split()]
moment = int(input("Is there a moment? (give the moment)\n"))
distance = int(input("Give the distance between lugs (for single lug input 0)\n"))

F_1 = moment / distance
P_1 = math.sqrt((forces[1] + F_1)**2 + forces[2]**2)
P_2 = math.sqrt((forces[1] - F_1)**2 + forces[2]**2)


