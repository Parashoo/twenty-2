import numpy as np

T_MAIN = 425
T_RCT = 25
MIN_MASS = 845.55
MASS_SOLAR_PANEL = 34.62
DISTANCE_CG = 2.37348

a_translation = T_RCT / MIN_MASS
a_main_engine = T_MAIN/ MIN_MASS

ax = a_translation
ay = -a_main_engine # y is downward in the coordinate system
az = a_translation

a = np.array([ax, ay, az])
reaction_force = MASS_SOLAR_PANEL * a
r = np.array([0, 0, DISTANCE_CG])
reaction_moment = np.cross(reaction_force, r)
