import numpy as np
from scipy import optimize

pi = np.pi
rows = 2
lug_strength = 276000000 # Heat treated Al-6061 (T651 Temper)
bolt_strength = 123391231
lug_density = 1234
bolt_weight = 0.13

applied_load = 5000

"""
w: width of the lug
lh: lug height
fh: flange height
t1: flange thickness
t2: lug thickness
D1: flange hole diameter
D2: fastener hole diameter
n_bolts: number of bolts
"""

variables = [w, lh, fh, t1, t2, D1, D2, n_bolts]

constraints = [{"type": "ineq", "fun": hole_spacing_check}
def opt_lug(load, material_strength):
    width, height = 

def calculate_weight(l):
    return ((l[0]*l[1] - pi*(l[6]/2)**2*l[7])*l[4] + ((l[2]-l[0]/2)*l[0]+0.5*pi*(l[0]/2)**2-pi*(l[5]/2)**2)*2*l[3])*lug_density  

def hole_spacing_check(l):
    return (w - 3 * l[6]) - 3 * l[6] * l[7]/2

def pull_through_check(l):
    
