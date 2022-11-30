"""
This code analyzes the minimum bolt size for fasteners in a joint. Here there are 2 
classes, one is the fastener, the other is the joint. A specific amount of fasteners
can be added to the joint at different positions in a normal carthesian coordinate system.
Aditonally, the force and moment can be added to the joint, and the margins of safety
can be computed. Also the bolt diameters can be optimized with a specific safety factor 
in mind.

Rules:
    - A fitting factor of at least 1.15 shall be used (Niu page 274)
    - The efficiency of the joint should be greater than that of the structure (Niu page 275)
    - Fasteners must not have different sizes (Niu page 276)
    - A safety factor of 1.5 should be used (Niu page 277)
    - The spacing of fasteners should be at least 4 times their diameter (Niu page 278)
    - The fasteners should be put at least 2 times their diameter of the edge (Niu page 278)
    
    - Protuding head rivets will be used in this anaylsis
"""
import numpy as np


class Fastener():
    def __init__(self, diameter, outer_diamter, yield_stress, x, y):
        self.diameter = diameter
        self.outer_diameter = outer_diamter
        self.area = np.pi * diameter**2 / 4
        self.position = np.array([x, y, 0])
        self.yield_stress = yield_stress
        
        # Niu page 281
        # If you compare with the table of page 287 for diameter = 1/8 inch
        # and yield strength double of the shear strength they give, this is
        # a conservative estimate (531 lbs they vs 502 lbs formula after
        # conversions
        
        self.max_shear_load = 0.5*yield_stress * (np.pi*diameter**2)/4
    
    def update(self):
        self.area = np.pi * self.diameter**2 / 4
        self.max_shear_load = 0.5*self.yield_stress * np.pi*self.diameter**2/4

class Joint():
    def __init__(self, thickness, sigma_y):
        self.fasteners = []
        self.force = np.zeros(3)
        self.moment = np.zeros(3)
        self.thickness = thickness
        self.sigma_y = sigma_y
        
        # Printing is disabled during optimization
        self.printing = True
        
        
    def add_fastener(self, fastener):
        self.fasteners.append(fastener)
    
    def update_centroid(self):
        # sum(A*d**2) needed for torque shear, scalar
        self.total_SMOI = 0
        # sum(A*d) vectorially for centroid determination
        self.total_FMOI = np.zeros(3)
        
        self.total_area = 0
        for fastener in self.fasteners:
            # sum(A*d)
            self.total_FMOI += fastener.position * fastener.area
            self.total_area += fastener.area
        self.centroid = self.total_FMOI * self.total_area**(-1)
        
        for fastener in self.fasteners:
            distance_centroid = np.linalg.norm(self.centroid - fastener.position)
            self.total_SMOI += distance_centroid**2 * fastener.area
        
    
    def add_load(self, force, moment):
        self.force = force
        self.moment = moment
    
    def calculate_ms(self):
        torque = self.moment[2]
        minimum_ms = 1e99 #some arbitrary high number
        for n,fastener in enumerate(self.fasteners):
            # Calculating the shear load on the fastener due to shear
            shear_load_shear = self.force*(fastener.area/self.total_area)
            shear_load_shear[2] = 0 # only looking at the planar case
            
            # Calculating the shear load on the fastener due to torque
            centroid_distance = fastener.position - self.centroid 
            r = np.linalg.norm(centroid_distance)
            
            # Getting the torque shear force direction
            shear_torque_magnitude = torque * fastener.area * r / self.total_SMOI
            shear_torque_direction = np.cross(centroid_distance, np.array([0,0,torque]))
            # Getting the torque shear force magnitutde
            direction_mag = np.linalg.norm(shear_torque_direction)
            if direction_mag != 0:
                shear_torque_direction /= direction_mag
            # Getting the torque shear force total vector
            shear_load_torque = shear_torque_magnitude*shear_torque_direction
            
            # Getting total shear load by combining both the shear due to shear
            # and the shear due to torque
            total_shear_load = shear_load_torque + shear_load_shear
            shear_magnitude = np.linalg.norm(total_shear_load)
            
            # bearing check
            sigma_br = shear_magnitude / (fastener.diameter*self.thickness)
            
            # First get the force in and out of the fasteners
            normal_force = self.force[2] * fastener.area/self.total_area
            # Force due to the y moment
            forcex = -self.moment[1] * fastener.area * centroid_distance[0] / self.total_SMOI
            # Force due to the x moment
            forcey = self.moment[0] * fastener.area * centroid_distance[1] / self.total_SMOI
            # Total force
            total_normal_force = normal_force + forcex + forcey
            
            tau_pull_trough = total_normal_force/(np.pi * fastener.diameter * self.thickness)
            sigma_z_pull_trough = total_normal_force/((np.pi/4) * (fastener.outer_diameter**2-fastener.diameter**2))
            ms_shear_fastener = fastener.max_shear_load/shear_magnitude - 1
            ms_bearing_stress = self.sigma_y/sigma_br - 1
            ms_shear_stress_plate = (0.5*self.sigma_y)/tau_pull_trough - 1
            ms_sigma_z_pull_trough = self.sigma_y/sigma_z_pull_trough - 1
            minimum_ms_fastener = min(ms_sigma_z_pull_trough, ms_shear_fastener, ms_bearing_stress, ms_shear_stress_plate)
            if minimum_ms_fastener < minimum_ms:
                minimum_ms = minimum_ms_fastener
            if self.printing == True:
                print("----------------------------------------------------------------------")
                print("rivet nr:                                           ", n+1)
                print("bolt diameter [mm]                                  ", fastener.diameter*1000)
                
                print("ms shear load / ultimate shear load fastener:       ", ms_shear_fastener)
                print("ms pull trough load on the plate:                   ", ms_sigma_z_pull_trough)
                print("ms bearing stress / ultimate bearing stress sheet:  ", ms_bearing_stress)
                print("ms shear stress sheet / ultimate shear stress sheet:", ms_shear_stress_plate)
                print("minimum margin of safety:                           ", minimum_ms)
        return minimum_ms
    def optimize(self, required_safety_factor, diameter_ratio, tolerance, alpha):
        print("Optimizing bolt diameter...     ")
        print("Inital diameter:                ", self.fasteners[0].diameter)
        self.printing=False
        iterations = 1
        for i in range(0,1000):
            # Its just gradient descent
            factor = self.calculate_ms() + 1 - required_safety_factor
            if -tolerance < factor < tolerance:
                break
            for n, fastener in enumerate(self.fasteners):
                self.fasteners[n].diameter = self.fasteners[n].diameter - alpha*self.fasteners[n].diameter*factor
                self.fasteners[n].outer_diameter = diameter_ratio*self.fasteners[n].diameter
                self.fasteners[n].update()
            self.update_centroid()
            iterations += 1
        self.printing=True
        if iterations == 1001:
            print("Optimization did not converge, try a different alpha or tolerance")
        else:
            print(f"Optimzation converged after {iterations} iterations")
            print("final diameter:              ", self.fasteners[0].diameter)
"""
Launch case
"""
sigma_y_bolt = 400e6
sigma_y_plate = 400e6
plate_thickness = 0.4e-3
D = 1e-3
launch_joint = Joint(plate_thickness, sigma_y_plate)


# force and moment on the joint
force = np.array([749.9, 0, 749.9])
moment = np.array([0, 61.4, 0])

# setting up the coordinates of each bolt
n = 1
bL = 0.1
hL = 0.1
B = 0.2
local_joint_x = np.array([0, bL])
local_joint_y = np.array([0, hL])
joint_x = np.array(local_joint_x)

for i in range(1,n):
    joint_x = np.append(joint_x, local_joint_x + joint_x[-1] + B)
joint_y = np.array([0, hL])

# adding all the bolts to the joint
for x in joint_x:
    for y in joint_y:
        rivet = Fastener(D, 1.1*D, sigma_y_bolt, x, y)
        launch_joint.add_fastener(rivet)

launch_joint.add_load(force, moment)
launch_joint.update_centroid()
launch_joint.optimize(1.5, 1.2, 1e-5, 0.1)
launch_joint.calculate_ms()
