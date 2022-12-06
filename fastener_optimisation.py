"""
This code analyzes the minimum bolt size for fasteners in a joint. Here there are 2 
classes, one is the fastener, the other is the joint. A specific amount of fasteners
can be added to the joint at different positions in a normal carthesian coordinate system.
Aditonally, the force and moment can be added to the joint, and the margins of safety
can be computed. Also the bolt diameters can be optimized with a specific safety factor 
in mind.
"""
import numpy as np
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

# Thermal expansion coefficient for the wall the joint is attached to
wall_alpha = 22.68e-6

def induced_load_thermal(alpha_c,alpha_b,delta_T,E_b,A_sm,phi):
    F_delta_T=(alpha_c-alpha_b)*delta_T*E_b*A_sm*(1-phi)
    return(F_delta_T)

class Fastener():
    def __init__(self, diameter, outer_diamter, sigma_y, E, alpha, x, y):
        self.diameter = diameter
        self.outer_diameter = outer_diamter
        self.area = np.pi * diameter**2 / 4
        self.position = np.array([x, y, 0])
        self.sigma_y = sigma_y
        self.E = E
        self.alpha = alpha
        
        # Niu page 281
        # If you compare with the table of page 287 for diameter = 1/8 inch
        # and yield strength double of the shear strength they give, this is
        # a conservative estimate (531 lbs they vs 502 lbs formula after
        # conversions
        
        self.max_shear_load = 0.5*sigma_y * (np.pi*diameter**2)/4
    
    def update(self):
        self.area = np.pi * self.diameter**2 / 4
        self.max_shear_load = 0.5*self.sigma_y * np.pi*self.diameter**2/4

class Joint():
    def __init__(self, thickness, wall_thickness, sigma_y, E, alpha):
        self.fasteners = []
        self.force = np.zeros(3)
        self.moment = np.zeros(3)
        self.thickness = thickness
        self.wall_thickness = wall_thickness
        self.sigma_y = sigma_y
        self.E = E
        self.alpha = alpha
        # Printing is disabled during optimization
        self.printing = True
        
    # Function to add a fastener
    def add_fastener(self, fastener):
        self.fasteners.append(fastener)
    
    # Calculate important values for later computation of loading on each bolt
    def update_centroid(self):
        # sum(A*d**2) needed for shear loading for torque and normal loading due
        # to Mx and My
        self.J = 0   # Polar moment of area of the fastener group
        self.Ixx = 0
        self.Iyy = 0
        self.Ixy = 0
        
        # sum(A*d) vectorially for centroid determination
        self.total_FMOI = np.zeros(3)
        
        self.total_area = 0
        for fastener in self.fasteners:
            self.total_FMOI += fastener.position * fastener.area
            self.total_area += fastener.area
        self.centroid = self.total_FMOI * self.total_area**(-1)
        
        for fastener in self.fasteners:
            distance_centroid_vector = self.centroid - fastener.position
            distance_centroid = np.linalg.norm(distance_centroid_vector)
            self.J += distance_centroid**2 * fastener.area
            self.Ixx += distance_centroid_vector[1]**2 * fastener.area
            self.Iyy += distance_centroid_vector[0]**2 * fastener.area
            self.Ixy += distance_centroid_vector[0] * distance_centroid_vector[1] * fastener.area
        
    # Add force vector, moment vector and temperature difference to the joint
    def add_load(self, force, moment, deltaTs_joint, deltaTs_wall):
        self.force = force
        self.moment = moment
        self.deltaTs_joint = deltaTs_joint
        self.deltaTs_wall = deltaTs_wall
    
    def calculate_ms(self):
        torque = self.moment[2]
        minimum_ms = 1e99 #some arbitrary high number
        for i,fastener in enumerate(self.fasteners):
            """First calculating all the loads on every bolt"""
            # Calculating the shear load on the fastener due to shear
            shear_load_shear = self.force*(fastener.area/self.total_area)
            shear_load_shear[2] = 0 # only looking at the planar case
            
            # Calculating the shear load on the fastener due to torque
            centroid_distance = fastener.position - self.centroid 
            r = np.linalg.norm(centroid_distance)
            
            # Getting the torque shear force direction
            shear_torque_magnitude = torque * fastener.area * r / self.J
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
            
            # Only the magnitude is important
            shear_magnitude = np.linalg.norm(total_shear_load)
            
            # Calculating thermal additional loads
            phi = 0 # Worst case
            thermal_additional_loads = np.zeros(4)
            n = 0
            for delta_T in deltaTs_joint:
                thermal_load = (self.alpha-fastener.alpha)*delta_T*fastener.E*fastener.area*(1-phi)
                thermal_additional_loads[n] = thermal_load
                n+=1
            for delta_T in deltaTs_wall:
                thermal_load = (wall_alpha-fastener.alpha)*delta_T*fastener.E*fastener.area*(1-phi)
                thermal_additional_loads[n] = thermal_load
                n+=1
                
            # Getting the total shear magnitude by comparing the maximum of each
            # Combination of the thermal loading
            shear_magnitude = max(abs(thermal_additional_loads + shear_magnitude))
            
            # Normal force due to external normal force
            normal_force = self.force[2] * fastener.area/self.total_area
            
            # Shortening some variables to make reading the equation easier
            
            Mx = self.moment[0]
            # In the derivation in the book they take M_y opposite of the postive
            # internal direction
            My = -self.moment[1] 
            x = centroid_distance[0]
            y = centroid_distance[1]
            Ixx = self.Ixx
            Iyy = self.Iyy
            Ixy = self.Ixy
            denom = Ixx * Iyy - Ixy**2
            
            # Getting the force components due to assymetrical bending
            forcey = ((Mx * Iyy - My * Ixy) / denom ) * fastener.area * y
            forcex = ((My * Ixx - Mx * Ixy) / denom ) * fastener.area * x

            # Total force
            total_normal_force = abs(normal_force + forcex + forcey)
            
            """Now computing the stresses on the bolt..."""
            total_thickness = self.thickness+self.wall_thickness
            
            # Bending stress bolt
            sigma_z_bending_bolt = 16*shear_magnitude*total_thickness/(np.pi*fastener.diameter**3)
            
            # Shear stress bolt
            tau_bolt = shear_magnitude / (fastener.area)
            
            # Normal stress bolt
            sigma_z_bolt = total_normal_force/fastener.area
            
            # Total yield criterion bolt
            yield_criterion_bolt = np.sqrt((sigma_z_bolt + sigma_z_bending_bolt)**2 + 3*tau_bolt**2)
            
            """Now computing the stresses on the plate..."""
            # Bearing stress plate
            sigma_br_plate = shear_magnitude / (fastener.diameter*self.thickness)
            
            # Pull trough shear stress plate
            tau_pull_trough_plate = total_normal_force/(np.pi * fastener.diameter * self.thickness)
            
            # Pull trough normal stress plate and wall
            sigma_z_pull_trough = total_normal_force/((np.pi/4) * (fastener.outer_diameter**2-fastener.diameter**2))
            
            # Total yield criterion plate on the side without bearing stresses
            yield_criterion_plate_no_bearing = np.sqrt(sigma_z_pull_trough**2 + 3*tau_pull_trough_plate**2)
            
            # Total yield criterion plate on the side with bearing stresses
            yield_criterion_plate_bearing = np.sqrt(0.5*(sigma_z_pull_trough**2 + (abs(sigma_z_pull_trough)-abs(sigma_br_plate))**2+sigma_br_plate**2) + 3*tau_pull_trough_plate**2)  
            
            """Now the stresses on the wall..."""
            # Bearing stress wall
            sigma_br_wall = shear_magnitude / (fastener.diameter*self.wall_thickness)
            
            # Pull trough shear stress wall
            tau_pull_trough_wall = total_normal_force/(np.pi * fastener.diameter * self.wall_thickness)
            
            # Total yield criterion wall on the side without bearing stresses
            yield_criterion_wall_no_bearing = np.sqrt(sigma_z_pull_trough**2 + 3*tau_pull_trough_wall**2)
            
            # Total yield criterion wall on the side with bearing stresses
            yield_criterion_wall_bearing = np.sqrt(0.5*(sigma_z_pull_trough**2 + (abs(sigma_z_pull_trough)-abs(sigma_br_wall))**2+sigma_br_wall**2) + 3*tau_pull_trough_wall**2)  
            
            # Computing the margins of safety
            ms_bolt = fastener.sigma_y/yield_criterion_bolt - 1
            ms_plate_no_bearing = self.sigma_y / yield_criterion_plate_no_bearing - 1
            ms_plate_bearing = self.sigma_y / yield_criterion_plate_bearing - 1
            ms_wall_no_bearing = self.sigma_y / yield_criterion_wall_no_bearing - 1
            ms_wall_bearing = self.sigma_y / yield_criterion_wall_bearing - 1
            minimum_ms_fastener = min(ms_bolt, ms_plate_no_bearing, ms_plate_bearing, ms_wall_no_bearing, ms_wall_bearing)
            
            # Updating the minimum margin of safety if required
            if minimum_ms_fastener < minimum_ms:
                minimum_ms = minimum_ms_fastener
                
            # Print everything out if function is called externally
            if self.printing == True:
                print("----------------------------------------------------------------------")
                print("rivet nr:                                            ", i+1)
                print("rivet diameter [mm]                                  ", fastener.diameter*1000)
                
                print("ms bolt yield criterion                              ", ms_bolt)
                print("ms plate yield criterion bearing side:               ", ms_plate_bearing)
                print("ms plate yield criterion no bearing:                 ", ms_plate_no_bearing)
                print("ms wall yield criterion bearing side:                ", ms_plate_bearing)
                print("ms wall yield criterion no bearing:                  ", ms_plate_no_bearing)
                print("minimum margin of safety:                            ", minimum_ms_fastener)
        return minimum_ms
    def optimize(self, required_safety_factor, diameter_ratio, tolerance, alpha):
        print("Optimizing bolt diameter...     ")
        print("Inital diameter:                ", self.fasteners[0].diameter)
        
        # Disable printing to prevent console spam
        self.printing=False
        iterations = 1
        for i in range(0,1000):
            # Its simple gradient descent, where the additional factor is based
            # of both the ratio of stresses and the safety factor
            factor = self.calculate_ms() + 1 - required_safety_factor
            
            # Stop the loop if the factor is within the tolerance
            if -tolerance < factor < tolerance:
                break
            
            # Update every fastener diameter based on the factor
            for fastener in self.fasteners:
                fastener.diameter -= alpha*fastener.diameter*factor
                fastener.outer_diameter = diameter_ratio*fastener.diameter
                fastener.update()
            # Update the centroid of the joint
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

Material Wall:  Aluminum 6061 T6
Material Joint: Aluminum 6061 T6
Material rivet: Aluminum 6061 T6

(All the same due to extreme thermal stresses)
"""
# Material properties bolt
D = 1e-3
sigma_y_bolt = 276e6
E_bolt = 68e9
alpha_bolt = 22.68e-6

# Material properties plate
plate_thickness = 2e-3
sigma_y_plate = 276e6
E_plate = 68e9
alpha_plate = 22.68e-6

# Material properties wall
wall_thickness = 8e-3
sigma_y_wall = 276e6
E_wall = 68e9

launch_joint = Joint(plate_thickness, wall_thickness, sigma_y_plate, E_plate, alpha_plate)


# force and moment on the joint
force = np.array([749.9, 0, 749.9])
moment = np.array([0, 61.4, 0])

assembly_temp=288 # temperature during assembly
max_temp_joint=500
min_temp_joint=350
max_temp_wall=282
min_temp_wall=225

delta_T_max_joint=max_temp_joint-assembly_temp
delta_T_min_joint=min_temp_joint-assembly_temp
delta_T_max_wall=max_temp_wall-assembly_temp
delta_T_min_wall=min_temp_wall-assembly_temp

deltaTs_joint = [delta_T_min_joint, delta_T_max_joint]
deltaTs_wall  = [delta_T_min_wall, delta_T_max_wall]

# setting up the coordinates of each bolt
n = 1
bL = 0.06
hL = 0.03
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
        rivet = Fastener(D, 1.1*D, sigma_y_bolt, E_bolt, alpha_bolt, x, y)
        launch_joint.add_fastener(rivet)

launch_joint.add_load(force, moment, deltaTs_joint, deltaTs_wall)
launch_joint.update_centroid()
launch_joint.optimize(1.5, 1.2, 1e-5, 0.1)
launch_joint.calculate_ms()
