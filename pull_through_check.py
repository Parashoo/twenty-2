import numpy as np
import math

#4.8

F_y=1000        #fill in or assign from earlier script
n_f=4        #fill in or assign from earlier script
M_z=500        #fill in or assign from earlier script
r_i=0.005          #fill in or assign from earlier script      #Inner radius of the fasteners (assuming equal sized fasteners!)

#THIS HAS TO BE FILLED IN AS INPUT:
coords=np.array([[0.05,0.02],[0.05,-0.02],[-0.05,-0.02],[-0.05,0.02]])          #Coordinates of each fastener with origin at center of the lug. Amount of coords changes based on the number of fasteners (obviously)

A_fi=math.pi*r_i**2*np.ones(n_f)             #Inner Area of each fastener (same for all fasteners)
r=np.zeros(n_f)
for i in range(n_f):
    r[i]=((coords[i][0])**2+(coords[i][1])**2)**(1/2)
F_pi=F_y/n_f
F_pMz=M_z*A_fi*r/sum(A_fi*r**2)

F_yi=np.zeros(n_f)
for i in range(n_f):
    F_yi[i]=F_pMz[i]+F_pi

F_yi_ascending=np.sort(F_yi)
order_ascending=np.argsort(F_yi)+np.ones(n_f)
print("The fasteners with load in ascending order are:")
for i in range(n_f):
    print(i+1,":","fastener",order_ascending[i],"with total out of plane force:",F_yi_ascending[i])

#4.9

sigma_yield=267    #fill in or assign from earlier script      #Yield stress of fastener material im MPa
tau_yield=0      #fill in or assign from earlier script      #Shear yield stress of fastener material in MPa (fill in 0 if material property is unknown for that material)
D_fo=np.array([0.2,0.2,0.3,0.4])       #fill in or assign from earlier script      #Diameter of Fastener Head (in m)
t_2=0.001                            #Fill in or assign from earlier script         #Thickness of COMBINED plates (So total thickness of the 2 plates on top of eachother) (in m)

D_fi=np.ones(n_f)*r_i*2

sigma_yield=sigma_yield*1000000
if tau_yield==0:
    tau_yield=0.577*sigma_yield

A_fo=math.pi*(D_fo/2)**2        #Area of Fastener Head
A_force=A_fo-A_fi              #Area where F_y acts (part where the arrows are in figure 4.7 in reader)
A_shear=math.pi*D_fi*t_2

sigma_y=F_yi/A_force           #Tensile Stress in the 'Attached parts'
tau_y=F_yi/A_shear             #Shear Stress in the Fastener

Y=((sigma_y**2)+3*tau_y**2)**(1/2)        #Tension Yield Stress (from eq 4.8 in the reader)

#Check if they exceed yield stress and increase thickness of plates if they don't:
while max(abs(Y))>sigma_yield:
    t_2+=0.001
    A_shear=math.pi*D_fi*t_2
    tau_y=F_yi/A_shear
    Y=((sigma_y**2)+3*tau_y**2)**(1/2)
    print("Lug will fail, as Stress in one of the fasteners exceeds the yield stress, thickness t_2 has been increased to:",round(t_2,4),"; recalculating...")   
print("None of the fasteners will fail, as the maximum stress does not exceed the yield stress in any of them. Final thickness t_2=", round(t_2,4))
print("Reminder: t_2 is the total thickness of both plates in [m], not of a single plate")
print("The total stress in each fastener is:",Y)
