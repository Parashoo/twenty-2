import numpy as np

F_y=        #fill in or assign from earlier script
n_f=        #fill in or assign from earlier script
M_z=        #fill in or assign from earlier script
A=          #fill in or assign from earlier script              #Area of each fastener
r=          #fill in or assign from earlier script              #Radial distance of the center of each fastener to the fastener's c.g.

F_pi=F_y/n_f
F_pMz=M_z*A*r/sum(A*r**2)

F_yi=np.zeros(n_f)
for i in range(n_f):
    F_yi[i]=F_pMz[i]+F_pi

F_yi_ascending=np.sort(F_yi)
order_ascending=np.argsort(F_yi)+np.ones(n_f)
print("The fasteners with load in ascending order are:")
for i in range(n_f):
    print(i+1,":","fastener",order_ascending[i],"with total out of plane force:",F_yi_ascending[i])

